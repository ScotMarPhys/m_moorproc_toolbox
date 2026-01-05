%%  Code for the quality control of Signature 55 ADCP

% K Burmeister, S Jones 12/2025


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PARAMETER PRAEMBLE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% Need to be updated for:
% - directories and file names
% - add additional deployment periods (add moorX = 'filename')
% - turn off/on checkplots
% - edit start/end of total time period
% - data version
% - set depth of shallowest instrument (idepth)
% - preamble for despiking
% - if you add a deployment period, you need to edit steps 1-3
%---------
close all
clear

addpath(genpath('../functions'));

% TO UPDATE <------------------------
% in- and output directories
pc_name = getenv('COMPUTERNAME');
if strcmp(pc_name,'SA07KB-3JN9YY2');
    basedir = 'C:\Users\sa07kb\Projects\Moor_Data_Proc\';
    dataindir = [basedir,'moor_examples\osnap\data\moor\raw\jc238\s55\'];
    pathgit = [basedir 'm_moorproc_toolbox\'];
    figureoutdir = [basedir,'RHADCP\figures\'];
    filename = 'S200044A008_RHADCP_2020';
    addpath(genpath('D:\Work_computer_sync\MATLAB_functions')); % General functions
elseif strcmp(pc_name,'SA01SJ-G9WC2J3')
    dataindir = 'E:\OSNAP\RHADCP\DY181\S200044A012_RHAD2_JC238\conversion2\';
    pathgit = 'D:\Work_computer_sync\OSNAP_postdoc\Python\m_moorproc_toolbox\';    
    figureoutdir = ['D:\Work_computer_sync\OSNAP_postdoc\Mooring\RHADCP\plots\'];
    addpath(genpath('D:\Work_computer_sync\MATLAB_functions')); % General functions
    filename = 'S200044A012_RHAD2_JC238';
else
    error('Please add your path above')
end


%% Load
load([dataindir filename '.mat']);


%% Build post-processing flag array

% 0: QC_NOT_EVALUATED
% 1: QC_GOOD
% 2: QC_UNKNOWN
% 3: QC_PROBABLY_BAD
% 4: QC_BAD
% 5: QC_CHANGED
% 6: QC_UNSAMPLED
% 7: QC_INTERPOLATED
% 8: ?
% 9: QC_MISSING

QC_vel = 0*double(Data.Average_VelEast);
QC_1D = 0*double(Data.Average_Pressure);

% NOTE: Do we want the ability to interactively flag and zoom on these early plots?
% E.g. It's a good opportunity to remove deployment and recovery based on
% pressure etc.  A bit of automatic flagging for demo in Fig. 1

%% figure settings
set(findall(gcf, '-property', 'FontSize'), 'FontUnits', 'points', 'FontSize', 14);

%% STAGE 1.  Nortek suggested quality control steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timemin = Data.Average_Time(1);
timemax = Data.Average_Time(end);

%% 1. Pressure %%%%%%%%%%%%%%%%%%%%%%%%
y = Data.Average_Pressure; 
x = Data.Average_Time;

P_median = median(Data.Average_Pressure);
P_std = std(Data.Average_Pressure);
n_std = 3;

yl = (P_median+n_std*[-P_std P_std]);

% Check for deployment and recovery time
S = suggestTrimForShallowEdges_simple(y, min(yl),x);

figure(1);
title([filename ' pressure'],'Interpreter','none');
hold on; grid on;


% Plot data
h.data = plot(Data.Average_Time,Data.Average_Pressure,'k.');
yline(P_median,'r')


[h, nAbove, nBelow] = markClippedPoints(h, y, yl, x);

datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');
xlim([timemin timemax]);
ylim(yl);
xlabel('Date');
ylabel('Pressure (db)');
set(gca,'YDir','reverse');

% Create a neat text box on the figure (normalized figure coordinates)
txt = {['ylim = median \pm ' sprintf('%d·std',n_std)], ...
       ['median = ' sprintf('%4.0f (db)', P_median)], ...
       ['std    = ' sprintf('%.3g (db)', P_std)],...
       ['See prombt to set QC']};

% Position: upper-right of axes (tweak if needed)
ax = gca;
axPos = get(ax, 'Position'); % in normalized figure units
% Compute a small textbox inside the axes at top-right
bx = axPos(1) + 0.66*axPos(3);
by = axPos(2) + 0.66*axPos(4);
bw = 0.32*axPos(3);
bh = 0.28*axPos(4);

hBox = annotation('textbox', [bx by bw bh], 'String', txt, ...
    'FitBoxToText', 'on', 'EdgeColor', 'none', ...
    'Interpreter', 'tex', 'HorizontalAlignment', 'left');

% Optional: give the box a subtle background for readability
set(hBox, 'BackgroundColor', [1 1 1 0.85]); % if MATLAB supports RGBA; otherwise use [1 1 1]

% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',[figureoutdir filename '_f3_pressure_QC']);

while true
    prompt = sprintf('Found %d shallow values at start and %d at end.\nFlag them as QC_BAD? Y/N [Y]: ', ...
                     S.nTrimStart, S.nTrimEnd);
    resp = input(prompt, 's');
    if isempty(resp), resp = 'Y'; end
    resp = upper(strtrim(resp));
    if strcmp(resp,'Y')
        doFlag = true;
        break
    elseif strcmp(resp,'N')
        doFlag = false;
        break
    else
        fprintf('Invalid input — try again.\n\n');
    end
end

if doFlag
   badind = ((S.maskStart+S.maskEnd)==1);
   QC_vel(badind,:) = 4;
   QC_1D(badind) = 4;
end

%% 2. Tilt / pitch %%%%%%%%%%%%%%%%%%%%%%%%

figure(2); 

title([filename ' pitch'],'Interpreter','none');
hold on; grid on;
y_trim = Data.Average_Pitch(QC_1D==4);
t_trim = Data.Average_Time(QC_1D==4);

% plot data
h.data = plot(Data.Average_Time,Data.Average_Pitch,'.');
plot(t_trim,y_trim,'oc');

datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');

% Suggested quality thresholds
line([timemin timemax],[0 0],'color','k');
line([timemin timemax],[10 10],'color','g');
line([timemin timemax],[-10 -10],'color','g');
line([timemin timemax],[30 30],'color','r');
line([timemin timemax],[-30 -30],'color','r');

text(timemin+50,12,'10 < 30: Post processing possible','color','g');
text(timemin+50,-12,'-10 > -30: Post processing possible','color','g');
text(timemin+50,32,'>30: Profiles probably bad','color','r');
text(timemin+50,-32,'<-30: Profiles probably bad','color','r');


xlim([timemin timemax]);
ylim([-40 40]);

xlabel('Date');
ylabel('Pitch (^o)');

% FLAG EXAMPLE: Automatically flag any bad timesteps
badind = find(Data.Average_Pitch > 30 | Data.Average_Pitch < -30);
plot(Data.Average_Time(badind),Data.Average_Pitch(badind),'or');
QC_vel(badind,:) = 4;
disp([num2str(length(badind)) ' timesteps flagged bad (4) due to excessive pitch']);

[h, nAbove, nBelow] = markClippedPoints(h, Data.Average_Pitch, [-40 40], Data.Average_Time);

% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',[figureoutdir filename '_f1_pitch_QC']);

%% 3. heading %%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
title([filename ' heading'],'Interpreter','none');
hold on; grid on;

% Plot data
plot(Data.Average_Time,Data.Average_Heading,'g.');
plot(Data.Average_Time,Data.Average_Heading,'g');

datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');
xlim([timemin timemax]);
ylim([-0 360]);
xlabel('Date');
ylabel('Heading (^o)');

% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',[figureoutdir filename '_f2_heading_QC']);



%% fuctions

function S = suggestTrimForShallowEdges_simple(y, ymin,x)
% Suggest trim indices for leading/trailing consecutive y < ymin.
% S.startIdx = first index to keep
% S.endIdx   = last  index to keep
% S.nTrimStart, S.nTrimEnd, S.keepMask

y = y(:);
n = numel(y);
isShallow = y < ymin;        % NaN -> false

% Find first and last keepable samples (non-shallow and not NaN)
firstKeep = find(~isShallow & ~isnan(y), 1, 'first');
lastKeep  = find(~isShallow & ~isnan(y), 1, 'last');

if isempty(firstKeep)      % no non-shallow non-NaN values
    S.startIdx = n+1;
    S.endIdx   = 0;
    S.nTrimStart = nnz(isShallow);
    S.nTrimEnd   = 0;
    S.maskStart = isShallow;
    S.maskEnd = false(n,1);
    S.keepMask = false(n,1);
    return
end

% Leading trim: contiguous shallow run from index 1 up to first non-shallow
leadRunLen = find(~isShallow(1:firstKeep-1), 1, 'first'); %#ok
if isempty(leadRunLen)
    nTrimStart = firstKeep-1;        % all samples before firstKeep are shallow
else
    nTrimStart = leadRunLen-1;       % stops before first non-shallow inside that prefix
end

% Trailing trim: contiguous shallow run from lastKeep+1 to end
if lastKeep < n
    tailPrefix = isShallow(lastKeep+1:end);
    tailNonSh = find(~tailPrefix, 1, 'first');
    if isempty(tailNonSh)
        nTrimEnd = numel(tailPrefix);   % all after lastKeep are shallow
    else
        nTrimEnd = tailNonSh-1;         % stops at first non-shallow after lastKeep
    end
else
    nTrimEnd = 0;
end

startIdx = firstKeep - nTrimStart;
endIdx   = lastKeep + nTrimEnd;

% Build masks
maskStart = false(n,1); if nTrimStart>0, maskStart(1:firstKeep-1) = true; maskStart(~isShallow) = false; end
maskEnd   = false(n,1); if nTrimEnd>0, maskEnd(lastKeep+1:end) = true; maskEnd(~isShallow) = false; end
keepMask = false(n,1);
if startIdx <= endIdx
    keepMask(startIdx:endIdx) = true;
end

S = struct('startIdx', startIdx, 'endIdx', endIdx, ...
           'nTrimStart', nTrimStart, 'nTrimEnd', nTrimEnd, ...
           'maskStart', maskStart, 'maskEnd', maskEnd, 'keepMask', keepMask);
end

function [h, nAbove, nBelow] = markClippedPoints(h, y, yl, x)
% markClippedPoints Mark points of y outside y-limits yl.
%   [h, nAbove, nBelow] = markClippedPoints(y, yl)
%   [h, nAbove, nBelow] = markClippedPoints(y, yl, x)
%
% Inputs
%   y  - vector of data values
%   yl - two-element vector [ymin ymax]
%   x  - (optional) x positions (same length as y). Default 1:numel(y)
%
% Outputs
%   h      - struct with handles: h.data, h.above, h.below
%   nAbove - number of points above yl(2)
%   nBelow - number of points below yl(1)

if nargin < 3 || isempty(x)
    x = 1:numel(y);
end
y = y(:); x = x(:);
assert(numel(x) == numel(y), 'x and y must have same length');
assert(numel(yl) == 2, 'yl must be [ymin ymax]');

ymin = yl(1); ymax = yl(2);

% Find outside points
deeper = y > ymax;
shallower = y < ymin;
nAbove = nnz(deeper);
nBelow = nnz(shallower);

% Clip positions at edge for marker placement
y_clipped = y;
y_clipped(deeper) = ymax;
y_clipped(shallower) = ymin;

% Markers: triangles pointing up for above, down for below
if nAbove > 0
    h.above = scatter(x(deeper), y_clipped(deeper), 50, 'r', 'v', 'filled');
else
    h.above = gobjects(0);
end
if nBelow > 0
    h.below = scatter(x(shallower), y_clipped(shallower), 50, 'r', '^', 'filled');
else
    h.below = gobjects(0);
end

% Ensure y-limits are set as requested
ylim(yl);

% Optional annotations: counts at center of clipped x-range
if nAbove > 0
    xc = mean(x(deeper));
    text(xc, ymax, sprintf('%d deeper', nAbove), ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', 'r');
end
if nBelow > 0
    xc = mean(x(shallower));
    text(xc, ymin, sprintf('%d shallower ', nBelow), ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Color', 'r');
end

hold off
end