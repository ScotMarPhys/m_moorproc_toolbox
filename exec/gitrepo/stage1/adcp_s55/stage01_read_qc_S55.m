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

QC_NOT_EVALUATED = 0;
QC_GOOD = 1;
QC_UNKNOWN = 2;
QC_PROBABLY_BAD = 3;
QC_BAD = 4;
QC_CHANGED = 5;
QC_UNSAMPLED = 6;
QC_INTERPOLATED = 7;
QC_COMPASS_BAD = 8;
QC_MISSING = 9;

QC_vel = 0*double(Data.Average_VelEast);
QC_1D = 0*double(Data.Average_Pressure);

% NOTE: Do we want the ability to interactively flag and zoom on these early plots?
% E.g. It's a good opportunity to remove deployment and recovery based on
% pressure etc.  A bit of automatic flagging for demo in Fig. 1

%% figure settings
fs = 14
set(findall(gcf, '-property', 'FontSize'), 'FontUnits', 'points', 'FontSize', fs);

%% STAGE 1.  Nortek suggested quality control steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timemin = Data.Average_Time(1);
timemax = Data.Average_Time(end);

%% 1. Instrument orientation and pressure %%%%%%%%%%%%%%%%%%%%%%%%
y = Data.Average_Pressure; 
x = Data.Average_Time;

P_median = median(Data.Average_Pressure);
P_std = std(Data.Average_Pressure);
n_std = 1;

yl = (P_median+n_std*[-P_std P_std]);

% Check for deployment and recovery time
S = suggestTrimForShallowEdges_simple(y, min(yl),x);
disp(' ')
disp('***Pressure check***')
disp(sprintf('Found %d shallow values at start and %d at end.\nIndicates instrument was sinking/rising.', ...
                     S.nTrimStart, S.nTrimEnd));
disp('Flagged them as QC_BAD (4).')
disp(' ')

badind = ((S.maskStart+S.maskEnd)==1);
QC_vel(badind,:) = QC_BAD; QC_1D(badind) = QC_BAD;

figure(1);
subplot(3,1,1)
title([filename ' pressure'],'Interpreter','none');
hold on; grid on;

% Plot data
h.data = plot(Data.Average_Time,Data.Average_Pressure,'k.');
yline(P_median,'r')
[h, nAbove, nBelow] = markClippedPoints_pres(h, y, yl, x);

datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');
xlim([timemin timemax]);
ylim(yl);
ylabel('Pressure (db)');
set(gca,'YDir','reverse');

% Create a neat text box on the figure (normalized figure coordinates)
txt = {['ylim = median \pm ' sprintf('%d·std',n_std)], ...
       [sprintf('Found %d shallow values at start and %d at end. Flagged them as bad (%d). ', ...
                     S.nTrimStart, S.nTrimEnd,QC_BAD)]};

% Position: upper-right of axes (tweak if needed)
ax = gca;
axPos = get(ax, 'Position'); % in normalized figure units
% Compute a small textbox inside the axes at top-right
bx = axPos(1) + 0.33*axPos(3);
by = axPos(2) + 0.7*axPos(4);
bw = 0.32*axPos(3);
bh = 0.28*axPos(4);

hBox = annotation('textbox', [bx by bw bh], 'String', txt, ...
    'FitBoxToText', 'on', 'EdgeColor', 'none', ...
    'Interpreter', 'tex', 'HorizontalAlignment', 'left');

% Optional: give the box a subtle background for readability
% set(hBox, 'BackgroundColor', [1 1 1 0.85]); % if MATLAB supports RGBA; otherwise use [1 1 1]

%% 2. Tilt / pitch %%%%%%%%%%%%%%%%%%%%%%%%
y = Data.Average_Pitch;
y(QC_1D==QC_BAD) = NaN;
x = Data.Average_Time;

y_trim = Data.Average_Pitch;
y_trim(QC_1D==0) = NaN;

yl =  [-40 40];

subplot(3,1,2) 
title([filename ' pitch and roll'],'Interpreter','none');
hold on; grid on;

% plot data
h.data = plot(x,y,'.','DisplayName', 'pitch raw bottom');
plot(x,y_trim,'.c','DisplayName', 'pitch QC_BAD sinking/rising');

plot(Data.Average_Time(QC_1D==0),Data.Average_Roll(QC_1D==0),'k','DisplayName', 'roll raw bottom','LineWidth',0.5);

% [h, nAbove, nBelow] = markClippedPoints(h, y,yl, x);

datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');

% Suggested quality thresholds
hl1 = yline(0,'color','k'); hl1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hl1 = yline([-10 10],'color','g');  arrayfun(@(h) set(h.Annotation.LegendInformation, 'IconDisplayStyle', 'off'), hl1);
hl1 = yline([-30 30],'color','r');  arrayfun(@(h) set(h.Annotation.LegendInformation, 'IconDisplayStyle', 'off'), hl1);

xlim([timemin timemax]);
ylim(yl);

ys_data = [ 10, -10, 30, -30 ];
ys_labels = { '10 < 30: Post processing possible', ...
           '-10 > -30: Post processing possible', ...
           '>30: Profiles probably bad', ...
           '<-30: Profiles probably bad' };
ys_colors =  {'g','g','r','r'};

annotateYLabels(ys_data, ys_labels, ys_colors);

ylabel('Pitch (^o)');

% Automatically flag any bad timesteps: excessive pitch 
badind = find(y > 30 | y < -30);

if ~isempty(badind)
    plot(y(badind),y(badind),'or', 'DisplayName', 'QC_BAD pitch bottom');
    QC_vel(badind,:) = QC_BAD;
    QC_1D(badind) = QC_BAD;
end
disp('***Pitch check***')
disp([num2str(length(badind)) ' timesteps flagged bad (',num2str(QC_BAD),') due to excessive pitch (>abs(30))']);
disp(' ')
legend('Interpreter','none')

% Automatically flag any bad timesteps with postprocessing possible: high
% pitch
badind = find( (y > 10 & y < 30) | (y < -10 & y > -30) );

if ~isempty(badind)
    plot(y(badind),y(badind),'or', 'DisplayName', 'QC_PROBABLY_BAD pitch bottom');
    QC_vel(badind,:) = QC_PROBABLY_BAD;
    QC_1D(badind) = QC_PROBABLY_BAD;
end
disp([num2str(length(badind)) ' timesteps flagged probably bad (',num2str(QC_PROBABLY_BAD),') due to pitch between +(-) 10 and +(-) 30']);
disp('Post processing possible.');
disp(' ')
legend('Interpreter','none','Location','southeast')


%% 3. heading %%%%%%%%%%%%%%%%%%%%%%%%
yl = [-0 360];
y = Data.Average_Heading;
x = Data.Average_Time; 

y(QC_1D==4) = NaN;
y_trim =  Data.Average_Heading;
y_trim(QC_1D==0) = NaN;

subplot(3,1,3)
title([filename ' heading'],'Interpreter','none');
hold on; grid on;

% Plot data
plot(x,y,'k.','DisplayName', 'heading raw bottom');
plot(x,y_trim,'.c','DisplayName', 'QC_BAD sinking/rising/pitch');

datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');
xlim([timemin timemax]);

ylim(yl);
ylabel('Heading (^o)');
[h, nAbove, nBelow] = markClippedPoints(h, y,yl,x);

legend('Interpreter','none','Location','southeast')

%% add remarks
% place at top-right of figure (normalized units)
str = 'Pitch and heading should remain fairly constant during deployment';
an = annotation('textbox', [0 0.02 1 0.06], ...   % [x y w h] in normalized figure units
                'String', str, ...
                'FitBoxToText', 'on', ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', ...
                'Interpreter', 'none', ...
                'Color', [0 0.5 0]);  % change color if desired


% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',[figureoutdir filename '_f1_pressure_pitch_heading_QC']);


%% fuctions

% looks for deployment and recovery period in pressure
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

%marks if data is outside of y-axes lim for pressure (positive downward)
function [h, nAbove, nBelow] = markClippedPoints_pres(h, y, yl, x)
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

%marks if data is outside of y-axes lim (positive upward)
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
    h.above = scatter(x(deeper), y_clipped(deeper), 50, 'r', '^', 'filled');
else
    h.above = gobjects(0);
end
if nBelow > 0
    h.below = scatter(x(shallower), y_clipped(shallower), 50, 'r', 'v', 'filled');
else
    h.below = gobjects(0);
end

% Ensure y-limits are set as requested
ylim(yl);

% Optional annotations: counts at center of clipped x-range
if nAbove > 0
    xc = mean(x(deeper));
    text(xc, ymax, sprintf('%d below', nAbove), ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', 'r');
end
if nBelow > 0
    xc = mean(x(shallower));
    text(xc, ymin, sprintf('%d above ', nBelow), ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Color', 'r');
end

hold off
end

function annotateYLabels(ys_data, labels, colors)
% annotateYLabels Place staggered labels for horizontal reference lines.
%   annotateYLabels(ys_data, labels) places labels at the right edge (x=0.98
%   in normalized units) for the y positions given in data units in ys_data.
%   labels is a cell array of strings the same length as ys_data.
%
%   annotateYLabels(ys_data, labels, colors) lets you provide colors as a
%   cell array (or char array) the same length as ys_data.
%
%   Labels are placed in 'Units','normalized' and staggered to avoid overlap.
%
%   Example:
%     ys = [12, -12, 32, -32];
%     labs = {'10 < 30: Post processing possible', ...
%             '-10 > -30: Post processing possible', ...
%             '>30: Profiles probably bad', ...
%             '<-30: Profiles probably bad'};
%     annotateYLabels(ys, labs, {'g','g','r','r'});

if nargin < 2
    error('Need ys_data and labels.');
end
if ~iscell(labels)
    labels = cellstr(labels);
end
n = numel(ys_data);
if numel(labels) ~= n
    error('ys_data and labels must have the same length.');
end
if nargin < 3 || isempty(colors)
    colors = repmat({'k'}, 1, n);
end
% normalize colors to cell array
if ischar(colors) || (isstring(colors) && isscalar(colors))
    colors = repmat({char(colors)}, 1, n);
elseif isstring(colors)
    colors = cellstr(colors);
end

ax = gca;
% compute normalized coordinates
yn = (ys_data - ax.YLim(1)) / diff(ax.YLim);

% decide grouping: above middle and below middle, stagger independently
mid = 0.5;
upper_idx = find(yn >= mid);
lower_idx = find(yn <  mid);

stagger = 0.05;   % normalized units between successive labels (adjust if needed)

% sort within groups to keep order stable (closest to edge first)
[~, s1] = sort(yn(upper_idx), 'ascend');  upper_idx = upper_idx(s1);
[~, s2] = sort(yn(lower_idx), 'descend'); lower_idx = lower_idx(s2);

% apply staggering: move upper labels upward, lower labels downward
for k = 1:numel(upper_idx)
    idx = upper_idx(k);
    yn(idx) = yn(idx) + stagger;
end
for k = 1:numel(lower_idx)
    idx = lower_idx(k);
    yn(idx) = yn(idx) - 0.02 -stagger;
end

xn = 0.028;   % x position in normalized units (right edge)
for k = 1:n
    text(xn, yn(k), labels{k}, ...
         'Units', 'normalized', ...
         'HorizontalAlignment', 'left', ...
         'VerticalAlignment', 'middle', ...
         'Color', colors{k}, ...
         'Interpreter', 'none');
end
end


