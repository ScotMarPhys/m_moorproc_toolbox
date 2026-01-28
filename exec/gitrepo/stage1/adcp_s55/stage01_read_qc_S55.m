% stage01_read_qc_S55.mQC_vel
% Read and quality control Signature 55 ADCP data
% raw ADCP data in .mat format as exported by Signature Viewer
%
% required inputs: moor - mooring name e.g. 'rhadcp_01_2020'
%                  
% optional inputs:     dataindir = varargin{1};
%                      filename = varargin{2}
%                      infofile = varargin{3};
%                      logfile = varargin{4};
%                      outdir = varargin{5}; #output dir (fig, data, logs)
% functions called:    rodbload
%
% 

function stage01_read_qc_S55(moor, varargin)

if nargin==0
    help stage01_read_qc_S55
    return
end

if nargin==1
    global MOORPROC_G
    operator = MOORPROC_G.operator;
    pd = moor_inoutpaths('adcp_S55',moor);
    dataindir = pd.rawpath;
    infofile = pd.infofile;
    logfile = pd.stage1log;
    outdir = pd.stage1path;
    ouput_form = pd.stage1form;
else
    operator = getenv('COMPUTERNAME');    
    dataindir = varargin{1};
    filename = varargin{2};
    infofile = varargin{3};
    logfile = varargin{4};
    outdir = varargin{5};
    ouput_form = [moor '_%d.nc'];
end

if exist(infofile,"file")
    % ----- read infofile / open logfile  ------------------------------------
    infovar = 'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth'; 
    [id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd]  =  rodbload(infofile,infovar);
else
    disp('No info file given. Parameters set to NaN')
    [id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd] = deal(NaN);
end

if ~exist(outdir,'dir')
    mkdir(outdir)
end

fidlog   = fopen(logfile,'a');
fprintf(fidlog, '\n==== START ENTRY  =====\n');
fprintf(fidlog,'Read and quality control Signature 55 ADCP data. \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));
fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n\n\n',lon);

bg = datenum(datetime([s_d.' , s_t.' , 0])); %start
ed = datenum(datetime([e_d.' , e_t.' , 0])); %end

if isnan(id)
   fprintf('No serial number given, use default 200044');
   serial_nums=200044;
else
    vec=find((id>=319) & (id <=328)); % Possible ADCP codes - taken from IMP moorings package
    serial_nums=sn(vec);
end

%% Load
for i = 1:length(serial_nums)
    fprintf('Processing sn %d',serial_nums(i));
    if ~exist('filename',"var")
        filename = sprintf('%d_data',serial_nums(i));
    end
infile=fullfile(dataindir,[filename,'.mat']);
load(infile);
Data = cell2struct(cellfun(@double,struct2cell(Data),'uni',false),fieldnames(Data),1);

fprintf(fidlog,'infile : %s\n',infile);
fprintf(fidlog,'ADCP serial number  : %d\n',serial_nums(i));

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

%% Checking few configuration which should be the same each time step
if strcmp(getenv('COMPUTERNAME'),'SA07KB-3JN9YY2'); % Numunique in V2025a onwards, temporarily skip check
    if serial_nums(i)~=Config.SerialNo
        error('Entered Serial Number does not match Config-file')
    end
    prombt = ' is not equal for all time steps. Please check ';
    if numunique(Data.Average_NBeams,"rows") ~= 1
        text = ['Number of beams',prombt,'Data.Average_NBeams'];
        disp(text); fprintf(fidlog,[text,'\n']);
    elseif numunique(Data.Average_NCells,"rows") ~= 1
        text = ['Number of cells',prombt,'Data.Average_NCells'];
        disp(text); fprintf(fidlog,[text,'\n']);
    elseif numunique(Data.Average_BeamToChannelMapping,"rows") ~= 1
        text = ['Beam to channel mapping',prombt,'Data.Average_BeamToChannelMapping'];
        disp(text); fprintf(fidlog,[text,'\n']);
    elseif numunique(Data.Average_Error,'rows')~= 1 | Data.Average_Error(1)~=0
        text = ['Errors occured during measuring. Please check Data.Average_Error'];
        disp(text); fprintf(fidlog,[text,'\n']);
    elseif numunique(Data.Average_Soundspeed,'rows')~= 1
        text = ['Sound speed not constant. Please check Data.Average_Soundspeed'];
        disp(text); fprintf(fidlog,[text,'\n']);
    end
else

end

%% figure settings
fs = 14;
set(findall(gcf, '-property', 'FontSize'), 'FontUnits', 'points', 'FontSize', fs);

%% STAGE 1.  Nortek suggested quality control steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timemin = Data.Average_Time(1);
timemax = Data.Average_Time(end);

%% 1. Instrument orientation and pressure %%%%%%%%%%%%%%%%%%%%%%%%
% pressure
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

fprintf(fidlog,'***Pressure check***\n');
fprintf(fidlog,sprintf('Found %d shallow values at start and %d at end.\nIndicates instrument was sinking/rising.\n', ...
                     S.nTrimStart, S.nTrimEnd));
fprintf(fidlog,'Flagged them as QC_BAD (4).\n');

% if bg
ds_vec = datevec(x([S.startIdx,S.endIdx]));
bg_suggest = datenum([ds_vec(1,1:end-1),0]);
ed_suggest = datenum([ds_vec(2,1:end-1),0]);
if bg_suggest~=bg
    prombt = ['***',newline,'Suggest to set mooring start date and time in ',... 
            moor,'info.dat to ',datestr(bg_suggest),newline,'***',newline,newline];
    disp(prombt);
    fprintf(fidlog,prombt);
elseif ed_suggest~=ed
        prombt = ['***',newline,'Suggest to set mooring end date and time in ',... 
            moor,'info.dat to ',datestr(ed_suggest),newline,'***',newline,newline];
    disp(prombt);
    fprintf(fidlog,prombt);
end

badind = (S.keepMask==0);
QC_vel(badind,:) = QC_BAD; QC_1D(badind) = QC_BAD;

figure(1);clf
subplot(3,1,1)
title([filename ' pressure'],'Interpreter','none');
hold on; grid on;

% Plot data
h.data = plot(x,y,'k.');
yline(P_median,'r')
[h, nAbove, nBelow] = markClippedPoints_pres(h, y, yl, x);

datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');
xlim([timemin timemax]);
ylim(yl);
ylabel('Pressure (db)');
set(gca,'YDir','reverse');

% --- RIGHT AXIS: Battery ---
yyaxis right
hold on
h_bat = plot(Data.Average_Time, Data.Average_Battery, 'b-', 'LineWidth', 1); 
ylabel('Battery (V)');
ylim([0,20])
set(gca, 'YColor', 'b'); 
hold off

% Create a neat text box on the figure (normalized figure coordinates)
txt = {['ylim = median \pm ' sprintf('%d·std',n_std)], ...
       [sprintf('Found %d shallow values at start and %d at end.\nFlagged them as bad (%d). ', ...
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

% Tilt / pitch %%%%%%%%%%%%%%%%%%%%%%%%
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
fprintf(fidlog,'***Pitch check***\n');
fprintf(fidlog,[num2str(length(badind)) ' timesteps flagged bad (',...
    num2str(QC_BAD),') due to excessive pitch (>abs(30))\n']);

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
fprintf(fidlog,[num2str(length(badind)) ' timesteps flagged probably bad (',...
    num2str(QC_PROBABLY_BAD),') due to pitch between +(-)10 and +(-)30.\n']);
fprintf(fidlog,'Post processing possible.\n\n');

legend('Interpreter','none','Location','southeast')


% heading %%%%%%%%%%%%%%%%%%%%%%%%
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

% add remarks
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
print('-dpng',fullfile(outdir,[filename,'_f1_pressure_pitch_heading_QC.png']));


%% 2. Beam amplitudes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nCells = 1:Data.Average_NCells(1);  %number of cells
bl_dist = Config.Average_BlankingDistance; %blanking distance (m)
Cell_Size = Config.Average_CellSize; %cell size (m)
Dist2Instr_CellMidpoint = bl_dist+nCells.*Cell_Size; 

if ~isnan(wd)
    Nominal_CellDepth = wd-Dist2Instr_CellMidpoint;
else
    wd = 1080;
    while true
        prompt = sprintf('Please enter nominal depth of instrument in m (default %dm): ', wd);
        val = input(prompt);
        if isempty(val)                     % user accepted default
            wd = wd;
            break
        end
        if isnumeric(val) && isscalar(val) && isfinite(val) && val > 0
            wd = val;
            break
        end
        fprintf('Invalid entry. Enter a positive number or press Return for default.\n');
    end
    Nominal_CellDepth = wd-Dist2Instr_CellMidpoint;
end

fprintf(fidlog,'***Beams, surface bins, sidelobe check***\n');
fprintf(fidlog,sprintf('Nominal depth of instrument set as %d dbar.\n',wd));

y=Nominal_CellDepth;
x = Data.Average_Time;
cb_lim =[0,100];

% define surface bins
Amp1=Data.Average_AmpBeam1; Amp1(QC_1D==4,:)=NaN;
Amp2=Data.Average_AmpBeam2; Amp2(QC_1D==4,:)=NaN;
Amp3=Data.Average_AmpBeam3; Amp3(QC_1D==4,:)=NaN;
Cor1=Data.Average_CorBeam1; Cor1(QC_1D==4,:)=NaN;
Cor2=Data.Average_CorBeam2; Cor2(QC_1D==4,:)=NaN;
Cor3=Data.Average_CorBeam3; Cor3(QC_1D==4,:)=NaN;

% mask any Cor<50 as bad
mask_corr = (Cor1 < 50) | (Cor2 < 50) | (Cor3 < 50);
QC_vel(mask_corr)=QC_BAD;

U = Data.Average_VelEast;U(QC_1D==4,:)=NaN;
V = Data.Average_VelNorth;V(QC_1D==4,:)=NaN;
W = Data.Average_VelUp;W(QC_1D==4,:)=NaN;

% surface bin detection
Amp1_pro = nanmean(Amp1);SB1 = find(islocalmin(Amp1_pro),1,'last');
Amp2_pro = nanmean(Amp2);SB2 = find(islocalmin(Amp2_pro),1,'last');
Amp3_pro = nanmean(Amp3);SB3 = find(islocalmin(Amp3_pro),1,'last');
Cor1_pro = nanmean(Cor1);CB1 = find(Cor1_pro>=50,1,'last');
Cor2_pro = nanmean(Cor2);CB2 = find(Cor2_pro>=50,1,'last');
Cor3_pro = nanmean(Cor3);CB3 = find(Cor3_pro>=50,1,'last');
SB = min([SB1,SB2,SB3]);
CB = min([CB1,CB2,CB3]);

% side lobe interference
SM1 = find(islocalmax(Amp1_pro),1,'last');
SM2 = find(islocalmax(Amp2_pro),1,'last');
SM3 = find(islocalmax(Amp3_pro),1,'last');
SM = min([SM1,SM2,SM3]);
A(1) = gsw_z_from_p(nanmean(Data.Average_Pressure(QC_1D==0)),lat); % range based on pressure sensor
A(2) = Dist2Instr_CellMidpoint(SM); % range based on nominal cell distance
%slant angle 20 degree plus abs(pitch)
pitch = Data.Average_Pitch; pitch(QC_1D==4)=NaN;
theta = abs(pitch)+20;
% R contains buffer for upper range of cell
R =max(A)*cosd(theta)-Cell_Size; 
% find(Dist2Instr_CellMidpoint <= R(i), 1, 'last') for each time step i
R(QC_1D==4)=wd;
idx_valid = arrayfun(@(r) find(Dist2Instr_CellMidpoint <= r, 1, 'last'),...
    R);
idx_valid(QC_1D==4)=0;
R(QC_1D==4)=NaN;


%%%%%%%%%%%%%%%%%%%%%%%
f2 = figure(2);clf
ax(1) = subplot(2,3,1);hold on; imagesc(x, y, Amp1');
ax(2) = subplot(2,3,2);hold on; imagesc(x, y, Amp2');
ax(3) = subplot(2,3,3);hold on; imagesc(x, y, Amp3');

ax(4) = subplot(2,3,4);hold on; imagesc(x, y, Cor1');
ax(5) = subplot(2,3,5);hold on; imagesc(x, y, Cor2');
ax(6) = subplot(2,3,6);hold on; imagesc(x, y, Cor3');

for k = 1:numel(ax)
    axis(ax(k),'xy','ij');
    ylim(ax(k),[0,1080])
    ylabel(ax(k),'Nominal cell depth (m)')
    datetick(ax(k),'x','mmm-yyyy','keepticks','keeplimits');
    hLine = plot(ax(k),x(QC_1D==0), y(idx_valid(QC_1D==0)), 'r', 'LineWidth', 0.1);
    hLine.DisplayName = 'Sidelobe interference';
    xlim(ax(k),[min(x),max(x)])
end

for k=1:3
    axes(ax(k));                   % make axis current
    title(sprintf('Amplitude Beam %d',k))
    caxis(ax(k), cb_lim); 
    colorbar(ax(k))
end

lgd = legend(ax(1),'show', 'Location', 'none','Box','off','FontSize',10);
lgd.Units = 'normalized';
lgd.Position = [0.13 0.35 0.12 0.3];  % [x y width height]

% Parameters
threshold = 50;
ncolors = 256;                    % size of colormap
colorBelow = [1 1 1];       % RGB for <50 (dark green)
colorAboveMap = parula(ncolors);  % gradient for >=50 (choose any colormap)

% compute how many entries correspond to < threshold
fracBelow = max(0, min(1, (threshold - cb_lim(1)) / (cb_lim(2) - cb_lim(1))));
nbelow = max(1, round(ncolors * fracBelow));
nabove = ncolors - nbelow;

% build colormap: first nbelow rows = colorBelow, rest from colorAboveMap
cmap = [repmat(colorBelow, nbelow, 1); colorAboveMap(end-nabove+1:end,:)];

% Apply to the correlation axes (assume ax(4:6) are your second-row axes)
for k = 4:6
    axes(ax(k));                   % make axis current
    title(sprintf('Correlation Beam %d',k-3))
    axes(ax(k));                   % make axis current
    clim(cb_lim);                   % ensure CLim is consistent
    colormap(ax(k), cmap);         % set custom colormap for this axes
    cb = colorbar(ax(k));          % add colorbar
    cb.Ticks = [cb_lim(1), threshold, cb_lim(2)];        % tick at threshold
    cb.TickLabels = {num2str(cb_lim(1)), num2str(threshold), num2str(cb_lim(2))};
end

% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',fullfile(outdir,[filename,'_f2_beam_amplitude_correlation_QC.png']));
clear ax

%% suface bin detection
figure(3),clf
ax(1) = subplot(1,3,1);
plot(Amp1_pro,y),hold on,plot(Amp2_pro,y),plot(Amp3_pro,y)
yline(y(SB),'--')
yline(y([min(idx_valid(QC_1D==0)),max(idx_valid(QC_1D==0))]),'r--')
xlabel('Temporal mean amplitude [dB]')
ylabel('Nominal cell depth')
axis ij
title([num2str(SB),' valid bins before surface'])

ax(2) = subplot(1,3,2);
plot(Cor1_pro,y),hold on,plot(Cor2_pro,y),plot(Cor3_pro,y)
yline(y(CB),'--')
yline(y([min(idx_valid(QC_1D==0)),max(idx_valid(QC_1D==0))]),'r--')
xline(50,'--')
xlabel('Temporal mean correlation [%]')
title([num2str(CB),' valid bins before surface'])
lgd_text = {'Amp','Cor'};

subplot(1,3,3)
plot(nanmean(U),y),hold on,plot(nanmean(V),y),plot(nanmean(W),y),
ylabel('Nominal cell depth')
xlabel('Velocity [m/s]')
axis ij
grid on
xlim([-0.05 0.2])
yline(y(CB),'--')
yline(y([min(idx_valid(QC_1D==0)),max(idx_valid(QC_1D==0))]),'r--')
xline(0,'--')
legend('U','V','W','surface bins','Sidelobe range','Location','southeast')

for k=1:numel(ax)
    axes(ax(k));    
    ylabel('Nominal cell depth')
    labels = {'Beam1','Beam2','Beam3',[lgd_text{k},' bin range'],...
        'Sidelobe min/max range'};
    legend(labels,'Location','southwest')
    axis ij
    grid on
end
% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',fullfile(outdir,[filename,'_f3_beam_amplitude_correlation_QC.png']));
clear ax

srf_bins = min([SB,CB]);
bins_to_process=input(['\nAutodetected ', num2str(srf_bins),...
    ' valid bins out of ',num2str(max(nCells)),' from the sensor head.',...
    '\nWill flag bins >', num2str(srf_bins),' as bad (',num2str(QC_BAD),')',...
    ' \nDo you want to adjust valid bin number?',...
    ' \nEnter 0 for no (default) or new adjusted number: ']);
if (isempty(bins_to_process) || bins_to_process==0)
    bins_to_process=srf_bins;
end

QC_vel(:,bins_to_process+1:end)=QC_BAD;

fprintf(fidlog,sprintf('Flagged bin >%d as QC_BAD (%d).\n\n',...
    bins_to_process,QC_BAD));

%% sidelobe contamination

prompt = sprintf([ ...
  'Sidelobe contamination for each time step varies between %d and %d ' ...
  'out of %d bins.\nSee red line in figure 2.\n Would you like to flag ' ...
  ' sidelobe contaminated bins for each time step as QC_BAD (%d). Y/N [Y]: '], ...
  min(idx_valid(QC_1D==0))+1, max(idx_valid(QC_1D==0))+1, max(nCells), QC_BAD);

reply = input(prompt, 's');
if isempty(reply)
    reply = 'Y';
end

while true
    if strcmpi(reply, 'Y') || strcmpi(reply, 'YES')
        mask = nCells > idx_valid(:);
        QC_vel(mask) = QC_BAD;

        prombt = ['Sidelobe contamination for each time step varies between ' ...
            ' %d and %d out of %d bins.\nBins contaminated by' ...
            ' sidelobe interference flagged as QC_BAD (%d).\n'];
        fprintf(fidlog, prombt, ...
            min(idx_valid(QC_1D==0))+1, max(idx_valid(QC_1D==0))+1, max(nCells), QC_BAD);
        break

    elseif strcmpi(reply, 'N') || strcmpi(reply, 'NO')
        prombt = ['Sidelobe contamination for each time step varies between ' ...
            'bin %d and %d out of %d bins.\nOperator chose NOT to flag ' ...
            'bins contaminated by sidelobe interference.\n'];
        fprintf(1, prombt, ...
            min(idx_valid(QC_1D==0))+1, max(idx_valid(QC_1D==0))+1, max(nCells)); % to screen
        fprintf(fidlog, prombt, ...
            min(idx_valid(QC_1D==0))+1, max(idx_valid(QC_1D==0))+1, max(nCells));
        break

    else
        fprintf('Invalid entry. Enter Y or N (press Return for default).\n');
        reply = input('Y/N [Y]: ', 's');   % re-prompt and read again
        if isempty(reply)
            reply = 'Y';
        end
    end
end

%% Apply correlation threshold


%% spikes
% find spikes - e.g. fish schools (short lived)
SDc =3;
fprintf('Amplitude 1\n')
mask_spikes_Amp1 = detect_spikes_amp(Amp1, SDc);
fprintf('Amplitude 2\n')
mask_spikes_Amp2 = detect_spikes_amp(Amp2, SDc);
fprintf('Amplitude 3\n')
mask_spikes_Amp3 = detect_spikes_amp(Amp3, SDc);


%% Velocity
% raw
U = Data.Average_VelEast;
V = Data.Average_VelNorth;
W = Data.Average_VelUp;
CSPD = sqrt(U.^2+V.^2);

% QC
W_thr = 0.1;
CUR_thr = 1;

mU = mean(U, 1, 'omitnan'); sU = std(U, 0, 1, 'omitnan');
mV = mean(V, 1, 'omitnan'); sV = std(V, 0, 1, 'omitnan');
mS = mean(CSPD, 1, 'omitnan'); sS = std(CSPD, 0, 1, 'omitnan');

% 2. Create Masks: Checks if each pixel is > 3 STD from its specific depth-mean
% (MATLAB automatically broadcasts the 1xDepth vector across the Time rows)
mask_u = abs(U - mU) > 3.*sU;
mask_v = abs(V - mV) > 3.*sV;
mask_s = abs(CSPD - mS) > 3.*sS;
mask_w = (W > W_thr) | (W < -W_thr);

%% spikes in velocity
U_diff = diff(U, 1, 2);
V_diff = diff(V, 1, 2);
[T,D]=size(U);
mask_ud = [false(T, 1), (abs(U_diff)>CUR_thr)]; 
mask_vd = [false(T, 1), (abs(V_diff)>CUR_thr)]; 


mask_vel_comb = (mask_u | mask_v | mask_s | mask_w | mask_ud | mask_vd);
QC_vel(mask_vel_comb)= QC_BAD;

prombt = ['Velocity QC Summary: Flagged cells as QC_BAD (%d) based on:\n' ...
          ' - Horizontal spikes (|dU/dz|, |dV/dz|) > %0.2f m/s\n' ...
          ' - Statistical outliers > 3 standard deviations (per depth bin)\n' ...
          ' - Vertical velocity outliers (|W|) > %0.2f m/s\n'];

% 1. Print to the command window
fprintf(1, prombt, QC_BAD, CUR_thr, W_thr); 

% 2. Print to the log file
fprintf(fidlog, prombt, QC_BAD, CUR_thr, W_thr);


U_QC = U; U_QC(QC_vel==4)=NaN;
V_QC = V; V_QC(QC_vel==4)=NaN;
W_QC = W; W_QC(QC_vel==4)=NaN;
CSPD_QC = CSPD; CSPD_QC(QC_vel==4)=NaN;

f4 = figure(4); clf;

W_thresh = 0.1; 
cb_titles = {'Zonal Velocity (m/s)', 'Zonal Velocity (m/s)', ...
          'Meridional Velocity (m/s)', 'Meridional Velocity (m/s)',...
          'Vertical Velocity (m/s)', 'Vertical Velocity (m/s)', ...
          'Current Speed (m/s)', 'Current Speed (m/s)'};
data_list = {U,U_QC,V,V_QC,W,W_QC,CSPD,CSPD_QC};

% High-Contrast Colormap: [Cyan; Blue-White-Red; Yellow]
b_to_w = [linspace(0,1,11)', linspace(0,1,11)', ones(11,1)];
w_to_r = [ones(10,1), linspace(0.9,0,10)', linspace(0.9,0,10)'];
cmap = [[0 1 1]; b_to_w; w_to_r; [1 1 0]];

for k = 1:8
    ax = subplot(4,2,k);
    im = imagesc(ax, x, y, data_list{k}');
    hold(ax, 'on');
    
    % Visual Formatting
    set(ax, 'YDir', 'reverse', 'Color', [0.8 0.8 0.8]); % Gray for NaNs
    set(im, 'AlphaData', ~isnan(data_list{k}'));       % Transparency for NaNs
    datetick(ax, 'x', 'mmm-yyyy', 'keepticks', 'keeplimits');
    ylim(ax, [0, 1080]);
    xlim(ax, [min(x), max(x)]);

    colormap(ax, cmap);
    
    if k==1
        title('Raw velocities')
    elseif k==2
        title_str = {
    sprintf('QC velocities: |U|,|V|,|CSPD| < mean+3*STD; |dU/dz|,|dV/dz|<1 m/s'), ...
    sprintf('|W|<0.1 m/s, and previous QCs')
};
        title(title_str)
    end

    if k == 5 || k== 6% Vertical Velocity logic
        limit_W = W_thresh * 1.1;
        clim(ax, [-limit_W, limit_W]);
        cb = colorbar(ax);
        cb.Ticks = [-limit_W, -W_thresh, 0, W_thresh, limit_W];
        cb.TickLabels = {sprintf('<%.2f',-W_thresh), '-0.1', '0', '0.1', sprintf('>%.2f',W_thresh)};
    else
        clim(ax, [-1.1, 1.1]);
        cb = colorbar(ax);
        cb.Ticks = [-1.1, -1, 0, 1, 1.1];
        cb.TickLabels = {'<-1', '-1', '0', '1', '>1'};
    end
    ylabel(cb, cb_titles{k});
end

% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',fullfile(outdir,[filename,'_f4_velocity_and_speed_QC.png']));
clear ax
%%
end
fprintf(fidlog, '\n==== END ENTRY  =====\n');
fclose(fidlog);
end

%% nested fuctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% looks for deployment and recovery period in pressure
function S = suggestTrimForShallowEdges_simple(y, ymin,x)
% Suggest trim indices for leading/trailing consecutive y < ymin.
% S.startIdx = first index to keep
% S.endIdx   = last  index to keep
% S.keepMask

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

startIdx = firstKeep;
endIdx   = lastKeep;

% Build masks
maskStart = false(n,1); if nTrimStart>0, maskStart(1:firstKeep-1) = true; maskStart(~isShallow) = false; end
maskEnd   = false(n,1); if nTrimEnd>0, maskEnd(lastKeep+1:end) = true; maskEnd(~isShallow) = false; end
keepMask = false(n,1);
if firstKeep <= lastKeep
    keepMask(firstKeep:lastKeep) = true;
end

S = struct('startIdx', startIdx, 'endIdx', endIdx, ...
           'nTrimStart', nTrimStart, 'nTrimEnd', nTrimEnd, ...
           'keepMask', keepMask);
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

function [mu_time,s_time,B] = calculate_block_mean_std(A)
    % A: time_dim x depth
    block = 10;
    [T, D] = size(A);
    nBlocks = T/block;
    % reshape to (block, nBlocks, depth)
    B = reshape(A, block, nBlocks, D);
    % mean and std over the first dimension (within each 10-s burst)
    mu = median(B, 1);    % 1 x nBlocks x D
    s  = std(B, 0, 1);  % 1 x nBlocks x D  (default normalization N-1)
    
    % replicate to (block, nBlocks, D)
    mu_rep = repmat(mu, block, 1, 1);
    s_rep  = repmat(s,  block, 1, 1);
    
    % reshape back to (time_dim_trimmed, depth)
    mu_time = reshape(mu_rep, nBlocks*block, D);
    s_time  = reshape(s_rep,  nBlocks*block, D);
end

function mask = detect_spikes_amp(A, SDc)
    % A: time_dim x depth
    block = 10;
    [T, D] = size(A);
    nBlocks = T/block;

    % --- Criterion 1: Amplitude spike along Time Dimension ---
    [A_bm, A_bs, B] = calculate_block_mean_std(A);
    mask_amp = (A > (A_bm + SDc*A_bs)) | (A < (A_bm - SDc*A_bs));

    % --- Criterion 2: Sudden Amplitude Increase along depth (QARTOD) ---
    % B is (10 x nBlocks x D)
    B_md = median(B, 1);    % Result: 1 x nBlocks x D
    
    % Diff along Depth (Dim 3 of the B_md matrix)
    Z = diff(B_md, 1, 3);   % Result: 1 x nBlocks x (D-1)

    % Running stats along nBlocks (Time), which is Dim 2
    window_size = 3;
    Z_md_block = movmedian(Z, window_size, 2);
    Z_sd_block = movstd(Z, window_size, 0, 2);

    % Create the mask for Z (1 x nBlocks x D-1)
    mask_z_block = (Z>(Z_md_block+SDc*Z_sd_block)) | (Z<(Z_md_block-SDc*Z_sd_block));

    % Replicate to match the 10 samples per block
    % Result: 10 x nBlocks x D-1
    mask_z_rep = repmat(mask_z_block, block, 1, 1);

    % Reshape back to (T x D-1)
    mask_z_2d = reshape(mask_z_rep, T, D-1);

    % Pad with false at the first depth column to return to (T x D)
    mask_z = [false(T, 1), mask_z_2d]; 

    % --- Final Combined Mask ---
    mask = mask_amp | mask_z;

    % Statistics
    total_points = numel(A);
    num_outliers = sum(mask(:));
    if num_outliers>0
    fprintf('\nTotal Amplitude Spikes Found: %d (%.2f%% of data)\n', num_outliers, (num_outliers/total_points)*100);
    fprintf('Consider smoothing velocity ensembles\n\n')
    end
end

% plot U,W,V and current speed
function plot_UVW(parent_container, U, V, W, CSPD, x, y)
    % Create layout inside the Panel
    t = tiledlayout(parent_container, 2, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
    
    W_thresh = 0.1; 
    titles = {'Zonal Velocity (m/s)', 'Meridional Velocity (m/s)', ...
              'Vertical Velocity (m/s)', 'Current Speed (m/s)'};
    data_list = {U, V, W, CSPD};
    
    % High-Contrast Colormap: [Cyan; Blue-White-Red; Yellow]
    b_to_w = [linspace(0,1,11)', linspace(0,1,11)', ones(11,1)];
    w_to_r = [ones(10,1), linspace(0.9,0,10)', linspace(0.9,0,10)'];
    cmap = [[0 1 1]; b_to_w; w_to_r; [1 1 0]];

    for k = 1:4
        ax = nexttile(t);
        im = imagesc(ax, x, y, data_list{k}');
        hold(ax, 'on');
        
        % Visual Formatting
        set(ax, 'YDir', 'reverse', 'Color', [0.8 0.8 0.8]); % Gray for NaNs
        set(im, 'AlphaData', ~isnan(data_list{k}'));       % Transparency for NaNs
        datetick(ax, 'x', 'mmm-yyyy', 'keepticks', 'keeplimits');
        ylim(ax, [0, 1080]);
        xlim(ax, [min(x), max(x)]);

        colormap(ax, cmap);
        
        if k == 3 % Vertical Velocity logic
            limit_W = W_thresh * 1.1;
            clim(ax, [-limit_W, limit_W]);
            cb = colorbar(ax);
            cb.Ticks = [-limit_W, -W_thresh, 0, W_thresh, limit_W];
            cb.TickLabels = {sprintf('<%.2f',-W_thresh), '-0.1', '0', '0.1', sprintf('>%.2f',W_thresh)};
        else
            clim(ax, [-1.1, 1.1]);
            cb = colorbar(ax);
            cb.Ticks = [-1.1, -1, 0, 1, 1.1];
            cb.TickLabels = {'<-1', '-1', '0', '1', '>1'};
        end
        ylabel(cb, titles{k});
    end
end