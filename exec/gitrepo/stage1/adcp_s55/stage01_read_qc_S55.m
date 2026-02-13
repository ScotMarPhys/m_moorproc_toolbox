% STAGE01_READ_QC_S55 Read and perform initial QC on Signature 55 ADCP data.
%   This function processes Nortek Signature 55 .mat files, performing 
%   automated flagging for pressure, tilt, signal strength, sidelobe interference,
%   and plots magnetometer hard-iron distortions and generell sensor
%   diagnostics
%
%   USAGE:
%      stage01_read_qc_S55('rhadcp_01_2020')
%      stage01_read_qc_S55(moor, dataindir, filename, infofile, logfile, outdir)
%
%   INPUTS:
%      moor      - Mooring name (string), e.g., 'rhadcp_01_2020'
%      dataindir - (Optional) Path to raw .mat files
%      filename  - (Optional) Specific input filename
%      infofile  - (Optional) Path to mooring info file
%      logfile   - (Optional) Path for output log
%      outdir    - (Optional) Path for output figures, data, and logs
%
%   QC CHECKS PERFORMED:
%      1. Pressure: Validity check for deployment/retrieval periods.
%         Orientation: Tilt thresholds (<10 deg good, >30 deg fail).
%         Heading
%      2. Signal/Corr: Surface bin dectection and noise floor (Corr < 50 flagged).
%      3. Sidelobe: Range masking based on depth (H) and max tilt (beta).
%      4. Velocity: Horizontal Spikes, statistical outliers (U,V) and unrealistic
%                   velocities (U,V,W).
%      5. Compass: Magnetometer circle-fit diagnostics (Hard-Iron check).
%      5. Other sensor diagnostics.
%
%   OUTPUTS:
%      - Stage 1 NetCDF file ([moor_id]_[sn]_stage1.nc)
%      - 6 Diagnostic QC Figures
%
%   DEPENDENCIES:
%      - [m_moorproc_toolbox](URL_TO_YOUR_REPO) (for rodbload.m)
%      - [GSW Oceanographic Toolbox](https://www.teos-10.org) (v3.06+)
%
%   NOTE: This is a Stage 1 function. It FLAGS data only; physical 
%   rotations and corrections (e.g., Hard-Iron U/V rotation) occur in Stage 2. 



function stage01_read_qc_S55(moor, varargin)

if nargin == 0
    help(mfilename);
    return;
end

global MOORPROC_G

if nargin == 1 && ~isempty(MOORPROC_G)
    operator    = MOORPROC_G.operator;
    pd          = moor_inoutpaths('adcp_S55', moor);
    dataindir   = pd.rawpath;
    infofile    = pd.infofile;
    logfile     = pd.stage1log;
    outdir      = pd.stage1path;
    output_form = pd.stage1form;
else
    try
        operator = getenv('COMPUTERNAME');    
        dataindir = varargin{1};
        filename = varargin{2};
        infofile = varargin{3};
        logfile = varargin{4};
        outdir = varargin{5};
        ouput_form = [moor '_%d_stage1.nc'];
    catch
        error('Not enough manual arguments provided. Expected 5 additional arguments after "moor".');
    end
end

% Directory and Logfile Validation
if ~exist(outdir, 'dir'), mkdir(outdir); end

[~, log_name, log_ext] = fileparts(logfile);
expected_log = [moor '_ADCP_stage1'];

if ~strcmp([log_name log_ext], [expected_log '.log'])
    error('LOGFILE MISMATCH: Expected %s.log but got %s%s.', expected_log, log_name, log_ext);
end

% Metadata and Log Initialization
info_adcp = read_adcp_infofile(infofile);

% Timestamp conversion using simpler datetime logic
info_adcp.start_timestamp = datenum(datetime([info_adcp.s_d.' , info_adcp.s_t.' , 0])); %start date and time bg
info_adcp.end_timestamp   = datenum(datetime([info_adcp.e_d.' , info_adcp.e_t.' , 0])); %end date and time ed

% Open Logfile
fidlog = fopen(logfile, 'w');
if fidlog == -1, error('Permission denied: Could not open %s', logfile); end

fprintf(fidlog, '\n==== START ENTRY  =====\n');
fprintf(fidlog, 'Read and quality control Signature 55 ADCP data. \n');
fprintf(fidlog, 'Processing carried out by %s at %s\n\n\n', operator, datestr(clock));
fprintf(fidlog, 'Mooring   %s \n', moor);
fprintf(fidlog, 'Latitude  %6.3f \n', info_adcp.lat);
fprintf(fidlog, 'Longitude %6.3f \n\n\n', info_adcp.lon);

% Serial Number Logic
if isnan(info_adcp.id)
    fprintf('No serial number given, using default 200044\n');
    serial_nums = 200044;
else
    % Filters IDs between 319 and 328
    mask = (info_adcp.id >= 319) & (info_adcp.id <= 328);
    serial_nums = info_adcp.sn(mask);
end

%% Load data: Config, Data, Description
for i = 1:length(serial_nums)
fprintf('Processing sn %d\n', serial_nums(i));

% Handle filename if not pre-defined
if ~exist('filename', 'var') || isempty(filename)
    current_filename = sprintf('%d_data', serial_nums(i));
else
    current_filename = filename;
end

outfile = fullfile(outdir, sprintf(ouput_form, serial_nums(i)));
infile  = fullfile(dataindir, [current_filename, '.mat']);

% Load and clean structure
load(infile); 
clear Descriptions % Not helpful for post-processing

% Efficiently convert all struct fields to double (Improved syntax)
Data = structfun(@(x) double(x), Data, 'UniformOutput', false);

fprintf(fidlog, 'infile : %s\n', infile);
fprintf(fidlog, 'ADCP serial number  : %d\n', serial_nums(i));

% Build post-processing flag array
% Reference values (Consistent with original)
QC_NOT_EVALUATED = 0; 
QC_GOOD = 1; 
QC_UNKNOWN = 2; 
QC_PROBABLY_BAD = 3;
QC_BAD = 4; 
QC_CHANGED = 5; 
QC_UNSAMPLED = 6; 
QC_INTERPOLATED = 7;
%QC_COMPASS_BAD = 8; 
QC_MISSING = 9;

flag_vals = [0, 1, 2, 3, 4, 5, 6, 7, 9];
flag_mean = ['1=QC_NOT_EVALUATED 1=QC_GOOD ' ...
    '3=QC_UNKNOWN QC_PROBABLY_BAD 4=QC_BAD QC_CHANGED ' ...
    '5=QC_UNSAMPLED 6=QC_INTERPOLATED ' ...
    '9=QC_MISSING'];

% Initialise QC masks using zeros (More efficient than 0 * double)
Data.mask_QC_3D = zeros(size(Data.Average_VelEast));
Data.mask_QC_1D = zeros(size(Data.Average_Pressure));

% Checking Serial Number Consistency
if serial_nums(i) ~= Config.SerialNo
    error('Entered Serial Number %d does not match Config %d', ...
        serial_nums(i), Config.SerialNo);
end

% Checking Configuration Consistency
check_fields = {'Average_NBeams', 'Average_NCells', 'Average_BeamToChannelMapping', ...
                'Average_AmbiguityVel', 'Average_NominalCorrelation', 'Average_Soundspeed'};
check_labels = {'Number of beams', 'Number of cells', 'Beam to channel mapping', ...
                'Ambiguity velocity', 'Nominal Correlation', 'Sound speed'};

msg_suffix = ' is not equal for all time steps. Please check ';

for k = 1:length(check_fields)
    if size(unique(Data.(check_fields{k}), 'rows'), 1) ~= 1
        prombt_2 = [check_labels{k}, msg_suffix, 'Data.', check_fields{k}];
        disp(prombt_2); fprintf(fidlog, '%s\n', prombt_2);
    end
end

% Specific Check for Average Error
if size(unique(Data.Average_Error, 'rows'), 1) ~= 1 || Data.Average_Error(1) ~= 0
    prombt_2 = 'Errors occurred during measuring. Please check Data.Average_Error';
    disp(prombt_2); fprintf(fidlog, '%s\n', prombt_2);
end

%% STAGE 1.  Nortek suggested quality control steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Load and Rename Variables ---
timemin = Data.Average_Time(1);
timemax = Data.Average_Time(end);

T     = Data.Average_Time;
P     = Data.Average_Pressure; 
B     = Data.Average_Battery;
Pitch = Data.Average_Pitch;
Roll  = Data.Average_Roll;
H     = Data.Average_Heading;

% pressure stats and ylims
P_median = median(P);
P_std    = std(P);
n_std    = 1;

%% figure settings
fs = 14;
set(findall(gcf, '-property', 'FontSize'), 'FontUnits', 'points', 'FontSize', fs);

% ylims
yl_P     = (P_median + n_std * [-P_std P_std]);
yl_Pitch =  [-40 40];
yl_H = [-0 360];

%% 1. Instrument orientation and pressure %%%%%%%%%%%%%%%%%%%%%%%%
% --- 1. Pressure Check (Sinking/Rising) ---
S = suggestTrimForShallowEdges_simple(P, min(yl_P));

% Update Masks for Pressure
Data.mask_QC_1D(S.keepMask == 0) = QC_BAD; 
Data.mask_QC_3D(S.keepMask == 0, :) = QC_BAD;

z_mean = gsw_z_from_p(mean(P(Data.mask_QC_1D==0)),info_adcp.lat);

% Pitch and Roll
Pitch_Pgood = Pitch; Pitch_Pbad = Pitch;
Pitch_Pgood(S.keepMask == 0) = NaN;Pitch_Pbad(S.keepMask == 1) = NaN;

Roll_Pgood = Roll; %Roll_Pbad = Roll;
Roll_Pgood(S.keepMask == 0) = NaN; %Roll_Pbad(S.keepMask == 0) = NaN;

% Output logging for Pressure
prompt = sprintf(['\n***Pressure check***\n', ...
    'Found %d shallow values at start and %d at end flagged bad (%d).\n\n'], ...
    S.nTrimStart, S.nTrimEnd, QC_BAD);
fprintf(1, prompt); % Print to Command Window
fprintf(fidlog, prompt); % Print to Log File
clear prombt

% --- 2. Information File Check ---
ds_vec = datevec(T([S.startIdx, S.endIdx])); 
bg_suggest = datenum([ds_vec(1, 1:end-1), 0]);
ed_suggest = datenum([ds_vec(2, 1:end-1), 0]);

if bg_suggest~=info_adcp.start_timestamp
    prombt = ['***',newline,'Suggest to set mooring start date and time in ',... 
            moor,'info.dat to ',datestr(bg_suggest),newline,'***',newline,newline];
    disp(prombt); fprintf(fidlog,prombt); clear prombt
elseif ed_suggest~=info_adcp.end_timestamp
        prombt = ['***',newline,'Suggest to set mooring end date and time in ',... 
            moor,'info.dat to ',datestr(ed_suggest),newline,'***',newline,newline];
    disp(prombt); fprintf(fidlog,prombt);clear prombt
end

% --- 3. Pitch and Roll Check ---
% Identify Excessive Pitch/Roll (> 30)
badind_tilt_severe = find((Pitch>30 | Pitch<-30 | Roll>30 | Roll<-30) ...
    & Data.mask_QC_1D == 0);
Data.mask_QC_1D(badind_tilt_severe) = QC_BAD;
Data.mask_QC_3D(badind_tilt_severe, :) = QC_BAD;

% Identify Moderate Pitch/Roll (10 to 30)
badind_tilt_mod = find(((abs(Pitch) > 10 & abs(Pitch) <= 30) | ...
               (abs(Roll) > 10 & abs(Roll) <= 30)) & Data.mask_QC_1D == 0);
Data.mask_QC_1D(badind_tilt_mod) = QC_PROBABLY_BAD;
Data.mask_QC_3D(badind_tilt_mod, :) = QC_PROBABLY_BAD;

% Logging for Tilt
prompt = sprintf(['***Pitch check***\n', ...
    '%d timesteps flagged bad (%d) due to excessive pitch (>abs(30))\n', ...
    '%d timesteps flagged probably bad (%d) due to pitch between 10 and 30\n', ...
    'Post processing possible.\n\n'], ...
    length(badind_tilt_severe), QC_BAD, ...
    length(badind_tilt_mod), QC_PROBABLY_BAD);
fprintf(1, prompt); fprintf(fidlog, prompt); clear prombt

% Heading mark values with bad pitch/roll/pressure
H_trim = H;
H_trim(Data.mask_QC_1D == 0) = NaN;

% --- Visualization ---
%% 1. Instrument orientation and pressure %%%%%%%%%%%%%%%%%%%%%%%%
figure(1);clf

% Pressure %%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,1)
title([filename ' pressure'],'Interpreter','none');
hold on; grid on;

h.data = plot(T,P,'k.');
yline(P_median,'r')
[h, nAbove, nBelow] = markClippedPoints_pres(h, P, yl_P, T);

datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');
xlim([timemin timemax]);
ylim(yl_P);
ylabel('Pressure (db)');
set(gca,'YDir','reverse');

% --- RIGHT AXIS: Battery ---
yyaxis right
hold on
h_bat = plot(T, B, 'b-', 'LineWidth', 1); 
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

bx = axPos(1) + 0.05*axPos(3);
by = axPos(2) + 0.05*axPos(4);
bw = 0.32*axPos(3);
bh = 0.28*axPos(4);

hBox = annotation('textbox', [bx by bw bh], 'String', sprintf('Z_{mean} = %4.0f m', z_mean), ...
    'FitBoxToText', 'on', 'EdgeColor', 'none', ...
    'Interpreter', 'tex', 'HorizontalAlignment', 'left');


%%% Roll and pitch %%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2) 
title([filename ' pitch and roll'],'Interpreter','none');
hold on; grid on;

% plot data
h.data = plot(T,Pitch_Pgood,'.','DisplayName', 'pitch raw bottom');
plot(T,Pitch_Pbad,'.c','DisplayName', 'pitch QC_BAD sinking/rising');
plot(T,Roll_Pgood,'k','DisplayName', 'roll raw bottom','LineWidth',0.5);

datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');

% Suggested quality thresholds
hl1 = yline(0,'color','k'); hl1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hl1 = yline([-10 10],'color','g');  arrayfun(@(h) set(h.Annotation.LegendInformation, ...
    'IconDisplayStyle', 'off'), hl1);
hl1 = yline([-30 30],'color','r');  arrayfun(@(h) set(h.Annotation.LegendInformation, ...
    'IconDisplayStyle', 'off'), hl1);

xlim([timemin timemax]);
ylim(yl_Pitch);

ys_data = [ 10, -10, 30, -30 ];
ys_labels = { '10 < 30: Post processing possible', ...
           '-10 > -30: Post processing possible', ...
           '>30: Profiles bad', ...
           '<-30: Profiles bad' };
ys_colors =  {'g','g','r','r'};

annotateYLabels(ys_data, ys_labels, ys_colors);

ylabel('Pitch (^o)');
legend('Interpreter','none','Location','southeast')


% heading %%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,3)
title([filename ' heading'],'Interpreter','none');
hold on; grid on;

% Plot data
plot(T,H,'k.','DisplayName', 'heading raw bottom');
plot(T,H_trim,'.c','DisplayName', 'QC_BAD sinking/rising/pitch');

datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');
xlim([timemin timemax]);

ylim(yl_H);
ylabel('Heading (^o)');

legend('Interpreter','none','Location','southeast')

% add remarks
% place at top-right of figure (normalized units)
str = 'Pitch and heading should remain fairly constant during deployment';
annotation('textbox', [0 0.02 1 0.06], ...   % [x y w h] in normalized figure units
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

if ~isnan(info_adcp.wd)
    Nominal_CellDepth = info_adcp.wd-Dist2Instr_CellMidpoint;
else
    info_adcp.wd = 1080;
    while true
        prompt = sprintf('Please enter nominal depth of instrument in m (default %dm): ', info_adcp.wd);
        val = input(prompt);
        if isempty(val)                     % user accepted default
            info_adcp.wd = info_adcp.wd;
            break
        end
        if isnumeric(val) && isscalar(val) && isfinite(val) && val > 0
            info_adcp.wd = val;
            break
        end
        fprintf('Invalid entry. Enter a positive number or press Return for default.\n');
    end
    Nominal_CellDepth = info_adcp.wd-Dist2Instr_CellMidpoint;
end

fprintf(fidlog,'***Beams, surface bins, sidelobe check***\n');
fprintf(fidlog,sprintf('Nominal depth of instrument set as %d dbar.\n',info_adcp.wd));

y=Nominal_CellDepth;
cb_lim =[0,100];

% define surface bins
Amp1=Data.Average_AmpBeam1; Amp1(Data.mask_QC_1D~=0,:)=NaN;
Amp2=Data.Average_AmpBeam2; Amp2(Data.mask_QC_1D~=0,:)=NaN;
Amp3=Data.Average_AmpBeam3; Amp3(Data.mask_QC_1D~=0,:)=NaN;
Cor1=Data.Average_CorBeam1; Cor1(Data.mask_QC_1D~=0,:)=NaN;
Cor2=Data.Average_CorBeam2; Cor2(Data.mask_QC_1D~=0,:)=NaN;
Cor3=Data.Average_CorBeam3; Cor3(Data.mask_QC_1D~=0,:)=NaN;

% mask any Cor<50 as bad
mask_corr = (Cor1 < 50) | (Cor2 < 50) | (Cor3 < 50);
Data.mask_QC_3D(mask_corr)=QC_BAD;

U = Data.Average_VelEast;U(Data.mask_QC_1D~=0,:)=NaN;
V = Data.Average_VelNorth;V(Data.mask_QC_1D~=0,:)=NaN;
W = Data.Average_VelUp;W(Data.mask_QC_1D~=0,:)=NaN;

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
A(1) = gsw_z_from_p(nanmean(Data.Average_Pressure(Data.mask_QC_1D==0)),info_adcp.lat); % range based on pressure sensor
A(2) = Dist2Instr_CellMidpoint(SM); % range based on nominal cell distance

%all angle in degrees
gamma = 20; %slant angle
alpha = [0, 120,240]; % beam 1 - aligned with x-axis, beam 2, beam 3
pitch = Data.Average_Pitch; pitch(Data.mask_QC_1D~=0)=NaN; % y-axis rotation, in deg
roll = Data.Average_Roll; roll(Data.mask_QC_1D~=0)=NaN; % x-axis rotation, in deg

% Pre-calculate trig terms (using 'd' versions for degrees)
cp = cosd(pitch); sp = sind(pitch);
cr = cosd(roll);  sr = sind(roll);
cg = cosd(gamma); sg = sind(gamma);
ca = cosd(alpha); sa = sind(alpha); % alpha is 1x3

% Angle to vertical for each beam (in degrees)
beta = acosd(-sp .* sg .* ca + cp .* sr .* sg .* sa + cp .* cr .* cg);

% Find the minimum angle to vertical across all three beams for each timestamp
max_beta = max(beta, [], 2);

% R contains buffer for upper range of cell
R =max(A)*cosd(max_beta)-Cell_Size; 
% find(Dist2Instr_CellMidpoint <= R(i), 1, 'last') for each time step i
R(Data.mask_QC_1D~=0)=info_adcp.wd;
idx_valid = arrayfun(@(r) find(Dist2Instr_CellMidpoint <= r, 1, 'last'),...
    R);
idx_valid(Data.mask_QC_1D~=0)=0;
R(Data.mask_QC_1D~=0)=NaN;


% Visualisation Amplides and Correlation
f2 = figure(2);clf
ax(1) = subplot(2,3,1);hold on; imagesc(T, y, Amp1');
ax(2) = subplot(2,3,2);hold on; imagesc(T, y, Amp2');
ax(3) = subplot(2,3,3);hold on; imagesc(T, y, Amp3');

ax(4) = subplot(2,3,4);hold on; imagesc(T, y, Cor1');
ax(5) = subplot(2,3,5);hold on; imagesc(T, y, Cor2');
ax(6) = subplot(2,3,6);hold on; imagesc(T, y, Cor3');

for k = 1:numel(ax)
    axis(ax(k),'xy','ij');
    ylim(ax(k),[0,1080])
    ylabel(ax(k),'Nominal cell depth (m)')
    datetick(ax(k),'x','mmm-yyyy','keepticks','keeplimits');
    hLine = plot(ax(k),T(Data.mask_QC_1D==0), y(idx_valid(Data.mask_QC_1D==0)), 'r', 'LineWidth', 0.1);
    hLine.DisplayName = 'Sidelobe interference';
    xlim(ax(k),[min(T),max(T)])
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

%% 3. suface bin detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3),clf
ax(1) = subplot(1,3,1);
plot(Amp1_pro,y),hold on,plot(Amp2_pro,y),plot(Amp3_pro,y)
yline(y(SB),'--','LineWidth',1.5)
yline(y([min(idx_valid(Data.mask_QC_1D==0)),max(idx_valid(Data.mask_QC_1D==0))]),'r--')
xlabel('Temporal mean amplitude [dB]')
ylabel('Nominal cell depth')
axis ij
title([num2str(SB),' valid bins before surface'])

ax(2) = subplot(1,3,2);
plot(Cor1_pro,y),hold on,plot(Cor2_pro,y),plot(Cor3_pro,y)
yline(y(CB),'--','LineWidth',1.5)
yline(y([min(idx_valid(Data.mask_QC_1D==0)),max(idx_valid(Data.mask_QC_1D==0))]),'r--')
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
yline(y([min(idx_valid(Data.mask_QC_1D==0)),max(idx_valid(Data.mask_QC_1D==0))]),'r--')
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

%% prombts after plotting
% surface bin
srf_bins = min([SB,CB]);
bins_to_process=input(['\n\n***Surface bins and sidelobe interference***', ...
    '\nAutodetected ', num2str(srf_bins),...
    ' valid bins out of ',num2str(max(nCells)),' from the sensor head.',...
    '\nWill flag bins >', num2str(srf_bins),' as bad (',num2str(QC_BAD),')',...
    ' \nDo you want to adjust valid bin number?',...
    ' \nEnter 0 for no (default) or new adjusted number: ']);
if (isempty(bins_to_process) || bins_to_process==0)
    bins_to_process=srf_bins;
end

Data.mask_QC_3D(:,bins_to_process+1:end)=QC_BAD;

fprintf(fidlog,sprintf('Flagged bin >%d as QC_BAD (%d).\n\n',...
    bins_to_process,QC_BAD));

% sidelobe contamination

prompt = sprintf([ ...
  '\n\nSidelobe contamination for each time step varies between %d and %d ' ...
  'out of %d bins.\nSee red line in figure 2.\n Would you like to flag ' ...
  ' sidelobe contaminated bins for each time step as QC_BAD (%d). Y/N [Y]: '], ...
  min(idx_valid(Data.mask_QC_1D==0))+1, max(idx_valid(Data.mask_QC_1D==0))+1, max(nCells), QC_BAD);

reply = input(prompt, 's');
if isempty(reply)
    reply = 'Y';
end

while true
    if strcmpi(reply, 'Y') || strcmpi(reply, 'YES')
        mask = nCells > idx_valid(:);
        Data.mask_QC_3D(mask) = QC_BAD;

        prombt = ['Sidelobe contamination for each time step varies between ' ...
            ' %d and %d out of %d bins.\nBins contaminated by' ...
            ' sidelobe interference flagged as QC_BAD (%d).\n'];
        fprintf(fidlog, prombt, ...
            min(idx_valid(Data.mask_QC_1D==0))+1, max(idx_valid(Data.mask_QC_1D==0))+1, max(nCells), QC_BAD);
        break

    elseif strcmpi(reply, 'N') || strcmpi(reply, 'NO')
        prombt = ['Sidelobe contamination for each time step varies between ' ...
            'bin %d and %d out of %d bins.\nOperator chose NOT to flag ' ...
            'bins contaminated by sidelobe interference.\n'];
        fprintf(1, prombt, ...
            min(idx_valid(Data.mask_QC_1D==0))+1, max(idx_valid(Data.mask_QC_1D==0))+1, max(nCells)); % to screen
        fprintf(fidlog, prombt, ...
            min(idx_valid(Data.mask_QC_1D==0))+1, max(idx_valid(Data.mask_QC_1D==0))+1, max(nCells));
        break

    else
        fprintf('Invalid entry. Enter Y or N (press Return for default).\n');
        reply = input('Y/N [Y]: ', 's');   % re-prompt and read again
        if isempty(reply)
            reply = 'Y';
        end
    end
end


% spikes
% find spikes - e.g. fish schools (short lived)
% SDc =3;
% fprintf('Amplitude 1\n')
% mask_spikes_Amp1 = detect_spikes_amp(Amp1, SDc);
% fprintf('Amplitude 2\n')
% mask_spikes_Amp2 = detect_spikes_amp(Amp2, SDc);
% fprintf('Amplitude 3\n')
% mask_spikes_Amp3 = detect_spikes_amp(Amp3, SDc);


%% 4. Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% spikes in velocity
U_diff = diff(U, 1, 2);
V_diff = diff(V, 1, 2);
[Tidx,Didx]=size(U);
mask_ud = [false(Tidx, 1), (abs(U_diff)>CUR_thr)]; 
mask_vd = [false(Tidx, 1), (abs(V_diff)>CUR_thr)]; 


mask_vel_comb = (mask_u | mask_v | mask_s | mask_w | mask_ud | mask_vd);
Data.mask_QC_3D(mask_vel_comb)= QC_BAD;

prombt = ['\n\nVelocity QC Summary: Flagged cells as QC_BAD (%d) based on:\n' ...
          ' - Horizontal spikes (|dU/dz|, |dV/dz|) > %0.2f m/s\n' ...
          ' - Statistical outliers > 3 standard deviations (per depth bin)\n' ...
          ' - Vertical velocity outliers (|W|) > %0.2f m/s\n'];
fprintf(1, prombt, QC_BAD, CUR_thr, W_thr); 
fprintf(fidlog, prombt, QC_BAD, CUR_thr, W_thr);


U_QC = U; U_QC(Data.mask_QC_3D~=0)=NaN;
V_QC = V; V_QC(Data.mask_QC_3D~=0)=NaN;
W_QC = W; W_QC(Data.mask_QC_3D~=0)=NaN;
CSPD_QC = CSPD; CSPD_QC(Data.mask_QC_3D~=0)=NaN;

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
    im = imagesc(ax, T, y, data_list{k}');
    hold(ax, 'on');
    
    % Visual Formatting
    set(ax, 'YDir', 'reverse', 'Color', [0.8 0.8 0.8]); % Gray for NaNs
    set(im, 'AlphaData', ~isnan(data_list{k}'));       % Transparency for NaNs
    datetick(ax, 'x', 'mmm-yyyy', 'keepticks', 'keeplimits');
    ylim(ax, [0, 1080]);
    xlim(ax, [min(T), max(T)]);

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
        cb.TickLabels = {sprintf('<%.2f',-W_thresh), '-0.1', '0', '0.1', ...
sprintf('>%.2f',W_thresh)};
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

%% Magnetometer Horizontal Intensity and Circle Fitting %%%%%%%%%%%%%%%%%%
Data_corr = s55_hard_iron_compass_correction2(Data,0,false);
err_min = min(Data_corr.hard_iron_CCW_angle,[],'omitnan');
err_max = max(Data_corr.hard_iron_CCW_angle,[],'omitnan');

% Logfile output
deg = char(176); 
prombt = ['\n\n***Hard-Iron compass correction***\n' ...
    'Simple circle fit shows error angle varies between %3.0f' deg ... 
    'and %3.0f' deg '\n'];
fprintf(1, prombt, err_min, err_max); 
fprintf(fidlog, prombt, err_min, err_max); 

% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',fullfile(outdir,[filename,'_f5_horizontal_magnetometer_QC.png']));
%% 6. Sensor diagnostics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f6 = figure(6); clf;

% --- TOP ROW: ALL TEMPERATURES GROUPED ---
ax(1) = subplot(4,1,1);
hold on;
plot(T,Data.Average_Temperature, 'LineWidth', 1.5, 'DisplayName', 'Water Temp');
plot(T,Data.Average_RTCTemperature, 'DisplayName', 'RTC (Internal)');
plot(T,Data.Average_MagnetometerTemperature, 'DisplayName', 'Mag Temp');
plot(T,Data.Average_PressureSensorTemperature, 'DisplayName', 'Pressure Temp');
hold off;
title('Temperature Sensors');
ylabel('°C');
legend('Location', 'north','NumColumns',4);
grid on;

% --- MIDDLE ROW: MECHANICAL/ENVIRONMENTAL ---
ax(2) = subplot(4,1,2);
plot(T,Data.Average_Magnetometer);
title('Magnetometer'); grid on;
ylabel('mG');

ax(3) = subplot(4,1,3);
plot(T,Data.Average_Accelerometer);
title('Accelerometer'); grid on;
ylabel('g')

% --- BOTTOM ROW: SENSOR STATUS ---
ax(4) = subplot(4,1,4);
plot(T,Data.Average_TransmitEnergy);
title('Transmit Energy'); grid on;
ylabel('Joules');

for k = 1:numel(ax)
    axis(ax(k));
    datetick(ax(k),'x','mmm-yyyy','keepticks','keeplimits');
    xlim(ax(k),[min(T)-1e1,max(T)+1e1])
end

% Improve spacing
sgtitle('ADCP Sensor Diagnostic');
% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',fullfile(outdir,[filename,'_f6_sensor_diagnostics.png']));
clear ax
%%
m_u=median(U_QC,"omitnan" )*1e2;m_v=median(V_QC,"omitnan" )*1e2;
m_w=median(W_QC,"omitnan" )*1e2;
m_Amp1=median(Amp1,"omitnan" ); m_Amp2=median(Amp2,"omitnan" ); 
m_Amp3=median(Amp3,"omitnan" );
m_Cor1=median(Cor1,"omitnan" ); m_Cor2=median(Cor2,"omitnan" ); 
m_Cor3=median(Cor3,"omitnan" ); 
dir_current = mod(atan2d(U_QC, V_QC), 360);
m_spd=median(CSPD_QC,'omitnan')*1e2;
m_dir=median(dir_current,'omitnan');

total_pings = size(Data.mask_QC_3D, 1);
valid_pings = sum(Data.mask_QC_3D == 0, 1); 
m_coverage = (valid_pings ./ total_pings) * 100;

for k=1:srf_bins-1;
fprintf(fidlog,'\nBin %d : nominally %3.2f m - %3.2f m from sensor head\n\n', k, Dist2Instr_CellMidpoint(k)-Cell_Size,Dist2Instr_CellMidpoint(k)+Cell_Size);
fprintf(fidlog,'Data coverage [%%]                               : %4.1f%%\n', m_coverage(k));
fprintf(fidlog,'Median velocity u / v / w [cm/s]                : %4.1f  %4.1f  %4.1f\n',m_u(k), m_v(k), m_w(k));
fprintf(fidlog,'Median Amp Beam1 / Beam2 / Beam3 [db]           : %3.0f  %3.0f  %3.0f\n',m_Amp1(k), m_Amp2(k), m_Amp3(k));
fprintf(fidlog,'Median correlation Beam1 / Beam2 / Beam3 [%%]    : %3.0f  %3.0f  %3.0f\n',m_Cor1(k), m_Cor2(k), m_Cor3(k));
fprintf(fidlog,'Median speed [cm/s]                             : %4.1f\n',m_spd(k));
fprintf(fidlog,'Median direction [deg]                         : %5.2f\n',m_dir(k));
end

save_signature_to_nc(Data, Config, nCells, Dist2Instr_CellMidpoint, ...
    Nominal_CellDepth,flag_vals,flag_mean,operator,moor,info_adcp, outfile);

end
fprintf(fidlog, '\n==== END ENTRY  =====\n');
fclose(fidlog);
end

%% nested fuctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% looks for deployment and recovery period in pressure
function S = suggestTrimForShallowEdges_simple(y, ymin)
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

function save_signature_to_nc(Data, Config, nCells, Dist2Instr_CellMidpoint, Nominal_CellDepth,flag_vals,flag_mean,operator,moor,info_adcp, filename)
% SAVE_SIGNATURE_TO_NC Saves full Signature ADCP Data structure and Config attributes to NetCDF
%
% USAGE:
%   save_signature_to_nc(Data, Config, nCells, Dist2Instr, Nominal_Depth, 'filename.nc')

    if exist(filename, 'file'), delete(filename); end

    % 1. Determine Dimensions
    [nTime, nBins] = size(Data.Average_VelEast);
    
    % 2. Create Dimensions
    nccreate(filename, 'time', 'Dimensions', {'time', nTime});
    nccreate(filename, 'cell', 'Dimensions', {'cell', nBins});
    nccreate(filename, 'xyz',  'Dimensions', {'xyz', 3});
    nccreate(filename, 'beam_map', 'Dimensions', {'beam_map', 4});

    % 3. Write Data Variables and Assign Units
    all_fields = fieldnames(Data);

    for i = 1:numel(all_fields)
        f = all_fields{i};
        if strcmp(f, 'Average_Time'),continue;end
        val = Data.(f);
        [rows, cols] = size(val);

        % --- Determine Dimensions and Create Variable ---
        if rows == nTime && cols == 1
            nccreate(filename, f, 'Dimensions', {'time', nTime});
        elseif rows == nTime && cols == nBins
            nccreate(filename, f, 'Dimensions', {'time', nTime, 'cell', nBins});
        elseif rows == nTime && cols == 3
            nccreate(filename, f, 'Dimensions', {'time', nTime, 'xyz', 3});
        elseif rows == nTime && cols == 4
            nccreate(filename, f, 'Dimensions', {'time', nTime, 'beam_map', 4});
        else
            continue; % Skip if dimensions don't match known patterns
        end

        % --- Write Data ---
        ncwrite(filename, f, double(val));

        % --- Assign Units Attribute ---
        unit_str = '';
        if contains(f, 'Vel'), unit_str = 'm/s';
        elseif contains(f, 'Amp'), unit_str = 'dB';
        elseif contains(f, 'Cor'), unit_str = '%';
        elseif contains(f, 'Temp'), unit_str = 'degC';
        elseif (contains(f, 'Pressure') || contains(f, 'PressureSensor')) && ~contains(f, 'Temp')
            unit_str = 'dbar';
        elseif contains(f, 'Magnetometer') && ~contains(f, 'Temp'), unit_str = 'mG';
        elseif contains(f, 'Accelerometer'), unit_str = 'g';
        elseif contains(f, 'Energy'), unit_str = 'Joules';
        elseif contains(f, 'Heading') || contains(f, 'Pitch') || contains(f, 'Roll')
            unit_str = 'degrees';
            if contains(f, 'mask_QC')              
                ncwriteatt(filename, f, 'flag_values', int8(flag_vals));
                ncwriteatt(filename, f, 'flag_meanings', flag_mean);
            end
        elseif contains(f, 'Time'), unit_str = 'datenum';
        elseif contains(f, 'mask') || contains(f, 'Status') || contains(f, 'Error'), unit_str = 'flag';
        elseif contains(f, 'Soundspeed'), unit_str = 'm/s';
        elseif contains(f, 'Battery'), unit_str = 'V';
        end
        
        if ~isempty(unit_str)
            ncwriteatt(filename, f, 'units', unit_str);
        end
    end

    % 4. Write Coordinates (Dimensions)
    ncwrite(filename, 'time', Data.Average_Time);
    ncwriteatt(filename, 'time', 'units', 'datenum');

    ncwrite(filename, 'cell', double(nCells));
    ncwriteatt(filename, 'cell', 'long_name', 'Cell Index');

    nccreate(filename, 'Dist2Instr_CellMidpoint', 'Dimensions', {'cell', nBins});
    ncwrite(filename, 'Dist2Instr_CellMidpoint', double(Dist2Instr_CellMidpoint));
    ncwriteatt(filename, 'Dist2Instr_CellMidpoint', 'units', 'm');

    nccreate(filename, 'Nominal_CellDepth', 'Dimensions', {'cell', nBins});
    ncwrite(filename, 'Nominal_CellDepth', double(Nominal_CellDepth));
    ncwriteatt(filename, 'Nominal_CellDepth', 'units', 'm');

    % Latitude Variable
    nccreate(filename, 'latitude', 'Datatype', 'double');
    ncwrite(filename, 'latitude', double(info_adcp.lat));
    ncwriteatt(filename, 'latitude', 'units', 'degrees_north');
    ncwriteatt(filename, 'latitude', 'long_name', 'Latitude');
    
    % Longitude Variable
    nccreate(filename, 'longitude', 'Datatype', 'double');
    ncwrite(filename, 'longitude', double(info_adcp.lon));
    ncwriteatt(filename, 'longitude', 'units', 'degrees_east');
    ncwriteatt(filename, 'longitude', 'long_name', 'Longitude');

    % Final Global Metadata
    ncwriteatt(filename, '/', 'history', ['Created on ', datestr(now)]);
    ncwriteatt(filename, '/', 'Operator', operator);
    ncwriteatt(filename, '/', 'Mooring', moor);
    ncwriteatt(filename, '/', 'nominal_water_depth_m', double(info_adcp.wd));
    % 5. Write FULL Config structure as Global Attributes
    conf_fields = fieldnames(Config);
    for i = 1:numel(conf_fields)
        
        cf = conf_fields{i};
        cval = Config.(cf);
        
        if isempty(cval), continue; end
        
        % Special Case: Matrix Beam2xyz cannot be a global attribute
        if strcmp(cf, 'Average_Beam2xyz')
            nccreate(filename, 'Config_Average_Beam2xyz', 'Dimensions', {'row', 3, 'col', 3});
            ncwrite(filename, 'Config_Average_Beam2xyz', double(cval));
            ncwriteatt(filename, 'Config_Average_Beam2xyz', 'description', 'Beam to XYZ transformation matrix');
        elseif isnumeric(cval) && isscalar(cval)
            % Write numeric scalars directly
            ncwriteatt(filename, '/', cf, double(cval));
        else
            % Convert strings, booleans, and multi-element chars to char arrays
            ncwriteatt(filename, '/', cf, char(string(cval)));
        end
    end
    fprintf('NetCDF file "%s" created successfully with full attributes.\n', filename);
end

