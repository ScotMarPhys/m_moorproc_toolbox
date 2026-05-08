% stage02_proc_S55(moor,varargin)
% basic preprocessing for stage 1 Signature 55 ADCP data (nc-format)
% Please make sure raw ADCP data in .mat format as exported by Signature
% Viewer were QC by calling stage01_read_qc_S55 before calling this
% function
%
% required inputs: moor - mooring name e.g. 'rhadcp_01_2020'
%                  
% optional inputs:     dataindir = varargin{1};
%                      infofile = varargin{2};
%                      logfile = varargin{3};
%                      outdir = varargin{4}; #output dir (fig, data, logs)
% functions called:    rodbload, gsw_z_from_p
% features
% 
% TODO
% - add option for MC pressure/speed of sound+profile_from_SCOTIA-clim for
%   depth correction during post-processing

function stage02_proc_S55(moor, varargin)

if nargin == 0
    help(mfilename);
    return;
end

global MOORPROC_G

if nargin == 1 && ~isempty(MOORPROC_G)
    operator = MOORPROC_G.operator;
    pd = moor_inoutpaths('adcp_S55',moor);
    dataindir = pd.stage1path;
    infofile = pd.infofile;
    logfile = pd.stage2log;
    outdir = pd.stage1path;
    infile_form = pd.stage1form;
    ouput_form = pd.stage2form;
else
    try
        operator = getenv('COMPUTERNAME');    
        dataindir = varargin{1};
        infofile = varargin{2};
        logfile = varargin{3};
        outdir = varargin{4};
        infile_form =  [moor '_%d_stage1.nc'];
        ouput_form = [moor '_%d_stage2.nc'];
    catch
        error('Not enough manual arguments provided. Expected 5 additional arguments after "moor".');
    end
end

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

% Directory and Logfile Validation
if ~exist(outdir, 'dir'), mkdir(outdir); end

[~, log_name, log_ext] = fileparts(logfile);
expected_log = [moor '_ADCP_stage2'];

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


% Load data
for i = 1:length(serial_nums)
fprintf('Processing sn %d',serial_nums(i));

outfile = fullfile(outdir,sprintf(ouput_form,serial_nums(i)));
[~, filename, ~] = fileparts(outfile);
infile = fullfile(dataindir,sprintf(infile_form,serial_nums(i)));
DS = s55_load_nc_as_struct(infile);

fprintf(fidlog,'infile : %s\n',infile);
fprintf(fidlog,'ADCP serial number  : %d\n\n',serial_nums(i));

%% Compass correction
% 1. Fix the Hard-Iron (Compass) error first
[DS, h_fig] = s55_hard_iron_compass_correction(DS, 0, true, true, true,outdir,filename);
err_min = min(DS.hard_iron_CCW_angle,[],'omitnan');
err_max = max(DS.hard_iron_CCW_angle,[],'omitnan');
t_lim = datestr(DS.calibration_window);

% Logfile output
deg = char(176); 
prombt = ['\n\n***Hard-Iron compass correction***\n' ...
    'Simple circle fitted to data from\n %s to %s.\n' ...
    'Horizontal velocity data corrected for error angle varying between %3.0f' deg ... 
    'and %3.0f' deg '\n'];
fprintf(1, prombt, t_lim(1,:),t_lim(2,:),err_min, err_max); 
fprintf(fidlog, prombt,t_lim(1,:),t_lim(2,:), err_min, err_max); 

% Save figure
set(h_fig,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print(h_fig,'-dpng',fullfile(outdir,[filename,'_f1_hard_iron_correction.png']));


% 2. Now rotate the corrected U/V to True North
magdev = info_adcp.magdev;
U_final = DS.U_hard_iron_corrected .* cosd(magdev) - DS.V_hard_iron_corrected .* sind(magdev);
V_final = DS.U_hard_iron_corrected .* sind(magdev) + DS.V_hard_iron_corrected .* cosd(magdev);

% Logfile output
deg = char(176); 
prombt = ['\n\n***Magnetic deviation***\n' ...
    'Horizontal velocity data corrected for magnetic deviation of %3.0f' deg '\n'];
fprintf(1, prombt, magdev); 
fprintf(fidlog, prombt,magdev); 

%% plot data to check rotation is reasonable
U_QC = DS.Average_VelEast; U_QC(DS.mask_QC_3D~=0)=NaN;
V_QC = DS.Average_VelNorth; V_QC(DS.mask_QC_3D~=0)=NaN;
U_HI = DS.U_hard_iron_corrected; U_HI(DS.mask_QC_3D~=0)=NaN;
V_HI = DS.V_hard_iron_corrected; V_HI(DS.mask_QC_3D~=0)=NaN;
U_final(DS.mask_QC_3D~=0)=NaN;
V_final(DS.mask_QC_3D~=0)=NaN;
y = DS.Nominal_CellDepth;

f2 = figure(2); clf;
% Change to 3 rows, 2 columns
tlo_main = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'tight'); 

% Reorder data so Row 1 = QC, Row 2 = HI, Row 3 = Final
data_list = {U_QC, V_QC, ...      % Row 1
             U_HI, V_HI, ...      % Row 2
             U_final, V_final};   % Row 3

titles = {'Uncorrected U', 'Uncorrected V', ...
          'Hard-Iron U', 'Hard-Iron V', ...
          'U Final', 'V Final'};

% High-Contrast Colormap
b_to_w = [linspace(0,1,11)', linspace(0,1,11)', ones(11,1)];
w_to_r = [ones(10,1), linspace(0.9,0,10)', linspace(0.9,0,10)'];
cmap = [[0 1 1]; b_to_w; w_to_r; [1 1 0]];

for k = 1:6
    % Nested layout remains the same
    tlo_sub = tiledlayout(tlo_main, 1, 10); 
    tlo_sub.Layout.Tile = k;
    
    % --- 1. The Main Contour Plot ---
    ax_cont = nexttile(tlo_sub, 1, [1 7]);
    im = imagesc(ax_cont, DS.time, y, data_list{k}');
    set(ax_cont, 'YDir', 'reverse', 'Color', [0.8 0.8 0.8]);
    set(ax_cont, 'YAxisLocation', 'right');
    set(im, 'AlphaData', ~isnan(data_list{k}'));
    datetick(ax_cont, 'x', 'mm-yy', 'keeplimits');
    ylim(ax_cont, [0, 1080]);
    colormap(ax_cont, cmap);
    clim(ax_cont, [-.52, .52]);
    title(ax_cont, titles{k});
    
    cb = colorbar(ax_cont, 'westoutside');
    cb.Ticks = [-.5, 0, .5];
    title(cb, 'm/s')
    
    % --- 2. The Slim Profile Plot ---
    ax_prof = nexttile(tlo_sub, 8, [1 3]);
    mean_vel = mean(data_list{k}, 1, 'omitnan');
    plot(ax_prof, mean_vel, y, 'k', 'LineWidth', 1);
    
    set(ax_prof, 'YDir', 'reverse', 'YTickLabel', [], 'Box', 'on');
    grid(ax_prof, 'on');
    ylim(ax_prof, [0, 1080]);
    xlim(ax_prof, [-0.1 0.1]); 
    xlabel(ax_prof, 'm/s');
    title(ax_prof,'mean')
end

% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',fullfile(outdir,[filename,'_stage2_f2_UV_hard_iron_mag_dev_correction.png']));
clear ax

%% Calculate depth matrix

% Use internal ADCP Sensors
pres = DS.Average_Pressure; 
pres(DS.mask_QC_1D ~= 0) = NaN; 
instr_depth = gsw_z_from_p(pres, DS.latitude);

% Prompt user: 'y' for MicroCAT, anything else for ADCP
user_choice = input('\n\nUse MicroCAT data for pressure/speed of sound correction? [y/n] (default n): ', 's');

use_MC = strcmpi(user_choice, 'y');

if use_MC
    try
        % 1. Determine folder
        indir_microcat = fullfile(fileparts(dataindir), 'microcat');
        if ~isfolder(indir_microcat)
            indir_microcat = input('Microcat path not found. Please enter path: ', 's');
        end
        if ~isfolder(indir_microcat), error('Path invalid'); end

        % 2. Find microcat closest to ADCP
        [mc_num, mc_sn] = find_microcat_index(info_adcp, serial_nums(i));

        mc_dlm_fn = sprintf('%s_%3.3d.microcat', moor, mc_num);
        mc_st2_fn = sprintf('%s_%4.4d.use', moor, mc_sn);

        if isfile(fullfile(indir_microcat, mc_dlm_fn))
            mc_infile = fullfile(indir_microcat, mc_dlm_fn);;
        elseif isfile(fullfile(indir_microcat, mc_st2_fn))
            mc_infile = fullfile(indir_microcat, mc_st2_fn);
        else
            error(sprintf('Neither %s nor %s file found.',mc_dlm_fn,mc_st2_fn));
        end
        
        % Updated prompt including the extension found
        prombt = sprintf('\n Using MicroCAT data (%s) & speed of sound correction.\n', ...
            strrep(mc_infile, '\', '\\'));

    catch ME
        % If any of the above fails, force use_MC to false and notify user
        warning('MicroCAT failed: %s. Falling back to ADCP sensors.', ME.message);
        use_MC = false;
    end
end

% This 'if' handles the original 'else' case AND the 'catch' fallback
if ~use_MC
    prombt = "\n Using internal ADCP sensors. No speed of sound correction.\n";
end


fprintf(fidlog,prombt);
fprintf(prombt);


if use_MC
        
    %% 1. Load and Interpolate MicroCAT data onto ADCP timeline
    MC = load_microcat_data(mc_infile);

    R_old = DS.Dist2Instr_CellMidpoint;
    c_old = DS.Average_Soundspeed;
    
    % Interpolate all variables to ADCP's Average_Time
    % 'pchip' or 'linear' are usually best; 'extrap' is risky but keeps the code running
    mc_t_interp  = interp1(MC.time, MC.t,  DS.time, 'linear', 'extrap');
    mc_p_interp  = interp1(MC.time, MC.p,  DS.time, 'linear', 'extrap');
    mc_sp_interp = interp1(MC.time, MC.SP, DS.time, 'linear', 'extrap');
    
    % 2. Flag ADCP data outside MicroCAT time range
    bad_in = (DS.time < min(MC.time, [], 'omitnan')) | ...
             (DS.time > max(MC.time, [], 'omitnan'));
    
    if any(bad_in)
        DS.mask_QC_1D(bad_in) = QC_BAD;
        DS.mask_QC_3D(bad_in, :) = QC_BAD;
    end

    mc_t_interp(DS.mask_QC_1D~=0)  = NaN;
    mc_p_interp(DS.mask_QC_1D~=0)  = NaN;
    mc_sp_interp(DS.mask_QC_1D~=0)  = NaN;
    pres(DS.mask_QC_1D ~= 0) = NaN; 
    c_old(DS.mask_QC_1D~=0)  = NaN;

    % apply offset to MC data
    offset = mean(pres-mc_p_interp,1,"omitmissing");
    
    % 3. Calculate Speed of Sound (SoS) using GSW
    % Use the interpolated values so the dimensions match DS.Average_Time
    [SA, ~] = gsw_SA_from_SP(mc_sp_interp, mc_p_interp, info_adcp.lon, info_adcp.lat);
    CT      = gsw_CT_from_t(SA, mc_t_interp, mc_p_interp);
    c_true   = gsw_sound_speed(SA, CT, mc_p_interp); 
    z_new = gsw_z_from_p(mc_p_interp+offset, info_adcp.lat); % MC Depth
    
    % scotia_fn = fullfile('C:\Users\sa07kb\OneDrive - SAMS\data\data_SCOTIA\SCOTIA_monthly_clim_V8.nc');
    scotia_fn = fullfile('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\OSNAPi_development_2023\Matlab\Climatology\SCOTIA_monthly_clim_V8.nc');

    scotia = load_scotia_at_location(scotia_fn, info_adcp.lon);
    scotia.z = gsw_z_from_p(scotia.pres, info_adcp.lat); % MC Depth
    
    %%
    scotia_sos_interp = interp1(-scotia.z, scotia.sos, ...
        DS.Nominal_CellDepth, 'linear', 'extrap');

    whos scotia_sos_interp c_old R_old
    R_true = correct_adcp_range(c_true, c_old, R_old);
    R_true_scotia = correct_adcp_range(scotia_sos_interp,1500*ones(1,12), R_old);
    
    %% 5. Correct Bin Range using the SoS Profile function
    % R_nominal is DS.Dist2Instr_CellMidpoint
    % c_inst is DS.Average_Soundspeed
    R_true = correct_adcp_range_by_profile(DS.Dist2Instr_CellMidpoint, ...
                                           DS.Average_Soundspeed, ...
                                           instr_depth, z_ref, c_ref, 'up');
else
    
    % Use nominal range (no profile correction)
    R_true = DS.Dist2Instr_CellMidpoint(:)';
end

% --- Final Calculation ---
% Since instr_depth and z_ref are negative (e.g. -500m) 
% and R_true is positive distance from head (e.g. +100m)
% For an UPWARD looking ADCP: Depth = instr_depth + R_true
cell_depth = instr_depth + R_true; 

fig = figure(1);
clf(fig); % Clear old plots from the previous ADCP loop
set(fig, 'Name', sprintf('Cell Depth - SN %d', serial_nums(i)), 'NumberTitle', 'off');
plot(DS.time,cell_depth')
datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits');
axis ij

%% % side lobe interference
%all angle in degrees
% gamma = 20; %slant angle
% alpha = [0, 120,240]; % beam 1 - aligned with x-axis, beam 2, beam 3
% pitch = Data.Average_Pitch; pitch(Data.mask_QC_1D~=0)=NaN; % y-axis rotation, in deg
% roll = Data.Average_Roll; roll(Data.mask_QC_1D~=0)=NaN; % x-axis rotation, in deg
% 
% % Pre-calculate trig terms (using 'd' versions for degrees)
% cp = cosd(pitch); sp = sind(pitch);
% cr = cosd(roll);  sr = sind(roll);
% cg = cosd(gamma); sg = sind(gamma);
% ca = cosd(alpha); sa = sind(alpha); % alpha is 1x3
% 
% % Angle to vertical for each beam (in degrees)
% beta = acosd(-sp .* sg .* ca + cp .* sr .* sg .* sa + cp .* cr .* cg);
% 
% % Find the minimum angle to vertical across all three beams for each timestamp
% max_beta = max(beta, [], 2);
% 
% % R contains buffer for upper range of cell
% R =R_true*cosd(max_beta)-Cell_Size; 
% % find(Dist2Instr_CellMidpoint <= R(i), 1, 'last') for each time step i
% R(Data.mask_QC_1D~=0)=info_adcp.wd;
% idx_valid = arrayfun(@(r) find(Dist2Instr_CellMidpoint <= r, 1, 'last'),...
%     R);
% idx_valid(Data.mask_QC_1D~=0)=0;
% R(Data.mask_QC_1D~=0)=NaN;

% plot it
% hLine = plot(time(Data.mask_QC_1D==0), y(idx_valid(Data.mask_QC_1D==0)), 'r', 'LineWidth', 0.1);
% hLine.DisplayName = 'Sidelobe interference';

% sidelobe contamination
% 
% prompt = sprintf([ ...
%   '\n\nSidelobe contamination for each time step varies between %d and %d ' ...
%   'out of %d bins.\nSee red line in figure 2.\n Would you like to flag ' ...
%   ' sidelobe contaminated bins for each time step as QC_BAD (%d). Y/N [Y]: '], ...
%   min(idx_valid(Data.mask_QC_1D==0))+1, max(idx_valid(Data.mask_QC_1D==0))+1, max(nCells), QC_BAD);
% 
% reply = input(prompt, 's');
% if isempty(reply)
%     reply = 'Y';
% end
% 
% while true
%     if strcmpi(reply, 'Y') || strcmpi(reply, 'YES')
%         mask = nCells > idx_valid(:);
%         Data.mask_QC_3D(mask) = QC_BAD;
% 
%         prombt = ['Sidelobe contamination for each time step varies between ' ...
%             ' %d and %d out of %d bins.\nBins contaminated by' ...
%             ' sidelobe interference flagged as QC_BAD (%d).\n'];
%         fprintf(fidlog, prombt, ...
%             min(idx_valid(Data.mask_QC_1D==0))+1, max(idx_valid(Data.mask_QC_1D==0))+1, max(nCells), QC_BAD);
%         break
% 
%     elseif strcmpi(reply, 'N') || strcmpi(reply, 'NO')
%         prombt = ['Sidelobe contamination for each time step varies between ' ...
%             'bin %d and %d out of %d bins.\nOperator chose NOT to flag ' ...
%             'bins contaminated by sidelobe interference.\n'];
%         fprintf(1, prombt, ...
%             min(idx_valid(Data.mask_QC_1D==0))+1, max(idx_valid(Data.mask_QC_1D==0))+1, max(nCells)); % to screen
%         fprintf(fidlog, prombt, ...
%             min(idx_valid(Data.mask_QC_1D==0))+1, max(idx_valid(Data.mask_QC_1D==0))+1, max(nCells));
%         break
% 
%     else
%         fprintf('Invalid entry. Enter Y or N (press Return for default).\n');
%         reply = input('Y/N [Y]: ', 's');   % re-prompt and read again
%         if isempty(reply)
%             reply = 'Y';
%         end
%     end
% end


end
fprintf(fidlog, '\n==== END ENTRY  =====\n');
fclose(fidlog);
end

% correct range for speed of sound
function R_true = correct_adcp_range(c_true, c_old, R_old)
% CORRECT_ADCP_RANGE Calculates true acoustic range using SOS ratio.
% Works if c_true is (Time x 1), (1 x Depth), or (Time x Depth).
%
% INPUTS:
%   c_true : [T x 1], [1 x D], or [T x D] Correct sound speed
%   c_old  : [T x 1] Sound speed used by instrument
%   R_old  : [D x 1] or [1 x D] Nominal bin ranges

    % 1. Standardize inputs to ensure correct orientation for broadcasting
    c_old = c_old(:);       % Force to [T x 1]
    R_old = R_old(:)';      % Force to [1 x D] (Row vector)
    
    % 2. Handle c_true orientation based on its size
    [rows, cols] = size(c_true);
    nT = length(c_old);
    nD = length(R_old);

    if isvector(c_true)
        if length(c_true) == nD
            c_true = c_true(:)'; % Force to [1 x D]
        else
            c_true = c_true(:);  % Force to [T x 1]
        end
    elseif rows == nD && cols == nT
        c_true = c_true';        % Flip [D x T] to [T x D]
    end

    % 3. Calculate using Implicit Expansion
    % MATLAB automatically expands (T x 1), (1 x D), or (T x D) 
    % to match the final result matrix [T x D].
    R_true = (c_true ./ c_old) .* R_old;

end

function scotia = load_scotia_at_location(scotia_fn, adcp_lon)
    % 1. Load longitude grid to find the nearest point
    lon_grid = ncread(scotia_fn, 'lon');
    
    % Handle potential 0-360 vs -180-180 differences
    adj_adcp_lon = adcp_lon;
    if any(lon_grid < 0) && adcp_lon > 180, adj_adcp_lon = adcp_lon - 360; end
    if any(lon_grid > 180) && adcp_lon < 0, adj_adcp_lon = adcp_lon + 360; end

    % Find index of the closest longitude
    [~, idx] = min(abs(lon_grid - adj_adcp_lon));
    scotia.lon_match = lon_grid(idx);
    
    % 2. Load 1D variables
    scotia.pres  = ncread(scotia_fn, 'pres');  % [251 x 1]
    scotia.month = ncread(scotia_fn, 'month'); % [12 x 1]
    
    % 3. Extract 3D slices for the specific longitude point (idx)
    % NC Dimensions are [pres, distance, month] -> [251, 652, 12]
    start_idx = [1, idx, 1];
    count_idx = [inf, 1, inf];
    
    % Use squeeze to get [Pressure x Month] matrices
    scotia.CT = squeeze(ncread(scotia_fn, 'CT', start_idx, count_idx));
    scotia.SA = squeeze(ncread(scotia_fn, 'SA', start_idx, count_idx));

    all_nan_rows = all(isnan(scotia.CT), 2); 
    
    scotia.pres  = scotia.pres(~all_nan_rows);
    scotia.CT    = scotia.CT(~all_nan_rows, :);
    scotia.SA    = scotia.SA(~all_nan_rows, :);
    scotia.sos = gsw_sound_speed(scotia.SA, scotia.CT, scotia.pres); 

    fprintf('Matched ADCP Lon %.3f to SCOTIA Lon %.3f (Index %d)\n', ...
            adcp_lon, scotia.lon_match, idx);
end

function [idx, closest_sn] = find_microcat_index(info_adcp,sn)
    % Find ADCP and Microcat indices based on ID ranges
    is_adcp = info_adcp.sn == sn;
    is_mc   = info_adcp.id >= 333 & info_adcp.id <= 337;
    
    % Get depth of the ADCP (assumes one ADCP; if multiple, uses the first)
    adcp_z = info_adcp.z(find(is_adcp, 1));
    
    % Get depths and serial numbers of all Microcats
    mc_zs = info_adcp.z(is_mc);
    mc_sns = info_adcp.sn(is_mc);
    
    % Find which Microcat is closest to the ADCP depth
    [~, closest_pos] = min(abs(mc_zs - adcp_z));
    target_z = mc_zs(closest_pos);
    closest_sn = mc_sns(closest_pos);
    
    % Count Microcats shallower than (or equal to) the target depth
    % This returns the "n-th" Microcat from the surface
    idx = sum(mc_zs <= target_z);
end