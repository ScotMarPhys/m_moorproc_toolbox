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
[DS, h_fig] = s55_hard_iron_compass_correction(DS, 0, true, true, true);
err_min = min(DS.hard_iron_CCW_angle,[],'omitnan');
err_max = max(DS.hard_iron_CCW_angle,[],'omitnan');
t_lim = datestr(DS.calibration_window);

% Logfile output
deg = char(176); 
prombt = ['\n\n***Hard-Iron compass correction***\n' ...
    'Simple circle fitted to data from\n %s to %s.\n' ...
    'Velocity data corrected for error angle varying between %3.0f' deg ... 
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
    'Velocity data corrected for magnetic deviation of %3.0f' deg '\n'];
fprintf(1, prombt, t_lim(1),t_lim(2),err_min, err_max); 
fprintf(fidlog, prombt,t_lim(1),t_lim(2), err_min, err_max); 

%% Calculate depth matrix

% Prompt user: 'y' for MicroCAT, anything else for ADCP
user_choice = input('\n\nUse MicroCAT data for \nPressure/SoS correction? [y/n] (default n): ', 's');

if strcmpi(user_choice, 'y')
    use_MC = true;
    prombt = ['\n\nUsing MicroCAT sensors data and apply speed of sound\n' ...
        ' correction for depth calculation.\n'];
else
    use_MC = false;
    prombt = ['\n\nUsing internal ADCP pressure sensors and \nno speed of sound' ...
        ' correction for depth calculation.\n'];
    
end

fprintf(fidlog,prombt);
fprintf(prombt);

% Check if use_MC is active and intercept if it's not operational
if use_MC
    warning(['The option to use MicroCat data is currently not operational. ' ...
             'Falling back to use internal ADCP sensors for pressure.']);
    use_MC = false; % Force the flag to false
end

if use_MC
    % 1. Load MicroCAT (MC) data placeholder
    % pd = moor_inoutpaths('microcat',moor)
      % case 'microcat'
      %   moor = loc;
      %   pd.rawpath = fullfile(mg.moordatadir, 'raw', mg.cruise, 'microcat');
      %   pd.infofile = fullfile(mg.moordatadir, 'proc', moor, [moor 'info.dat']);
      %   pd.stage1path = fullfile(mg.moordatadir, 'proc', moor, 'microcat');
      %   pd.stage1form = [moor '_%4.4d.raw'];
      %   pd.stage1log = fullfile(pd.stage1path,'stage1_log');
      %   pd.stage2path = fullfile(mg.moordatadir, 'proc', moor, 'microcat');
      %   pd.stage2form = [moor '_%4.4d.use'];
      %   pd.stage2log = fullfile(pd.stage2path, ['stage2_log_' moor,'.log']);
      %   pd.stage2figpath = fullfile(mg.reportdir, 'figs');
    % mooringpath  = [pathosnap '/data/moor/proc']; use pd for mooringpath
    % mc_path = [mooringpath,':',moor,':microcat:[',num2str(mc_ind),']'];  
    % [yy_mc,mm,dd,hh,t,c,p,sn_mc,depth_mc] = ...
    %     rodbload(mc_path,'yy:mm:dd:hh:t:c:p:SerialNumber:Instrdepth');

    % 2. Interpolate MC Pressure onto ADCP Time
    % interp1(source_time, source_data, target_time)
    pres = interp1(MC.time, MC.pres, DS.Average_Time, 'linear', 'extrap');
    
    % 3. Calculate Speed of Sound (SoS) profile from MC T and S
    % If you have a string of MicroCATs, c_ref and z_ref would be vectors.
    % If you only have one, it acts as a constant offset.
    c_ref = gsw_speed_of_sound(MC.salt, MC.temp, MC.pres); % True SoS
    z_ref = gsw_z_from_p(MC.pres, DS.latitude);           % MC Depth
    
    % 4. Calculate Instrument Depth (Negative for GSW consistency)
    instr_depth = gsw_z_from_p(pres, DS.latitude);
    
    % 5. Correct Bin Range using the SoS Profile function
    % R_nominal is DS.Dist2Instr_CellMidpoint
    % c_inst is DS.Average_Soundspeed
    R_true = correct_adcp_range_by_profile(DS.Dist2Instr_CellMidpoint, ...
                                           DS.Average_Soundspeed, ...
                                           instr_depth, z_ref, c_ref, 'up');
else
    % Use internal ADCP Sensors
    pres = DS.Average_Pressure; 
    pres(DS.mask_QC_1D ~= 0) = NaN; 
    instr_depth = gsw_z_from_p(pres, DS.latitude);
    
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
function [R_true] = correct_adcp_range_by_profile(R_nominal, c_inst, z_instr, z_ref, c_ref, orientation)
% CORRECT_ADCP_RANGE_BY_PROFILE Calculates true acoustic range using a SS profile.
%
% USAGE:
%   R_true = correct_adcp_range_by_profile(R_nominal, c_inst, z_instr, z_ref, c_ref, 'up')
%
% INPUTS:
%   R_nominal   : [1 x Cells] or [Cells x 1]Nominal distance from head to bin (positive m)
%   c_inst      : [Time x 1] Speed of sound used by ADCP (m/s)
%   z_instr     : [Time x 1] Depth of ADCP transducer (negative m, e.g., -500)
%   z_ref       : [M x 1] Reference depths for the SS profile (negative m)
%   c_ref       : [M x 1] Reference speed of sound values (m/s)
%   orientation : 'up' or 'down' (string)
%
% OUTPUTS:
%   R_true      : [Time x Cells] Physically corrected range from head (positive m)

    if ~isvector(R_nominal)
        error('R_nominal must be a vector (1D array).');
    end
    R_nominal = R_nominal(:)'; 

    nTime = length(c_inst);
    nCells = length(R_nominal);
    R_true = zeros(nTime, nCells);

    % 1. Convert nominal range to travel time (one-way)
    % This is the "Time Gate" the ADCP used to define each bin
    T_nominal = R_nominal ./ c_inst; 

    % 2. Iteratively integrate the profile for each time step
    for t = 1:nTime
        current_R = 0; % Distance from head starts at 0
        
        for k = 1:nCells
            % Time increment for this bin segment
            if k == 1
                dt = T_nominal(t, k);
            else
                dt = T_nominal(t, k) - T_nominal(t, k-1);
            end
            
            % Determine the depth of this specific bin segment
            if strcmpi(orientation, 'up')
                % Moving toward the surface: -500 + 10 = -490m
                current_Z = z_instr(t) + current_R; 
            else
                % Moving toward the bottom: -500 - 10 = -510m
                current_Z = z_instr(t) - current_R;
            end
            
            % Find local Speed of Sound at this segment's depth
            c_local = interp1(z_ref, c_ref, current_Z, 'linear', 'extrap');
            
            % Update "True" range
            current_R = current_R + (c_local * dt);
            R_true(t, k) = current_R;
        end
    end
end