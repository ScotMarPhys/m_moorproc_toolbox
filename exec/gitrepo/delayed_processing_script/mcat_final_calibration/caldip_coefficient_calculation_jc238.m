%   Calculate the different pre and post deployment calibration coefficient
%   from the cal_dip casts
%
%       by loh, 06/05/2016
%
% JC238 update - Lewis Drysdale, 2022

close all 
%========================================================================
% Calculation of calibration coefficient
% For cruise jc238: cast= [1,3,4,19,33,37,38,43]
p_insitucal.cruise           =  'jc238' ;
p_insitucal.cast             = 19;   
p_insitucal.depl_period      = 'osnap6';

% ---- parameters ----------------------------------------------------
p_insitucal.sensorselec      = 1;
p_insitucal.sensor_id        = [332 337];  % MicroCAT ID range (in info.dat)
p_insitucal.basedir          = pathosnap; % base directory for osnap mooring
p_insitucal.datadir          = [p_insitucal.basedir filesep 'data']; % data directory 
p_insitucal.coef_dir         = [pathgit '/data/moor/cal_coef/']; 
p_insitucal.apply_offset   = 'n'; % if offset == 'y'/'n'/'i', time offset between CTD and MC 
                     % will / will not be applied / individual offsets
                     % applied  
% Define the time interval in which data are selected for each bottle stop. The time selection is centered on the the bottle stop time defined in the rosfile                      
if strcmp(p_insitucal.cruise,'pe399') & (p_insitucal.cast == 17 | p_insitucal.cast == 38)   
    p_insitucal.interval_move  = [-100 -50];
elseif strcmp(p_insitucal.cruise,'dy078') | strcmp(p_insitucal.cruise,'ar304') | strcmp(p_insitucal.cruise,'dy120')
    p_insitucal.interval_move  = [-600 600]*2; % bottlestop can be up to 20 mins for 02 sensors   
elseif strcmp(p_insitucal.cruise,'jc238') % try 5 mins either side of bot time
    p_insitucal.interval_move  = [-300 300]; % bottlestop can be up to 20 mins for 02 sensors 
else
    % this is completely wrong, I think! Would always through up an error
    p_insitucal.interval_move  = [-280 -320];    
end

%[-220 -260] % originally: [-100 -50];  %[-100 -50] move0bottlestop averaging interval [begin end] [seconds]  
                        % to reach the best positioning of interval within the bottlestops  
                        % (should be deep, to reduce error  from long adjustment 
                        % times  required by older microcat conductivity sensors    
p_insitucal.ctd_latestart_offset = 0; % [s]; time adjustment for ctd in case that was started AFTER it went into the water,
                           % in which case the detection of time difference
                           % between the ctd and microcats will not work
                           % unless the parameter is adjusted (if used in
                           % conjunction with apply_offset = 'i' the same
                           % offset has to be added to interval_move)
p_insitucal.dp_tol           = 200; % if rms(dp) > dp_tol then try to replace MicroCAT by CTD pressure

% --- Graphics parameters (only for display purposes -------
%%%%%%%%%%%%%%%%%%%%%%%%%%% THIS IS NOT TRUE ^^^^^^ %%%%%%%%%%
% The parameter average_interval below directly affects the calculation of 
% calibration coeffecients

if strcmp(p_insitucal.cruise,'jc238') & p_insitucal.cast == 1
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2580];
    p_insitucal.average_interval = [1000 2580];% [1000 2150];     
elseif strcmp(p_insitucal.cruise,'jc238') & p_insitucal.cast == 3
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  1850];
    p_insitucal.average_interval = [1500 1800];% [1000 2150]; 
elseif strcmp(p_insitucal.cruise,'jc238') & p_insitucal.cast == 4
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2050];
    p_insitucal.average_interval = [1500 2050];% [1000 2150];  
elseif strcmp(p_insitucal.cruise,'jc238') & p_insitucal.cast ==19
    p_insitucal.c_interval     = [-0.03 0.03]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [500 1090];
    p_insitucal.average_interval = [800 1100];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'jc238') & p_insitucal.cast == 33
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  3000];
    p_insitucal.average_interval = [2000 3000];% [1000 2150]; 
elseif strcmp(p_insitucal.cruise,'jc238') & p_insitucal.cast == 37
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  3000];
    p_insitucal.average_interval = [2000 3000];% [1000 2150];  
elseif strcmp(p_insitucal.cruise,'jc238') & p_insitucal.cast == 38
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  3000];
    p_insitucal.average_interval = [1500 2950];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'jc238') & p_insitucal.cast == 43
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  3000];
    p_insitucal.average_interval = [1500 2950];% [1000 2150];
else
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  4000];
    p_insitucal.average_interval = [1000 4000];        
end
    
% % --- futhers parameters (unlikely to be changed) -------------------
% if strcmp(p_insitucal.cruise,'kn221-03')  & p_insitucal.cast == 1 
%     p_insitucal.bottlestop_tmin    = 30;  % minimum accepted duration of bottlestop [sec] 
%     p_insitucal.bottlestop_average = 2*60;  % default data average interval in seconds (during bottlestop)
% else
%     p_insitucal.bottlestop_tmin    = 20;  % minimum accepted duration of bottlestop [sec] 
%     p_insitucal.bottlestop_average = 5*60;  % default data average interval in seconds (during bottlestop)
% end                          

p_insitucal.bottlestop_tmin    = 60*4; % minimum accepted duration of bottlestop [sec] 
p_insitucal.bottlestop_average = 5*60; % default data average interval in seconds (during bottlestop)

p_insitucal.cond_threshold     = 30;  % [mS/cm] threshold to determine when instruments enter the water,
                          % use same threshold as in ctd_impact.m  
p_insitucal.impact_var         = 'c'; % conductivity as variable in ctd_impact.m
 
p_insitucal.bottlestop_tmin    = 60; %20 % minimum accepted range between 2 bottlestops [sec] 
p_insitucal.bottlestop_dpmin   = 30;  % minimum accepted pressure difference of bottlestops [dbar] 


p_insitucal.cnv_time_correction=  0;  % time correction for CTD .CNV file (rel. GMT) 
                          % (required for some casts of CD170 and RB0701 see below)
                          % times in cnv_cal do not match times of
                          % cnv,btl,ros or microcats. Normally 1 hour out.
                          
                          % !! Do not apply changes here, but below where input is defined 
p_insitucal.lat   = 58.0    % latitude to convert depths into pressures. Difference between 
                         % pressures using 23.5N or 28N are <2dbar (upto
                         % depths of 6000m)
                          
                         
%==========================================================================
%-------- calibration of the coefficient ----------------------------------
insitu_cal_osnap2
%--------------------------------------------------------------------------

%=========================================================================
% Then edit the microcat_temp.csv, microcat_cond.csv and microcat_pres.csv
% files in osnap/data/moor/cal_coef/
%=========================================================================
