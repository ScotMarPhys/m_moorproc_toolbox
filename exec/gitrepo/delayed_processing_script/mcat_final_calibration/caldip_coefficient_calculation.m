%   Calculate the different pre and post deployment calibration coefficient
%   from the cal_dip casts
%
%       by loh, 06/05/2016
%
% JC238 update - Lewis Drysdale, 2022

close all 

% caldip casts for ar304: 1, 2, 3, 5, 6, 45
% NEED to reprocess station 6 and 45, look for problem in detection of
% bottle stop (diff mcat-ctd = nan)


%========================================================================
% Calculation of calibration coefficient

p_insitucal.cruise           =  MOORPROC_G.cruise;
p_insitucal.cast = input('which cast number? ');
p_insitucal.depl_period = input('which deployment period (e.g. osnap6) ','s');

% ---- parameters ----------------------------------------------------
p_insitucal.sensorselect      = input('which CTD sensor? (1 or 2) ');
p_insitucal.sensor_id        = [332 337];  % MicroCAT ID range (in info.dat) %***
pd = moor_inoutpaths('cal_coef',p_insitucal.cast);
p_insitucal.datadir          = pd.datadir;
p_insitucal.coef_dir         = pd.coef_dir;
p_insitucal.apply_offset   = input('apply time offset between CTD and MC (y/n/i) ','s');
% if offset == 'y'/'n'/'i', time offset between CTD and MC 
                     % will / will not be applied / individual offsets
                     % applied  
% Define the time interval in which data are selected for each bottle stop.
% The time selection is centered on the the bottle stop time defined in the
% rosfile
if strcmp(p_insitucal.cruise,'pe399') & (p_insitucal.cast == 17 | p_insitucal.cast == 38)   
    p_insitucal.interval_move  = [-100 -50];
elseif strcmp(p_insitucal.cruise,'dy078') | strcmp(p_insitucal.cruise,'ar304') | strcmp(p_insitucal.cruise,'dy120')
    p_insitucal.interval_move  = [-600 600]*2; % bottlestop can be up to 20 mins for 02 sensors   
else
    p_insitucal.interval_move  = [-480 480];    
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
if strcmp(p_insitucal.cruise,'kn221') &   p_insitucal.cast == 4    
    p_insitucal.c_interval     = [-0.025 0.025]; % [-.06 .03]; % cond. dev.plot range
    p_insitucal.t_interval     = [-.015 .015]; % temp. dev.plot range 
    p_insitucal.dp_interval    = [-6 6];   % pres. dev. plot range 
    p_insitucal.p_interval       = [0  3200]; % pressure plot range
    p_insitucal.average_interval = [1000 3200];  % depth interval in which dc / dt are averaged 
                                     % (should be deep, to reduce error  from long adjustment 
                                     % times  required by older microcat conductivity sensors   
elseif strcmp(p_insitucal.cruise,'kn221') &   p_insitucal.cast == 5   
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-6 6];   
    p_insitucal.p_interval       = [0  2400];
    p_insitucal.average_interval = [1000 2400]; 
elseif strcmp(p_insitucal.cruise,'pe399') &   p_insitucal.cast == 4   
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2300];
    p_insitucal.average_interval = [1450 2300];      
elseif strcmp(p_insitucal.cruise,'pe399') &   p_insitucal.cast == 17   
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2800];
    p_insitucal.average_interval = [1000 2800];    
elseif strcmp(p_insitucal.cruise,'pe399') &   p_insitucal.cast == 3   
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2100];
    p_insitucal.average_interval = [1000 2100];     
elseif strcmp(p_insitucal.cruise,'pe399') &   p_insitucal.cast == 38  
    p_insitucal.c_interval     = [-0.4 0.08]; 
    p_insitucal.t_interval     = [-.3 .3];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2600];
    p_insitucal.average_interval = [1000 2600];   
elseif strcmp(p_insitucal.cruise,'pe400') & (p_insitucal.cast == 26  | p_insitucal.cast == 32)  
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-6 6];   
    p_insitucal.p_interval       = [0  3000];
    p_insitucal.average_interval = [1000 3000]; 
elseif strcmp(p_insitucal.cruise,'pe400') & p_insitucal.cast == 35  
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-6 6];   
    p_insitucal.p_interval       = [0  2200];
    p_insitucal.average_interval = [1000 2220];     
elseif strcmp(p_insitucal.cruise,'kn221-03') & p_insitucal.cast == 1  
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-6 6];   
    p_insitucal.p_interval       = [0  2200];
    p_insitucal.average_interval = [1000 2200];      
elseif strcmp(p_insitucal.cruise,'dy053') & p_insitucal.cast == 4 
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2400];
    p_insitucal.average_interval = [1000 2250];      
elseif strcmp(p_insitucal.cruise,'dy053') & p_insitucal.cast == 5 
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2400];
    p_insitucal.average_interval = [1000 2200];         
elseif strcmp(p_insitucal.cruise,'dy053') & p_insitucal.cast == 22 
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2950];
    p_insitucal.average_interval = [1000 2950]; 
elseif strcmp(p_insitucal.cruise,'dy078') & p_insitucal.cast == 2 
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2100];
    p_insitucal.average_interval = [300 800];      
elseif strcmp(p_insitucal.cruise,'dy078') & p_insitucal.cast == 3
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2150];
    p_insitucal.average_interval = [500 1000];% [1000 2150];    
elseif strcmp(p_insitucal.cruise,'dy078') & p_insitucal.cast == 27
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2150];
    p_insitucal.average_interval = [1500 2150];% [1000 2150];   
elseif strcmp(p_insitucal.cruise,'ar304') & p_insitucal.cast == 1
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2580];
    p_insitucal.average_interval = [1000 2580];% [1000 2150];  
elseif strcmp(p_insitucal.cruise,'ar304') & p_insitucal.cast == 2
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2580];
    p_insitucal.average_interval = [1050 2580];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'ar304') & p_insitucal.cast == 3
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  3000];
    p_insitucal.average_interval = [1050 3000];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'ar304') & p_insitucal.cast == 5
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  1850];
    p_insitucal.average_interval = [250 1850];% [1000 2150];   
elseif strcmp(p_insitucal.cruise,'ar304') & p_insitucal.cast == 6
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  1850];
    p_insitucal.average_interval = [1820 1850];% [1000 2150];    
elseif strcmp(p_insitucal.cruise,'ar304') & p_insitucal.cast == 45
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  3100];
    p_insitucal.average_interval = [1500 3100];% [1000 2150];  
elseif strcmp(p_insitucal.cruise,'dy120') & p_insitucal.cast == 1
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2580];
    p_insitucal.average_interval = [1000 2580];% [1000 2150];     
elseif strcmp(p_insitucal.cruise,'dy120') & p_insitucal.cast == 3
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  1850];
    p_insitucal.average_interval = [1500 1800];% [1000 2150]; 
elseif strcmp(p_insitucal.cruise,'dy120') & p_insitucal.cast == 4
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2050];
    p_insitucal.average_interval = [1500 2050];% [1000 2150];  
elseif strcmp(p_insitucal.cruise,'dy120') & p_insitucal.cast == 5
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  2050];
    p_insitucal.average_interval = [1500 2050];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'dy120') & p_insitucal.cast == 8
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  3000];
    p_insitucal.average_interval = [2000 3000];% [1000 2150]; 
elseif strcmp(p_insitucal.cruise,'dy120') & p_insitucal.cast == 9
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  3000];
    p_insitucal.average_interval = [2000 3000];% [1000 2150];  
elseif strcmp(p_insitucal.cruise,'dy120') & p_insitucal.cast == 10
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  3000];
    p_insitucal.average_interval = [1500 2950];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'jc238') & p_insitucal.cast == 19
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  1200];
    p_insitucal.average_interval = [1500 2950];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'dy181')  & p_insitucal.cast == 3
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0  1200];
    p_insitucal.average_interval = [1500 2950];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'dy181') & p_insitucal.cast == 4
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0 2500];
    p_insitucal.average_interval = [1000 1600];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'dy181') & p_insitucal.cast == 5
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0 2200];
    p_insitucal.average_interval = [1600 2200];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'dy181') & p_insitucal.cast == 46
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0 2000];
    p_insitucal.average_interval = [1200 1700];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'dy181') & p_insitucal.cast == 63
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0 1000];
    p_insitucal.average_interval = [0 1000];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'dy181') & p_insitucal.cast == 66
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0 2000];
    p_insitucal.average_interval = [100 2200];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'dy181') & p_insitucal.cast == 67
    p_insitucal.c_interval     = [-0.025 0.025]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0 3000];
    p_insitucal.average_interval = [1500 3100];% [1000 2150];
elseif strcmp(p_insitucal.cruise,'dy181') & p_insitucal.cast == 76
    p_insitucal.c_interval     = [-0.025 0.07]; 
    p_insitucal.t_interval     = [-.02 .02];
    p_insitucal.dp_interval    = [-10 10];   
    p_insitucal.p_interval       = [0 3000];
    p_insitucal.average_interval = [1200 3100];% [1000 2150];
else
    p_insitucal.c_interval     = [-0.025 0.025]; % plotting only
    p_insitucal.t_interval     = [-.02 .02]; % plotting only
    p_insitucal.dp_interval    = [-10 10];    % plotting only
    p_insitucal.p_interval       = [0  2000]; % plotting only
    p_insitucal.average_interval = [1000 2000];        
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
 
p_insitucal.bottlestop_tmin    = 60; %20 % minimum accepted separation between 2 bottlestops [sec] 
p_insitucal.bottlestop_dpmin   = 20;  % minimum accepted pressure difference between 2 bottlestops [dbar] 


p_insitucal.cnv_time_correction=  0;  % time correction for CTD .CNV file (rel. GMT) 
                          % (required for some casts of CD170 and RB0701 see below)
                          % times in cnv_cal do not match times of
                          % cnv,btl,ros or microcats. Normally 1 hour out.
                          
                          % !! Do not apply changes here, but below where input is defined 
p_insitucal.lat   = 58.0;    % latitude to convert depths into pressures. Difference between 
                         % pressures using 23.5N or 28N are <2dbar (upto
                         % depths of 6000m)
                          
                         
%==========================================================================
%-------- calibration of the coefficient ----------------------------------
insitu_cal_osnap2(p_insitucal)
%--------------------------------------------------------------------------

%=========================================================================
% Then edit the microcat_temp.csv, microcat_cond.csv and microcat_pres.csv
% files in osnap/data/moor/cal_coef/
%=========================================================================
