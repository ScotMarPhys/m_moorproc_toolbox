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
calp.cruise           =  MOORPROC_G.cruise;
calp.cast = input('which cast number? ');
calp.depl_period = input('which deployment period (e.g. osnap6) ','s');

% ---- parameters ----------------------------------------------------
calp.sensorselect      = input('which CTD sensor? (1 or 2) ');
calp.sensor_id        = [332 337];  % MicroCAT ID range (in info.dat) %***
pd = moor_inoutpaths('cal_coef',calp.cast);
calp.datadir          = pd.datadir;
calp.coef_dir         = pd.coef_dir;
calp.apply_offset   = input('apply time offset between CTD and MC (y/n/i) ','s');
% if offset == 'y'/'n'/'i', time offset between CTD and MC 
                     % will / will not be applied / individual offsets
                     % applied  
% Define the time interval in which data are selected for each bottle stop.
% The time selection is centered on the the bottle stop time defined in the
% rosfile
if strcmp(calp.cruise,'pe399') & (calp.cast == 17 | calp.cast == 38)   
    calp.interval_move  = [-100 -50];
elseif strcmp(calp.cruise,'dy078') | strcmp(calp.cruise,'ar304') | strcmp(calp.cruise,'dy120')
    calp.interval_move  = [-600 600]*2; % bottlestop can be up to 20 mins for 02 sensors
elseif strcmp(calp.cruise,'jc238')
    calp.interval_move  = [-300 300]; % [-600 600] for RHADCP_SJ
else
    calp.interval_move  = [-480 480];    
end
%[-220 -260] % originally: [-100 -50];  %[-100 -50] move0bottlestop averaging interval [begin end] [seconds]  
                        % to reach the best positioning of interval within the bottlestops  
                        % (should be deep, to reduce error  from long adjustment 
                        % times  required by older microcat conductivity sensors    
calp.ctd_latestart_offset = 0; % [s]; time adjustment for ctd in case that was started AFTER it went into the water,
                           % in which case the detection of time difference
                           % between the ctd and microcats will not work
                           % unless the parameter is adjusted (if used in
                           % conjunction with apply_offset = 'i' the same
                           % offset has to be added to interval_move)
calp.dp_tol           = 200; % if rms(dp) > dp_tol then try to replace MicroCAT by CTD pressure

% --- Graphics parameters (only for display purposes -------
if strcmp(calp.cruise,'kn221') &   calp.cast == 4    
    calp.c_interval     = [-0.025 0.025]; % [-.06 .03]; % cond. dev.plot range
    calp.t_interval     = [-.015 .015]; % temp. dev.plot range 
    calp.dp_interval    = [-6 6];   % pres. dev. plot range 
    calp.p_interval       = [0  3200]; % pressure plot range
    calp.average_interval = [1000 3200];  % depth interval in which dc / dt are averaged 
                                     % (should be deep, to reduce error  from long adjustment 
                                     % times  required by older microcat conductivity sensors   
elseif strcmp(calp.cruise,'kn221') &   calp.cast == 5   
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-6 6];   
    calp.p_interval       = [0  2400];
    calp.average_interval = [1000 2400]; 
elseif strcmp(calp.cruise,'pe399') &   calp.cast == 4   
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2300];
    calp.average_interval = [1450 2300];      
elseif strcmp(calp.cruise,'pe399') &   calp.cast == 17   
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2800];
    calp.average_interval = [1000 2800];    
elseif strcmp(calp.cruise,'pe399') &   calp.cast == 3   
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2100];
    calp.average_interval = [1000 2100];     
elseif strcmp(calp.cruise,'pe399') &   calp.cast == 38  
    calp.c_interval     = [-0.4 0.08]; 
    calp.t_interval     = [-.3 .3];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2600];
    calp.average_interval = [1000 2600];   
elseif strcmp(calp.cruise,'pe400') & (calp.cast == 26  | calp.cast == 32)  
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-6 6];   
    calp.p_interval       = [0  3000];
    calp.average_interval = [1000 3000]; 
elseif strcmp(calp.cruise,'pe400') & calp.cast == 35  
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-6 6];   
    calp.p_interval       = [0  2200];
    calp.average_interval = [1000 2220];     
elseif strcmp(calp.cruise,'kn221-03') & calp.cast == 1  
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-6 6];   
    calp.p_interval       = [0  2200];
    calp.average_interval = [1000 2200];      
elseif strcmp(calp.cruise,'dy053') & calp.cast == 4 
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2400];
    calp.average_interval = [1000 2250];      
elseif strcmp(calp.cruise,'dy053') & calp.cast == 5 
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2400];
    calp.average_interval = [1000 2200];         
elseif strcmp(calp.cruise,'dy053') & calp.cast == 22 
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2950];
    calp.average_interval = [1000 2950]; 
elseif strcmp(calp.cruise,'dy078') & calp.cast == 2 
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2100];
    calp.average_interval = [300 800];      
elseif strcmp(calp.cruise,'dy078') & calp.cast == 3
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2150];
    calp.average_interval = [500 1000];% [1000 2150];    
elseif strcmp(calp.cruise,'dy078') & calp.cast == 27
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2150];
    calp.average_interval = [1500 2150];% [1000 2150];   
elseif strcmp(calp.cruise,'ar304') & calp.cast == 1
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2580];
    calp.average_interval = [1000 2580];% [1000 2150];  
elseif strcmp(calp.cruise,'ar304') & calp.cast == 2
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2580];
    calp.average_interval = [1050 2580];% [1000 2150];
elseif strcmp(calp.cruise,'ar304') & calp.cast == 3
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  3000];
    calp.average_interval = [1050 3000];% [1000 2150];
elseif strcmp(calp.cruise,'ar304') & calp.cast == 5
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  1850];
    calp.average_interval = [250 1850];% [1000 2150];   
elseif strcmp(calp.cruise,'ar304') & calp.cast == 6
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  1850];
    calp.average_interval = [1820 1850];% [1000 2150];    
elseif strcmp(calp.cruise,'ar304') & calp.cast == 45
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  3100];
    calp.average_interval = [1500 3100];% [1000 2150];  
elseif strcmp(calp.cruise,'dy120') & calp.cast == 1
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2580];
    calp.average_interval = [1000 2580];% [1000 2150];     
elseif strcmp(calp.cruise,'dy120') & calp.cast == 3
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  1850];
    calp.average_interval = [1500 1800];% [1000 2150]; 
elseif strcmp(calp.cruise,'dy120') & calp.cast == 4
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2050];
    calp.average_interval = [1500 2050];% [1000 2150];  
elseif strcmp(calp.cruise,'dy120') & calp.cast == 5
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2050];
    calp.average_interval = [1500 2050];% [1000 2150];
elseif strcmp(calp.cruise,'dy120') & calp.cast == 8
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  3000];
    calp.average_interval = [2000 3000];% [1000 2150]; 
elseif strcmp(calp.cruise,'dy120') & calp.cast == 9
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  3000];
    calp.average_interval = [2000 3000];% [1000 2150];  
elseif strcmp(calp.cruise,'dy120') & calp.cast == 10
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  3000];
    calp.average_interval = [1500 2950];% [1000 2150];
elseif strcmp(calp.cruise,'jc238') & calp.cast == 1
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2580];
    calp.average_interval = [0 750];% [1000 2150];   
elseif strcmp(calp.cruise,'jc238') & calp.cast == 3
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  1850]; % [0  2200] for RHADCP_SJ
    calp.average_interval = [1500 1800];% [1200 2200] for RHADCP_SJ
elseif strcmp(calp.cruise,'jc238') & calp.cast == 4
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  2050];
    calp.average_interval = [1500 2050];% [1000 2150];  
elseif strcmp(calp.cruise,'jc238') & calp.cast ==19
    calp.c_interval     = [-0.03 0.03]; % [-0.025 0.025] for RHADCP_SJ
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [500 1090]; %[0  1200] for RHADCP_SJ
    calp.average_interval = [800 1100];% [1500 2950] for RHADCP_SJ
elseif strcmp(calp.cruise,'jc238') & calp.cast == 33
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  3000];
    calp.average_interval = [2000 3000];% [1000 2150]; 
elseif strcmp(calp.cruise,'jc238') & calp.cast == 37
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  3000];
    calp.average_interval = [2000 3000];% [1000 2150];  
elseif strcmp(calp.cruise,'jc238') & calp.cast == 38
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  3000];
    calp.average_interval = [1500 2950];% [1000 2150];
elseif strcmp(calp.cruise,'jc238') & calp.cast == 43
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  3000];
    calp.average_interval = [1500 2950];% [1000 2150];
elseif strcmp(calp.cruise,'dy181')  & calp.cast == 3
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0  1200];
    calp.average_interval = [1500 2950];% [1000 2150];
elseif strcmp(calp.cruise,'dy181') & calp.cast == 4
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0 2500];
    calp.average_interval = [1000 1600];% [1000 2150];
elseif strcmp(calp.cruise,'dy181') & calp.cast == 5
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0 2200];
    calp.average_interval = [1600 2200];% [1000 2150];
elseif strcmp(calp.cruise,'dy181') & calp.cast == 46
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0 2100];
    calp.average_interval = [950 2150]; %[1200 1700];% [1000 2150];
elseif strcmp(calp.cruise,'dy181') & calp.cast == 63
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0 1000];
    calp.average_interval = [0 1000];% [1000 2150];
elseif strcmp(calp.cruise,'dy181') & calp.cast == 66
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0 3000];
    calp.average_interval = [1500 2940];% [1000 2150]; 
elseif strcmp(calp.cruise,'dy181') & calp.cast == 67
    calp.c_interval     = [-0.025 0.025]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0 3000];
    calp.average_interval = [1500 3100];% [1000 2150];
elseif strcmp(calp.cruise,'dy181') & calp.cast == 76
    calp.c_interval     = [-0.025 0.07]; 
    calp.t_interval     = [-.02 .02];
    calp.dp_interval    = [-10 10];   
    calp.p_interval       = [0 3000];
    calp.average_interval = [1500 3100];% [1000 2150];
else
    calp.c_interval     = [-0.025 0.025]; % plotting only
    calp.t_interval     = [-.02 .02]; % plotting only
    calp.dp_interval    = [-10 10];    % plotting only
    calp.p_interval       = [0  2000]; % plotting only
    calp.average_interval = [1000 2000];        
end

    
% % --- futhers parameters (unlikely to be changed) -------------------
% if strcmp(calp.cruise,'kn221-03')  & calp.cast == 1 
%     calp.bottlestop_tmin    = 30;  % minimum accepted duration of bottlestop [sec] 
%     calp.bottlestop_average = 2*60;  % default data average interval in seconds (during bottlestop)
% else
%     calp.bottlestop_tmin    = 20;  % minimum accepted duration of bottlestop [sec] 
%     calp.bottlestop_average = 5*60;  % default data average interval in seconds (during bottlestop)
% end                          

calp.bottlestop_tmin    = 60*4; % minimum accepted duration of bottlestop [sec] 
calp.bottlestop_average = 5*60; % default data average interval in seconds (during bottlestop)

calp.cond_threshold     = 30;  % [mS/cm] threshold to determine when instruments enter the water,
                          % use same threshold as in ctd_impact.m  
calp.impact_var         = 'c'; % conductivity as variable in ctd_impact.m
 
calp.bottlestop_tmin    = 60; %20 % minimum accepted separation between 2 bottlestops [sec] 
calp.bottlestop_dpmin   = 20;  % minimum accepted pressure difference between 2 bottlestops [dbar] 


calp.cnv_time_correction=  0;  % time correction for CTD .CNV file (rel. GMT) 
                          % (required for some casts of CD170 and RB0701 see below)
                          % times in cnv_cal do not match times of
                          % cnv,btl,ros or microcats. Normally 1 hour out.
                          
                          % !! Do not apply changes here, but below where input is defined 
calp.lat   = 58.0;    % latitude to convert depths into pressures. Difference between 
                         % pressures using 23.5N or 28N are <2dbar (upto
                         % depths of 6000m)
                          
                         
%==========================================================================
%-------- calibration of the coefficient ----------------------------------
insitu_cal_osnap2(calp)
%--------------------------------------------------------------------------

%=========================================================================
% Then edit the microcat_temp.csv, microcat_cond.csv and microcat_pres.csv
% files in osnap/data/moor/cal_coef/
%=========================================================================
