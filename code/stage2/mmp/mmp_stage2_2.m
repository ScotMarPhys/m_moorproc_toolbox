% mmp_stage2.m
%
% Processes mmp data:  Conductivity Temperature Pressure, assumed to 
% be in rapid data format.  
%
% The user must specify an input directory where the file 
% mmp_ctd.txt is stored, by way of the mooring identifier. The 
% output directory is the same. The syntax for the file structure is 
% such that, if properly set up on each machine, only a prefix must be 
% changed to move from UNIX to Windows and back.  Forward slash
% was no problem on Windows XP.
%
% Data are read from a text file, mmp_ctd.txt, with the
% first 6 columns being timestamp (YYYY MM DD HH MM SS), then C, T, P, 
% and profile number.  All data from mmp_ctd.txt are held
% in a 3-D array alldata(maxrec,10,lastprof+1);     
% Dimension 1: max maxrec records per profile; 
% Dimension 2: 6 variables for time; c,t,p; profile number
% Dimension 3: Number of profiles
%
% MATLAB allows dynamic memory allocation, so one can go beyond 
% the size set initially.  But that makes the program extremely slow.  
% Therefore, one should initialise so large that all is covered.
%
% Some very basic extraction information is stored in a log file.
%
% mmp_stage2 is set up to work with arbitrary  start and end 
% profiles (firstprof and lastprof, respectively), to be set in
% the m-file.  Operator could not be bothered setting this up
% interactively.  Notice that the first mmp profile, the initial 
% downward profile, has number zero.  Since matlab doesn't like 
% index zero, this had to be accounted for in the loops through 
% profiles. But output is strictly according to original
% profile number. 
% 
% 16 - 17 & 26 March 2006 
% Jochem Marotzke, MPI-M (jym)

% HISTORY
% Version Date    Who     Changes
% _1      16,17,26/3/06 jym see introductory comments
% _2      30/5/06 scu     To speed programme now reads mmp_ctd_matlab.mat
%                         written by stage 1 mmp_read.m. This file is the
%                         3D file described above. Data are first decimated
%                         by a factor 5, then interpolated to 2db grid. The
%                         2db file is saved and programme ends with
%                         invitaiton to run plotting programme.

   
  format compact
  
  close all,  clear,    clear all,     
  moor     = 'eb2_4_200562';                
  operator = 'scu';
  disp(['Operator: ',operator])
  disp(['Start processing [YYYY MM DD HH MM SS]: ',num2str(fix(clock))])
% addpath(genpath ('C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data\exec\moor'))
% addpath(genpath ('C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data\exec\mfiles'))
% needs to be automated with startup

% Set prefix for Windows or UNIX:

%  prefix = 'C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data'  %jym's PC
%  prefix = '/local1/rapid/data'         % UNIX
%    prefix = '/data32/rapid/data'; % D304
     prefix = '/local/users/pstar/Data/rpdmoc/rapid/data' % D334 test
  inpath = [prefix,'/moor/proc/',moor,'/mmp/'];
  outpath = [prefix,'/moor/proc/',moor,'/mmp/'];
  infofile = [prefix,'/moor/proc/',moor,'/',moor,'info.dat']; 
  infile = [inpath,'mmp_ctd_matlab.mat'];
  outsave_2db = [outpath,'mmp_ctd_proc_2db.mat'];
  outdiary = [outpath,'mmp_stage2_log.txt'];
  diary (outdiary);
  disp(['Loading : ',infile])
  load(infile); 

  disp(['Decimating data'])
  
  firstprof = 0*7;     % First profile. NB: Could be profile zero!
  %lastprof = 156;      % eb2_4_200562
  lastprof = 156;
  %lastprof = 174;     % eb2_3_200516 
  maxrec = 14000;
  disp(['Processing profile ',num2str(firstprof),' to ', ...
      num2str(lastprof)])
  lprof = lastprof - firstprof +1; 
  cond = NaN(maxrec,lastprof+1);  % Initialise 2D matrix
  temp = NaN(maxrec,lastprof+1);
  press = NaN(maxrec,lastprof+1);
  
    for profnumber = firstprof : lastprof         % Loop through profiles

      prof = profnumber+1;
      disp(['Profile : ',num2str(profnumber)])
% Extract data from 3D array into a single vector to remove NANs
% for decimate function
        bad_data=ones(14000,1);
        i = find(isnan(alldata(:,7,prof)));bad_data(i)=NaN;
        i = find(isnan(alldata(:,8,prof)));bad_data(i)=NaN;
        i = find(isnan(alldata(:,9,prof)));bad_data(i)=NaN;
        i = find(~isnan(bad_data));
        a = decimate(alldata(i,7,prof),5,'FIR');
        cond(1:length(a),prof)=a;
        a = decimate(alldata(i,8,prof),5,'FIR');
        temp(1:length(a),prof)=a;
        a = decimate(alldata(i,9,prof),5,'FIR');
        press(1:length(a),prof)=a;
        runtime_r(:,prof) = datenum(alldata(:,1:6,prof)) - ...
          datenum('01-Jan-2006', 'dd-mmm-yyyy');
% find record number of pmax in each profile
        i = find(alldata(:,9,prof) == max(alldata(:,9,prof)));
        i = i(end); % sometimes more than one value of pmax - take the last
        profile_time(1:7,prof)=[runtime_r(i,prof) ...
            alldata(i,1:6,prof)];
    end
  
 disp(['... done it : ',num2str(fix(clock))])

% grid data to 2db
 disp(['Starting gridding to 2db']) 

% basic despiking limits
  c_lim = [20 70]; t_lim = [-2 40];p_lim = [0 6000]; s_lim = [30 40];

% Output array for 2db gridded data
  p2db = [80:2:2500]; % pressure grid
  lprof = lastprof - firstprof +1;      
  maxrec = length(p2db);          
  maxvar = 4;  % c t p prof_num
  alldata_2db = NaN(maxrec,maxvar,lastprof+1);

% interpolate
      for profnumber = firstprof : lastprof
        prof = profnumber+1;
        disp(['Profile : ',num2str(profnumber)])
        pt=press(:,prof);ct=cond(:,prof);tt=temp(:,prof);
        [pt,pII] = sort(pt);ct = ct(pII);tt = tt(pII);
        bad_data = ones(size(pt));
        i = find(pt<p_lim(1) & pt>p_lim(2)); bad_data(i) = NaN;
        i = find(tt<t_lim(1) & tt>t_lim(2)); bad_data(i) = NaN;
        i = find(ct<c_lim(1) & ct>c_lim(2)); bad_data(i) = NaN;
        i = find(isnan(pt)); bad_data(i) = NaN;
        i = find(isnan(tt)); bad_data(i) = NaN;
        i = find(isnan(ct)); bad_data(i) = NaN;
        i = find(~isnan(bad_data));
        ci  = interp1(pt(i),ct(i),p2db);
        ti  = interp1(pt(i),tt(i),p2db);
        % Write data to 3-D array:
        alldata_2db(1:maxrec,:,prof) = ...
            [ci' ti' p2db' profnumber*ones(size(p2db'))];
      end
 disp(['... done it. : ',num2str(fix(clock))])
 
 if 1 == 2
     % scu : I intend more processing of 2db data to go here
     % removing bad profiles, any drift corrections required
     % 
 end
 
 % Save 2db as binary
 disp(['Saving 2db file : ', outsave_2db])
 % file format profile time doy(rel to 1/1/06, yy mm dd hh mm ss
 % file alldata_2db c t p prof_num)
 
 save  (outsave_2db, 'alldata_2db', 'profile_time') 
 disp([''])
 disp(['Now run mmp_stage2_plot.m'])
 
