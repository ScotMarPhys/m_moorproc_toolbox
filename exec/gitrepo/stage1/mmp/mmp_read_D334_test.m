% mmp_read.m
%
% Loads mmp data:  Conductivity Temperature Pressure and converts to 
% rapid data format.  
%
% The user must specify an input directory where the Cnnnnnnn.txt
% ctd files are stored, and an output directory, by way of the 
% mooring identifier.  The syntax for the file structure is such
% that, if properly set up on each machine, only a prefix must be 
% changed to move from UNIX to Windows and back.  Forward slash
% was no problem on Windows XP.
%
% Data are written into a text file, mmp_ctd.txt, with the
% first 6 columns being timestamp (YYYY MM DD HH MM SS), then C, T, P, 
% and profile number.  The resulting file is very large, e.g., 
% O(100MB) for 174 profiles. All data written to mmp_ctd.txt are held
% in a 3-D array alldata(maxrec,10,lastprof+1);     
% Dimension 1: max maxrec records per profile; 
% Dimension 2: 6 variables for time; c,t,p; profile number
% Dimension 3: Number of profiles
%
% MATLAB allows dynamic memory allocation, so one can go beyond 
% the size set initially.  But that makes the program extremely slow.  
% Therefore, one should initialise so large that all is covered.
%
% A binary .mat file is also written with alldata and the Julian 
% Day stamps for beginning and end of profiles (jd_on and jd_off).
%
% Some very basic extraction information is stored in a log file.
%
% The timestamps for each data record are calculated from the recorded 
% "turn on" and "turn off" times of the Cnnnnnnn.txt files, by 
% dividing this time span by the number of records. The basic order
% parameter in the output file mmp_ctd.txt is record number or,
% equivalently, time progression.  No other processing has been
% done at this stage.
% 
% mmp_read is set up to work with arbitrary  start and end 
% profiles (firstprof and lastprof, respectively), to be set in
% the m-file.  Operator could not be bothered setting this up
% interactively.  Notice that the first mmp profile, the initial 
% downward profile, has number zero.  Since matlab doesn't like 
% index zero, this had to be accounted for in the loops through 
% profiles. But output is strictly according to original
% profile number.
%
% mmp_read creates a number of matlab figures and prints them into
% colour ps files.  The most useful for initial quality control are 
% expected to be the min/max pressures of all profiles, the runtimes 
% for the profiles, and the conductivity-temperature plots, both
% straight (on top of each other) and offset by 1 mmho/cm per profile.  
% The latter ("grid off" is crucial!) makes possible a purely visual
% detection of odd conductivities through irregularities in an 
% otherwise very regular pattern.
%
% The final stage of the processing must be done with all profiles,
% but the initial QC will suggest running a subset of the profiles,
% making it possible to identify the profiles where odd things happen.
% By default, MATLAB runs through seven colours in the plots, so 
% choosing multiples of seven for firstprof and lastprof makes 
% identifying odd profiles much easier. Also, the zoom tool in 
% MATLAB figures is expected to help greatly, especially in the 
% figure showing the pressure time series.
% 
% 11 - 16 March 2006 
% Jochem Marotzke, MPI-M (jym)

  close all,   clear,    clear all,   format compact
  moor     = 'eb2_4_200562';                
  operator = 'scu';
  disp(['Operator: ',operator])
  disp(['Start processing [YYYY MM DD HH MM SS]: ',num2str(fix(clock))])
  
% addpath(genpath ('C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data\exec\moor'))
% addpath(genpath ('C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data\exec\mfiles'))
% needs to be automated with startup

% Set prefix for Windows or UNIX:

%  prefix = 'C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data'  %jym's PC
  %prefix = '/local1/rapid/data'         % UNIX
  prefix = '/local/users/pstar/Data/rpdmoc/rapid/data' % D334 test
  inpath = [prefix,'/moor/raw/d304/mmp/'];
  outpath = [prefix,'/moor/proc/',moor,'/mmp/']
  infofile = [prefix,'/moor/proc/',moor,'/',moor,'info.dat'];  
  outfile = [outpath,'mmp_ctd.txt']
  outsave = [outpath,'mmp_ctd_matlab.mat']
  outdiary = [outpath,'mmp_stage1_log.txt']
  diary (outdiary)

  firstprof = 0*7 % First profile. NB: Could be profile zero!
  %lastprof = 156           % eb2_4_200562
  lastprof = 7           % eb2_4_200562
  %lastprof = 174          % eb2_3_200516 
                           % Colours cycle modulo 7. For identifying
                           % profiles, good to use multiples of 7 for
                           % firstprof and lastprof, so colours don't
                           % change
  
% Initialise data array.  NB: MATLAB allows dynamic memory allocation,
% so one can go beyond the size set here.  But that makes the program
% extremely slow.  Therefore, one should initialise so large that
% all is covered.

  lprof = lastprof - firstprof +1;      
  maxrec = 14000;          
  maxvar = 10;
  alldata = NaN(maxrec,10,lastprof+1);     
% Dimension 1: max maxrec records per profile; 
% Dimension 2: 6 variables for time; c,t,p; profile number
% Dimension 3: Number of profiles
  
  fid_ctd = fopen(outfile,'w');         % Open output file
  fprintf(fid_ctd,'%% YYYY MM DD HH MM  SSSS     COND     TEMP    PRESS   PROF \n');

  
  for profnumber = firstprof : lastprof         % Loop through profiles

% Matlab does not like index zero, so the operations are done with
% profile count up by one.  But output is strictly according to original
% profile number.
  prof = profnumber+1;              

  fprintf('\r Profile: %2.2d',profnumber)
  fprintf('\r')
   
% Get data (this is code originating from Torsten Kanzow)
% First, open input file and read entire contents into string "linenumber":

  fid = fopen([inpath,'/C',sprintf('%7.7d',profnumber),'.TXT'],'r');
  linenumber = fscanf(fid,'%c');
  fclose(fid);

% Next, look for linefeed (end of line) and decimal points in linenumber 
% Once their positions are known, they can be read into the data arrays. 

  ret = sprintf('\n');

  xx     = findstr(ret,linenumber);
  point  = findstr(linenumber,'.');

  beginning = find(xx<point(1));
  finish   = find(xx>point(end));
  data   = str2num(linenumber(xx(beginning(end)):xx(finish(1))));

  c      = data(:,1);  % cond
  t      = data(:,2);  % temp
  p      = data(:,3);  % pres 

  disp(['Active records: ',num2str(length(c))])
  nrec(prof)  = length(c);        % Number of active records
  
  if nrec(prof) > maxrec, 
      disp('Too many records, stop with ctrl+C and increase ')
      disp(    'first dimension of alldata'), 
      pause, 
  end
  
% Get cast start and end times:
  
  onI(1) = findstr(linenumber,'on at') + 5;
  III    = find(xx>onI(1));
  onI(2) = xx(III(1))-1;
  onstr  = linenumber(onI(1):onI(2));
  disp(['Time on: ',onstr])
  ontime(prof,:)  = datevec(datenum(onstr));
  
  offI(1) = findstr(linenumber,'ff at') + 5;
  III    = find(xx>offI(1));
  offI(2) = xx(III(1))-1;
  offstr  = linenumber(offI(1):offI(2));
  disp(['Time off: ',offstr])
  offtime(prof,:)  = datevec(datenum(offstr));

% Create timestamp for every data record, using utility datenum,
% which counts in days:  
  time_increment = (datenum(offstr) - datenum(onstr))/(nrec(prof)-1);
  disp(['Time increment: ',num2str(time_increment*86400),' s'])
  
  timestamp = [datenum(onstr):time_increment:datenum(offstr)]';
  timevec = datevec(timestamp);
   
% Write data to ascii outputfile:
  for nr=1:nrec(prof)
      fprintf(fid_ctd,'  %4.0f %2.2g %2.2g %2.2g %2.2g %6.2f % 8.4f % 8.4f % 8.3f %3.0f\n' ,...
          [timevec(nr,:) c(nr) t(nr) p(nr) profnumber]);
  end

% Write data to 3-D array:
  alldata(1:nrec(prof),:,prof) = ...
                [timevec c t p profnumber*ones(nrec(prof),1)];

            end                   % Loop through profiles

            
  jd_on = julian([ontime(firstprof+1:lastprof+1,1:3) ... 
      hms2h(ontime(firstprof+1:lastprof+1,4:6))]);
  jd_off = julian([offtime(firstprof+1:lastprof+1,1:3) ...
      hms2h(offtime(firstprof+1:lastprof+1,4:6))]);

  fclose(fid_ctd);              % Close output file

% Save binary (much easier than ascii!):
  save  (outsave, 'alldata', 'jd_on', 'jd_off')    

  
% --------- graphics --------------------

% First, plot profile travel times, using Julian Day calculation:            

  figure(1), clf reset  
  
% Remember: Calculations are done with profile count increased
% by one, to avoid index zero.  But labels count the profiles correctly. 
  if mod(firstprof,2)==0
      downprof = [firstprof:2:lastprof];
      upprof = [firstprof+1:2:lastprof];
      plot(downprof,(jd_off(1:2:lprof) - jd_on(1:2:lprof))*24,'ko')
      hold on
      plot(downprof,(jd_off(1:2:lprof) - jd_on(1:2:lprof))*24,'k')
      plot(upprof,(jd_off(2:2:lprof) - jd_on(2:2:lprof))*24,'rx')
      plot(upprof,(jd_off(2:2:lprof) - jd_on(2:2:lprof))*24,'r')
      hold off  
  else
      upprof = [firstprof:2:lastprof];
      downprof = [firstprof+1:2:lastprof];
      plot(downprof,(jd_off(2:2:lprof) - jd_on(2:2:lprof))*24,'ko')
      hold on
      plot(downprof,(jd_off(2:2:lprof) - jd_on(2:2:lprof))*24,'k')
      plot(upprof,(jd_off(1:2:lprof) - jd_on(1:2:lprof))*24,'rx')
      plot(upprof,(jd_off(1:2:lprof) - jd_on(1:2:lprof))*24,'r')
      hold off   
  end
  title([moor,', Time taken for profiles; black circles: ',...
      'down; red crosses: up'],'interpreter','none')
  xlabel('Profile Number')
  ylabel('hours')
  grid
  print ('-dpsc',[outpath,'profile_times.ps']) 
  
 
 % Now, use the data stored in alldata. Easiest to transfer them
 % back into 2-dimensional fields. Remember that 
 % data in second dimension in alldata are held in this order:
 % YYYY MM DD HH MM  SSSS     COND     TEMP    PRESS   PROF 
 %
 
  cond = NaN(maxrec,lastprof+1);        % Initialise conductivity matrix
  temp = NaN(maxrec,lastprof+1);        % Initialise temperature matrix
  press = NaN(maxrec,lastprof+1);       % Initialise pressure matrix
  temp_plot = NaN(maxrec,lastprof+1);   % Initialise temp offset matrix
  cond_plot = NaN(maxrec,lastprof+1);   % Initialise cond offset matrix
  runtime = NaN(maxrec,lastprof+1);     % Initialise time matrix
  
  for profnumber = firstprof:lastprof
      prof = profnumber + 1;
      runtime(:,prof) = datenum(alldata(:,1:6,prof));
      runtime(:,prof) = runtime(:,prof) - ...
          datenum('01-Jan-2005', 'dd-mmm-yyyy');
      cond(:,prof) = alldata(:,7,prof);
      temp(:,prof) = alldata(:,8,prof);
      press(:,prof) = alldata(:,9,prof);
% Introduce temperature and conductivity offsets to separate curves:
      toffset = 1;
      temp_plot(:,prof) = temp(:,prof) + ...
          toffset*(profnumber-firstprof);
      cond_plot(:,prof) = cond(:,prof) + ...
          toffset*(profnumber-firstprof);
  end
% The transfer to 2-dim fields may be done using the squeeze function
% scu 30/5/06 but not yet implemented
% cond = squeeze(alldata(:,7,:));
% temp = squeeze(alldata(:,8,:));
% press = squeeze(alldata(:,9,:));

% Plot maximum and minimum pressures of each profile:
  figure (2), clf reset
  maxpress = max(press);
  minpress = min(press);
  plot([firstprof:lastprof],maxpress(firstprof+1:lastprof+1),'ko',...
      [firstprof:lastprof],minpress(firstprof+1:lastprof+1),'rx')
  hold on
  plot([firstprof:lastprof],maxpress(firstprof+1:lastprof+1),'k',...
      [firstprof:lastprof],minpress(firstprof+1:lastprof+1),'r')
  set(gca,'Ydir','reverse'),           % Pressure increaes downward 
  title([moor,', Black circles: maximum depth; red crosses: ',...
      'minimum depth'],'interpreter','none')
  xlabel('Profile Number')
  ylabel('dbar')
  grid
  print ('-dpsc',[outpath,'max_min_press.ps']) 
  
% Plot pressure time series:
  figure (3), clf reset
  set(gcf, 'PaperOrientation', 'landscape', 'PaperPosition', [0 0 11 8.5]);  
  plot(runtime,press)
  set(gca,'Ydir','reverse'),           % Pressure increaes downward 
  title([moor,'; first profile: ',num2str(firstprof),', last profile: ',...
      num2str(lastprof),'; Pressure time series vs. time'],...
      'interpreter','none')
  xlabel('Yearday 2005')
  ylabel('dbar')
  grid
  print ('-dpsc',[outpath,'pressure_time_series.ps']) 
  
% Plot temperature profiles:
  figure (4), clf reset
  set(gcf, 'PaperOrientation', 'landscape', 'PaperPosition', [0 0 11 8.5]);  
  plot(temp_plot,press)
  set(gca,'Ydir','reverse'),           % Pressure increaes downward 
  title([moor,'; first profile: ',num2str(firstprof),', last profile: ',...
      num2str(lastprof),'; Temperature vs. pressure, T offset by ',...
      '1deg per profile'],'interpreter','none')
  xlabel('dec C')
  ylabel('dbar')
  grid
  print ('-dpsc',[outpath,'temperature_profiles.ps']) 
  
% Plot conductivity profiles:
  figure (5), clf reset
  set(gcf, 'PaperOrientation', 'landscape', 'PaperPosition', [0 0 11 8.5]);  
  plot(cond_plot,press)
  set(gca,'Ydir','reverse'),           % Pressure increaes downward 
  title([moor,'; first profile: ',num2str(firstprof),', last profile: ',...
      num2str(lastprof),'; Conductivity vs. pressure, '...
      'C offset by 1 mmho/cm per profile'],'interpreter','none')
  xlabel('mmho/cm')
  ylabel('dbar')
  grid
  print ('-dpsc',[outpath,'conductivity_profiles.ps']) 
  
% Plot conductivity vs. temperature:
  figure (6), clf reset
  plot(cond,temp)
  title([moor,'; first profile: ',num2str(firstprof),', last profile: ',...
    num2str(lastprof), '; Conductivity vs. temperature'],'interpreter','none')
  xlabel('mmho/cm')
  ylabel('deg C')
  grid
  print ('-dpsc',[outpath,'conductivity_temperature.ps']) 

% Plot offset conductivity vs. temperature:
  figure (7), clf reset
  plot(cond_plot,temp)
  title([moor,'; first profile: ',num2str(firstprof),', last profile: ',...
      num2str(lastprof), '; Conductivity vs. temperature, '...
      'C offset by 1 mmho/cm per profile'],'interpreter','none')
  xlabel('mmho/cm')
  ylabel('deg C')
  grid off                          % Important for visual inspection
  print ('-dpsc',[outpath,'conductivity_offset_temperature.ps']) 

  diary off