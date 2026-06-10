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
% mmp_stage2 creates a number of matlab figures and prints them into
% colour ps files.  
%
% By default, MATLAB runs through seven colours in the plots, so 
% choosing multiples of seven for firstprof and lastprof makes 
% identifying odd profiles much easier. Also, the zoom tool in 
% MATLAB figures is expected to help greatly.
%
% Work in progress.  Signifant portions of the file have been 
% deactivated by setting loops of the sort if 1==2; body; end,
% which formally are valid code but obviously never execute "body".
% Some of that stuff, such as pressure gridding and plotting, might
% one day still be useful.  Just dunno.
% 
% 16 - 17 & 26 March 2006 
% Jochem Marotzke, MPI-M (jym)

   
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
%   prefix = '/data32/rapid/data' % D304
  prefix = '/local/users/pstar/Data/rpdmoc/rapid/data' % D334 test
  inpath = [prefix,'/moor/proc/',moor,'/mmp/']
 % inpath = [prefix,'\moor\proc\',moor,'\mmp\']
  outpath = [prefix,'/moor/proc/',moor,'/mmp/']
  infofile = [prefix,'/moor/proc/',moor,'/',moor,'info.dat'];  
  infile = [inpath,'mmp_ctd.txt']
  outfile = [outpath,'mmp_ctd_proc.txt']
  outsave = [outpath,'mmp_ctd_proc.mat']
  outdiary = [outpath,'mmp_stage2_log.txt']
  diary (outdiary)

  alldata = load(infile); 
  
%   yyyy = alldata (:,1);
%   mm   = alldata (:,2);
%   dd   = alldata (:,3);
%   hh   = alldata (:,4); 
%   mm   = alldata (:,5);  
%   ssss = alldata (:,6);     
%   cond = alldata (:,7);     
%   temp = alldata (:,8);    
%   press= alldata (:,9);   
%   prof = alldata (:,10);
% 
% --- Calculate salinity:
  c3515   = sw_c3515;
  c_ratio = alldata(:,7)/c3515;
  alldata(:,11)  = sw_salt(c_ratio,alldata(:,8),alldata(:,9));
  
% --- Calculate potential temperature:
  alldata(:,12)  = sw_ptmp(alldata(:,7),alldata(:,8),alldata(9),0);

% --- Calculate potential density:
  pref = 2000;
  alldata(:,13)  = sw_pden(alldata(:,7),alldata(:,8),alldata(9),pref);


% Determine first and profiles:  
  firstprof = alldata(1,10)
  lastprof  = alldata(end,10)
  nprof  = lastprof - firstprof +1;

  
% Find beginning and end indices for profiles:  
  profbeg = [1, find(diff(alldata(:,10)'))+1];
  profend = [find(diff(alldata(:,10)')), length(alldata(:,10)')];  
  lrec   = (profend - profbeg) + 1;           % Length of each profile
  maxrec = max(lrec); 
  
% Rearrange data (back) into profiles:
% First, initialise overall 3-D array:
  profdata = NaN(maxrec,13,nprof);
  
% Then, loop through profiles and fill 3-D array.
% Matlab does not like index zero, so the operations are done with
% profile count up by one.  But output is strictly according to original
% profile number.

  for profnumber = firstprof : lastprof         % Loop through profiles
      prof = profnumber+1;      
      profdata (1:lrec(prof),:,prof) = ...
          alldata(profbeg(prof):profend(prof),:);
  end
  
  yyyy = profdata (:,1,:);
  mm   = profdata (:,2,:);
  dd   = profdata (:,3,:);
  hh   = profdata (:,4,:); 
  mm   = profdata (:,5,:);  
  ssss = profdata (:,6,:);     
  cond = profdata (:,7,:);     
  temp = profdata (:,8,:);    
  press= profdata (:,9,:);   
  prof = profdata (:,10,:);
  sal  = profdata (:,11,:);
  ptem = profdata (:,12,:);
  pden = profdata (:,13,:);
  
  clear profdata alldata

  yyyy = reshape (yyyy, maxrec, nprof);
  mm   = reshape (mm  , maxrec, nprof);
  dd   = reshape (dd  , maxrec, nprof);
  hh   = reshape (hh  , maxrec, nprof); 
  mm   = reshape (mm  , maxrec, nprof);  
  ssss = reshape (ssss, maxrec, nprof);     
  cond = reshape (cond, maxrec, nprof);     
  temp = reshape (temp, maxrec, nprof);    
  press= reshape (press, maxrec, nprof);   
  prof = reshape (prof, maxrec, nprof);
  sal  = reshape (sal , maxrec, nprof);

  
                                                            if 1==0s
% jym: What follows is stage 2 type processing. Must see through it.
% But that's for another day.

% sort profiles / remove multiple pressures 
  
  [p,pII] = sort(p);
  move   = find(diff(p)~=0);

  eliminate.end_values(prof) = length(p) - length(move);
  
  c = c(pII(move));
  t = t(pII(move));
  p = p(move);

% basic despiking 
  
  c_lim = [20 70];
  t_lim = [-2 40];
  p_lim = [0 6000];
  s_lim = [30 40];


  valp = find(p>p_lim(1) & p<p_lim(2)); 
  
  pelim = length(p) - length(valp);
  eliminate.p(prof) = pelim;
  
  p    = p(valp);
  t    = t(valp);
  c    = c(valp);
  
  valt = find(t>t_lim(1) & t<t_lim(2));
  valc = find(c>c_lim(1) & c<c_lim(2)); 

  telim = length(t) - length(valt);
  eliminate.t(prof) = telim;
  celim = length(c) - length(valc);
  eliminate.t(prof) = celim;
  

% interpolate and filter
 
  pgrid = [ 120:1:2500]; % pressure grid
  pstep = 5;       % pressure sampling frequency [1/dbar]
  pcut  = 1/2;     % pressure cut frequency [1/dbar]

  pii               = [min(p):1/pstep:max(p)]';
  lpi              = length(pii);
  pint(1:lpi,prof)   = pii;
  ci(1:lpi,prof)   = interp1(p(valc),c(valc),pii);
  ti(1:lpi,prof)   = interp1(p(valt),t(valt),pii);
  
  if lpi >  pstep / pcut *2 
    C(:,prof) = interp1(pii,auto_filt(ci(1:lpi,prof),pstep,pcut),pgrid');
    T(:,prof) = interp1(pii,auto_filt(ti(1:lpi,prof),pstep,pcut),pgrid');
  else 
     T(:,prof) = NaN;
     C(:,prof) = NaN;
  end  
     
                                                            end

  
% Hand-pick part of the profiles:  
  firstprof = 0*7          % First profile. NB: Could be profile zero!
%  lastprof = 174           % Last profile eb2_3_200516
  lastprof = 156            % eb2_4_200562
  nprof  = lastprof - firstprof +1;
                           % Colours cycle modulo 7. For identifying
                           % profiles, good to use multiples of 7 for
                           % firstprof and lastprof, so colours don't
                           % change
  
% Load plot arrays; introduce salinity offset to separate curves:

  temp_plot = temp(:,firstprof+1:lastprof+1);
  press_plot = press(:,firstprof+1:lastprof+1);
  sal_plot = sal(:,firstprof+1:lastprof+1);
  ptem_plot = ptem(:,firstprof+1:lastprof+1);
  pden_plot = pden(:,firstprof+1:lastprof+1);

  for prof = 1:nprof
      soffset = .1;
      pdenoffset = .1;
      sal_off(:,prof) = sal_plot(:,prof) + ...
          soffset*(prof-1);
      pden_off(:,prof) = pden_plot(:,prof) + ...
          pdenoffset*(prof-1);
  end
  
% Plot salinity profiles:
  figure (1), clf reset
  set(gcf, 'PaperOrientation', 'landscape', 'PaperPosition', [0 0 11 8.5]);  
  plot(sal_off,press_plot)
  set(gca,'Ydir','reverse'),           % Pressure increaes downward 
  title([moor,'; first profile: ',num2str(firstprof),', last profile: ',...
      num2str(lastprof), '; salinity vs. pressure, S offset by ',...
      num2str(soffset),' psu per profile'],'interpreter','none')
  xlabel('psu')
  ylabel('dbar')
  grid
  print ('-dpsc',[outpath,'salinity_profiles.ps']) 
  
% Plot potential density profiles:
  figure (2), clf reset
  set(gcf, 'PaperOrientation', 'landscape', 'PaperPosition', [0 0 11 8.5]);  
  plot(pden_off,press_plot)
  set(gca,'Ydir','reverse'),           % Pressure increaes downward 
  title([moor,'; first profile: ',num2str(firstprof),', last profile: ',...
      num2str(lastprof), '; sigma (',num2str(pref),...
      ' db) vs. pressure, pden offset by ',...
      num2str(pdenoffset),' kg/m^3 per profile'],'interpreter','none')
  xlabel('kg/m^3')
  ylabel('dbar')
  grid
  print ('-dpsc',[outpath,'potential_density_profiles.ps']) 
  
  
% Plot salinity vs. temperature:
  figure (3), clf reset
  plot(sal_plot,ptem_plot)
  title([moor,'; first profile: ',num2str(firstprof),', last profile: ',...
    num2str(lastprof), '; salinity vs. pot. temperature'],'interpreter','none')
  xlabel('psu')
  ylabel('deg C')
  hold on
% Plot isopycnals
  aa = axis;
  s1 = aa(1); s2 = aa (2); t1 = aa(3); t2 = aa(4);
  lvec = 100;
  svec = [s1:(s2-s1)/lvec:s2];
  tvec = [t1:(t2-t1)/lvec:t2]';
  s_contour = ones(lvec+1,1)*svec;
  t_contour = tvec*ones(1,lvec+1);
  pden_contour = sw_pden(s_contour,t_contour,2000,2000);
  contour(svec,tvec,pden_contour,[1030:.5:1040],'k-')
  contour(svec,tvec,pden_contour,[1034:.05:1038],'k:')
  grid off
  print ('-dpsc',[outpath,'salinity_pot_temperature.ps']) 
  
% Plot offset salinity vs. temperature:
  figure (4), clf reset
  set(gcf, 'PaperOrientation', 'landscape', 'PaperPosition', [0 0 11 8.5]);  
  plot(sal_off,temp_plot)
  title([moor,'; first profile: ',num2str(firstprof),', last profile: ',...
      num2str(lastprof), '; salinity vs. temperature, S offset by ',...
      num2str(soffset),' psu per profile'],'interpreter','none')
  xlabel('psu')
  ylabel('deg C')
  grid off                          % Important for visual inspection
  print ('-dpsc',[outpath,'salinity_offset_temperature.ps']) 
  
                                                if 1==0

% --- calculate salinity

[m,n]   = size(T);
c3515   = sw_c3515;
c_ratio = C/c3515;
S      = sw_salt(c_ratio,T,pgrid'*ones(1,n));
 

 tmean = nanmean(T,2);
 smean = nanmean(S,2);

TA = T - tmean * ones(1,n);
SA = S - smean * ones(1,n);

% save data (legacy from TK)


% save (outsave, 'ci', 'ti', 'pint', 'pgrid', 'TA', 'SA', 'S',...
%    'T', 'C', 'tmean', 'smean',  'jd_on', 'jd_off', 'eliminate')

fprintf(1,'\n')


% --------- graphics --------------------

 jd1 = jd_on(1:2:end);
 jd2 = jd_off(2:2:end);
 jd_deep = sort([jd1' jd2']);

                                                 if 0==1
 
 figure(2); clf reset
 subplot(2,1,2)
 hold on
 contourf(jd_on-jd_on(1),pgrid(1:10:end),SA(1:10:end,1:end),[-.1:0.01:.1])
% contourf(jd_on-jd_on(1),pgrid(1:10:end),S(1:10:end,1:end)-36,[-1.1:0.05:1.1])
 set(gca,'ydir','reverse') 
 ylabel('Pressure [dbar]')
 xlabel('Time [days]')
 title('Salinity anomaly')
 shading flat
 xlim([0 190])
 ylim([120 2500])
 orient landscape
 colorbar

 subplot(2,1,1)
 hold on
 contourf(jd_on-jd_on(1),pgrid(1:10:end),TA(1:10:end,1:end),[-.5:0.05:.5])
% contourf(jd_on-jd_on(1),pgrid(1:10:end),T(1:10:end,1:end))
 set(gca,'ydir','reverse') 
 ylabel('Pressure [dbar]')
 title('Temperature anomaly')
 xlabel('Time [days]')
 shading flat
 xlim([0 190])
 ylim([120 2500])
 orient landscape
 colorbar
                                                 end
 % print ('-depsc', outprint)
                                                end

 diary off
% generate path for mmp processing
% needs to be automated with startup
% jym 12 March 2006
% addpath(genpath ('C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data\exec\moor'))
% addpath(genpath ('C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data\exec\mfiles'))