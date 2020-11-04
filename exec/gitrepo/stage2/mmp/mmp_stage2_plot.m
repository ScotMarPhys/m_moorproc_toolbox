% By default, MATLAB runs through seven colours in the plots, so 
% choosing multiples of seven for firstprof and lastprof makes 
% identifying odd profiles much easier. Also, the zoom tool in 
% MATLAB figures is expected to help greatly.
%
  format compact
  
  close all,  clear,    clear all,     
  moor     = 'eb2_4_200562';                
  operator = 'scu';
  disp(['Operator: ',operator])
 
% addpath(genpath ('C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data\exec\moor'))
% addpath(genpath ('C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data\exec\mfiles'))
% needs to be automated with startup
% Set prefix for Windows or UNIX:
%  prefix = 'C:\Documents and Settings\m221001\My Documents\Brown 2006\rapid\data'  %jym's PC
%  prefix = '/local1/rapid/data'         % UNIX
%   prefix = '/data32/rapid/data'; % D304
 prefix = '/local/users/pstar/Data/rpdmoc/rapid/data' % D334 test
  inpath = [prefix,'/moor/proc/',moor,'/mmp/'];
  outpath = [prefix,'/moor/proc/',moor,'/mmp/'];
  infofile = [prefix,'/moor/proc/',moor,'/',moor,'info.dat']; 
  infile = [inpath,'mmp_ctd_proc_2db.mat'];

  disp(['Loading : ',infile])
  load(infile); 

 % file format profile time doy(rel to 1/1/06, yy mm dd hh mm ss
 % file alldata_2db c t p prof_num)
 % extract 2db data to 2D arrays
 disp(['Extracting c t p to 2D array'])
 press = squeeze(alldata_2db(:,3,:));
 temp = squeeze(alldata_2db(:,2,:));
 cond = squeeze(alldata_2db(:,1,:));                                                             

  disp(['Calculating salin, potemp, sigma2'])
% salinity:
  c3515   = sw_c3515;
  cratio =cond/c3515;
  salin = sw_salt(cratio,temp,press);
    
% potential temperature:
  potemp  = sw_ptmp(salin,temp,press,0);

% potential density:
 pref = 2000;
 pden = sw_pden(salin,temp,press,pref);

  
% Hand-pick part of the profiles:  
  firstprof = 0*7 ;         % First profile. NB: Could be profile zero!
%  lastprof = 174           % Last profile eb2_3_200516
  lastprof = 156;            % eb2_4_200562
  nprof  = lastprof - firstprof +1;
                           % Colours cycle modulo 7. For identifying
                           % profiles, good to use multiples of 7 for
                           % firstprof and lastprof, so colours don't
                           % change
  disp(['Plotting profiles ',num2str(firstprof),' to ',num2str(lastprof)])
% Load plot arrays; introduce salinity offset to separate curves:

  temp_plot = temp(:,firstprof+1:lastprof+1);
  press_plot = press(:,firstprof+1:lastprof+1);
  salin_plot = salin(:,firstprof+1:lastprof+1);
  potemp_plot = potemp(:,firstprof+1:lastprof+1);
  pden_plot = pden(:,firstprof+1:lastprof+1);

  for prof = 1:nprof
      soffset = .1;
      pdenoffset = .1;
      salin_off(:,prof) = salin_plot(:,prof) + ...
          soffset*(prof-1);
      pden_off(:,prof) = pden_plot(:,prof) + ...
          pdenoffset*(prof-1);
  end
  
% Plot salinity profiles:
  figure (1), clf reset
  set(gcf, 'PaperOrientation', 'landscape', 'PaperPosition', [0 0 11 8.5]);  
  plot(salin_off,press_plot)
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
  plot(salin_plot,potemp_plot)
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
  plot(salin_off,temp_plot)
  title([moor,'; first profile: ',num2str(firstprof),', last profile: ',...
      num2str(lastprof), '; salinity vs. temperature, S offset by ',...
      num2str(soffset),' psu per profile'],'interpreter','none')
  xlabel('psu')
  ylabel('deg C')
  grid off                          % Important for visual inspection
  print ('-dpsc',[outpath,'salinity_offset_temperature.ps']) 
  
  
% extract first and last times from profile_time
yys=profile_time(2,1);mns=profile_time(3,1);dds=profile_time(4,1);
hhs=profile_time(5,1);mms=profile_time(6,1);sss=profile_time(7,1);
yye=profile_time(2,end);mne=profile_time(3,end);dde=profile_time(4,end);
hhe=profile_time(5,end);mme=profile_time(6,end);sse=profile_time(7,end);
jd_s  = julian(yys,mns,dds,hhs(1)+mms/60+sss/60/60);  % start time
jd_e  = julian(yye,mne,dde,hhe+mme/60+sse/60/60);     % end time
% convert time to jday
jd = julian(profile_time(2,:),profile_time(3,:),profile_time(4,:),...
    profile_time(5,:)+profile_time(6,:)/60+profile_time(7,:)/60/60);
jdg=repmat(jd,length(press),1); % correct dimensions for plotting
plot_interval = [2005 11 01 00;   % start time of time axis on plot
		 2006 06 01 00];  % end time of time axis on plot
jd1 = julian(plot_interval(1,:));
jd2 = julian(plot_interval(2,:)); 


figure(5)
z_grid = [2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 ...
    10.5 11 11.5 12 13 14 15 17.5 20 22.5 25 27.5];
subplot(211),[c,h]=contourf(jdg-jd1,press,potemp,z_grid);
set(gca,'Ydir','reverse')           % Pressure increaes downward
colorbarf(c,h);
xlim([0 jd2-jd1]); timeaxis(plot_interval(1,1:3)); 
potemp_anom=potemp-repmat(nanmean(potemp,2),1,157);
subplot(212),[c,h]=contourf(jdg-jd1,press,potemp_anom);
set(gca,'Ydir','reverse')
colorbarf(c,h);
xlim([0 jd2-jd1]); timeaxis(plot_interval(1,1:3));

figure(6)
z_grid = [34.5 34.6 34.7 34.8 34.9 35 35.1 35.2 35.3 35.4...
    35.5 35.75 36 36.5];
subplot(211),[c,h]=contourf(jdg-jd1,press,salin,z_grid);
set(gca,'Ydir','reverse')           % Pressure increaes downward
colorbarf(c,h);
xlim([0 jd2-jd1]); timeaxis(plot_interval(1,1:3)); 
%mean_salin=repmat(nanmean(salin,2),1,157);
salin_anom=salin-repmat(nanmean(salin,2),1,157);
subplot(212),[c,h]=contourf(jdg-jd1,press,salin_anom);
set(gca,'Ydir','reverse')
colorbarf(c,h);
xlim([0 jd2-jd1]); timeaxis(plot_interval(1,1:3));

