% process_caldip.m
% this script reads a single SBE-9 data file formatted in the standard
% SeaBird '.cnv' format
% then looks for a text file containing a list of the microcat files to be
% checked against the SBE-9
% the files do not have to be in any particular location, but the list file
% should be in the same directory as the microcat files
% the script will generate four figures for each microcat, so if the
% computer is getting bogged down with too many open figure windows,
% comment out the plotting sections
% only figure number 2 will be written to a postscript format
% 
% Adam Houk
% 4/17/2011

c_outliers = 0.05 ; % largest ctd-mc cond difference allowed in mS/cm
nsd = 3; % standard deviation limit for second edit

pdiff = 0.4; % CTD pressure sensor is deeper in water relative to microcats

% set conductivity at s=35, t=15 for salt calculation
c3515 = sw_c3515();

diary on
% load 1-sec CTD cast file
[sbe9_filename,sbe9_pathname] = uigetfile('*.cnv','Load 1-sec Averaged CTD cast');
[sb,cfg] = load_sbe9([sbe9_pathname, sbe9_filename]);
	disp('****************************************')
	disp(['Loaded SBE-9 CTD cast: ' sbe9_filename]);
	disp('****************************************')
% filter ctd data to same period as microcats (1 sample every 10s)
% NB this bit is necessary before running interp1 as sometimes there is a
% duplicate time in the ctd file which can cause interp1 to fail
%disp(['Filtering CTD data using filtfilt WindowSize = 10'])
windowSize = 10;
filt.p_ctd = filtfilt(ones(1,windowSize)/windowSize,1,sb.prDM);
filt.t_ctd = filtfilt(ones(1,windowSize)/windowSize,1,sb.t190C); %temporarily using secondary temperature probe
filt.c_ctd = filtfilt(ones(1,windowSize)/windowSize,1,sb.c0S.*10);
filt.jd_ctd = filtfilt(ones(1,windowSize)/windowSize,1,sb.timeJ);

% set 10-sec time interval for interpolation
TI = [filt.jd_ctd(1):(10/86400):filt.jd_ctd(end)];

% interpolate CTD data to new time base
p_ctdi = interp1(filt.jd_ctd,filt.p_ctd,TI);
t_ctdi = interp1(filt.jd_ctd,filt.t_ctd,TI);
c_ctdi = interp1(filt.jd_ctd,filt.c_ctd,TI);

% load microcats using list file containing SBE37 filenames plus one header
% line
[sbe37_listfile,sbe37_pathname] = uigetfile('*.*','Load SBE-37 list file');
cd(sbe37_pathname)
fidl = fopen(sbe37_listfile,'r');
i=1;

% read header line containing caldip cast number and time adjustment (in seconds)
% for SBE37s relative to CTD profiler (which we assume has the correct time)
% for example: if SBE37 clock is fast by 60 sec, set value to -60
hdr = textscan(fgetl(fidl),'%s %f')
caldip_id = hdr{1}{1};
t_offset = hdr{2};

str = fgetl(fidl);
while ischar(str)
	disp('****************************************')
	disp(['Loading MicroCat: ' str]);
	disp('****************************************')
	[mc{i},mc_cfg{i}] = load_microcat(str);
	mc_datevec{i} = datevec(strcat(mc{i}{4},{' '},mc{i}{5}),'dd mmm yyyy HH:MM:SS');
	dec_hr{i} = mc_datevec{i}(:,4) + (mc_datevec{i}(:,5)/60) + (mc_datevec{i}(:,6)./3600);
	mc_jd{i} = julian(mc_datevec{i}(:,1),mc_datevec{i}(:,2),mc_datevec{i}(:,3),dec_hr{i}) ...
	- julian(mc_datevec{i}(1,1),1,0,0);
	mc_sd{i} = datenum(mc_datevec{i});
	mc_sal{i} = sw_salt(((mc{i}{2}.*10)./sw_c3515),mc{i}{1},mc{i}{3});
	
	mc{i}{2} = mc{i}{2} .* 10; % convert microcat conductivity to mS/cm
	
	% apply time correction
	if t_offset > 0 || t_offset < 0
		mc_jd{i} = mc_jd{i} + (t_offset/86400);
		mc_sd{i} = mc_sd{i} + (t_offset/86400);
	end
	
	% put data in matrix form
	eval(['sn' mc_cfg{i}.temp_sn ' = cell2mat([mc_sd(i) mc{i}(1:3) mc_sal(i)]);'])
	
	% interpolate microcat data to new time base
	p_mci{i} = interp1(mc_jd{i},mc{i}{3},TI);
	t_mci{i} = interp1(mc_jd{i},mc{i}{1},TI);
	c_mci{i} = interp1(mc_jd{i},mc{i}{2},TI);
	
	i=i+1;
	str = fgetl(fidl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate rate of change of pressure of CTD
n_ctd = length(p_ctdi);
dpdt(1)=NaN;
dpdt(n_ctd)=NaN;
for i=2:n_ctd-1
	dpdt(i)=(p_ctdi(i+1)-p_ctdi(i-1))/(86400*(TI(i+1)-TI(i-1)));
end

n_mc = length(mc); % number of microcats used
for i = 1:n_mc
	% create mask, keep only data where dp/dt~0 & remove clearly rubbish C
	bad_data{i}=ones(size(p_ctdi));
	kk = find(abs(dpdt) > 0.02);
	bad_data{i}(kk) = NaN;
	kk = find(c_ctdi < 30 & c_ctdi > 65);
	bad_data{i}(kk) = NaN;
	kk = find(c_mci{i} < 30 & c_mci{i} > 65);
	bad_data{i}(kk) = NaN;
	bad_data{i}(1:floor(length(bad_data{i})/3)) = NaN;

	% calculate ctd cond (mS/cm) - microcat cond (S/m)
	% K=c_ctd/cond and cat cond_corr = K*cond
	% take mean of K over all stations apply to data and replot to see effect
    i
    whos bad_data c_ctdi c_mci
	c_diff{i} = (c_ctdi - c_mci{i}) .* bad_data{i};
	cd_mn{i} = nanmean(c_diff{i});
	cd_sd{i} = nanstd(c_diff{i});
	K{i} = (c_ctdi ./ c_mci{i}) .* bad_data{i};
	K_mn{i} = nanmean(K{i});
	K_sd{i} = nanstd(K{i});

	% reject extravagent outliers
	kk = find(c_diff{i} > c_outliers);
	bad_data{i}(kk) = NaN;

	% recompute mean and standard deviations of c_diff and K
	c_diff{i} = (c_ctdi - c_mci{i}) .* bad_data{i};
	cd_mn{i} = nanmean(c_diff{i});
	cd_sd{i} = nanstd(c_diff{i});
	K{i} = (c_ctdi ./ c_mci{i}) .* bad_data{i};
	K_mn{i} = nanmean(K{i});
	K_sd{i} = nanstd(K{i});

	%  then remove data outside +/- nsd standard deviations of new mean

	%kk = find(abs(c_diff) > (cd_mn+(nsd*cd_sd)));  this is a logic bug..
	kk = find(c_diff{i} > cd_mn{i} + nsd*cd_sd{i} | c_diff{i} < cd_mn{i} - nsd*cd_sd{i}); %BJ
	bad_data{i}(kk) = NaN;

	% recompute mean and standard deviations of c_diff and K
	c_diff{i} = (c_ctdi - c_mci{i}) .* bad_data{i};
	cd_mn{i} = nanmean(c_diff{i});
	cd_sd{i} = nanstd(c_diff{i});
	K{i} = (c_ctdi ./ c_mci{i}) .* bad_data{i};
	K_mn{i} = nanmean(K{i});
	K_sd{i} = nanstd(K{i});

	% compute ctd-mc temperature and pressure differences
	t_diff{i} = (t_ctdi - t_mci{i}) .* bad_data{i};
	p_diff{i} = (p_ctdi - p_mci{i}) .* bad_data{i};
end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	
	disp('*************************************************************')
	disp(['Now carefully check time base alignment of CTD and microcat'])
	disp(['If they dont agree count the number of datacycles and enter this'])
	disp(['as mcat_tcorr and rerun programme. Use Figure(1).'])

for i = 1:n_mc
	% first get a data set without NaN for polyfit 
	ii = find(~isnan(bad_data{i}));
	p_mc_pf{i} = p_mci{i}(ii);
	t_mc_pf{i} = t_mci{i}(ii);
	p_ctd_pf{i} = p_ctdi(ii);
	t_ctd_pf{i} = t_ctdi(ii);
	p_diff_pf{i} = p_diff{i}(ii);
	t_diff_pf{i} = t_diff{i}(ii);

	n = length(p_diff_pf);
	disp(['Number of data points in calibration : ',num2str(n)])

	% degree 2 - quadratic fit for pressure hysteresis
	[pcoeff{i},sp{i}] = polyfit(p_mc_pf{i},p_diff_pf{i},2);
	p_fit{i} = polyval(pcoeff{i},p_mc_pf{i});

	disp('*************************************************************')
	disp([' Conductivity slope correction : ',num2str(K_mn{i})])
	disp(['Pressure fit : dpress(ctd-mc) = A + B*press_mc + C*press_mc^2'])

	C = pcoeff{i}(1); B = pcoeff{i}(2); A = pcoeff{i}(3);
	disp (['A = ',num2str(A),' B = ',num2str(B),' C  ',num2str(C)])

	% degree 2 - quadratic fit for temperature differences against ctd
	% temperature.
	[ptcoeff{i},st{i}] = polyfit(t_mc_pf{i},t_diff_pf{i},2);

	t_fit{i} = polyval(ptcoeff{i},t_mc_pf{i});
	disp(['Temperature fit : dtemp(ctd-mc) = A + B*temp_mc + C*temp_mc^2'])

	C = ptcoeff{i}(1); B = ptcoeff{i}(2); A = ptcoeff{i}(3);
	disp (['A = ',num2str(A),' B = ',num2str(B),' C  ',num2str(C)])

	% Now calculate corrected Microcat data
	% correct mc conductivity
	c_corr{i} = (c_mci{i} * nanmean(K{i})) .* bad_data{i}; 
	% correct mc pressure
	p_corr{i} = (p_mci{i} + (polyval(pcoeff{i},p_mci{i})));
	% correct mc temperature
	t_corr{i} = (t_mci{i} + (polyval(ptcoeff{i},t_mci{i})));

	% recompute differences with ctd to check corrections
	c_diff_corr{i} = (c_ctdi - c_corr{i}) .* bad_data{i};
	K_corr{i} = (c_ctdi ./ c_corr{i}).*bad_data{i};
	t_diff_corr{i} = (t_ctdi - t_corr{i}) .* bad_data{i};
	p_diff_corr{i} = (p_ctdi - p_corr{i}) .* bad_data{i};
end

	
		%%%%%%%%%%% FIGURE 1 %%%%%%%%%%%%%%%
for i = 1:n_mc
	% Figure(1) : show which portions of data are to be compared
	figure(1); orient tall;clf
	subplot(4,1,1),plot(TI,dpdt,'k')
	grid on; hold on
	plot(TI,dpdt.*bad_data{i},'r+')
	axis([min(TI) max(TI) -1.5 1.5]); ylabel('dp/dt (db/s)'); xlabel('time (s) from start of record')
	title(['Microcat serial number ', num2str(mc_cfg{i}.temp_sn)])

	subplot(4,1,2),plot(TI,p_ctdi,'ko-')
	hold on; grid on
	plot(TI,p_mci{i},'g-')
	plot(TI,p_mci{i} .* bad_data{i},'r+-')
	axis([min(TI) max(TI) 0 5000]); ylabel('press ctd(blk) mc(red)')

	subplot(4,1,3),plot(TI,t_ctdi,'ko-')
	hold on; grid on
	plot(TI,t_mci{i},'g-')
	plot(TI,t_mci{i} .* bad_data{i},'r+-')
	axis([min(TI) max(TI) 0 25]); ylabel('temp ctd(blk) mc(red)')

	subplot(4,1,4),plot(TI,c_ctdi,'ko-')
	grid on; hold on
	plot(TI,c_mci{i},'g-')
	plot(TI,c_mci{i} .* bad_data{i},'r+-')
	axis([min(TI) max(TI) 32 55]); ylabel('cond ctd(blk) mc(red)')
	%%%%%%%%%%% END FIGURE 1 %%%%%%%%%%%%%%%

	%%%%%%%%%%% FIGURE 2 %%%%%%%%%%%%%%%
	figure(2); orient tall;clf
	subplot(4,1,1),plot(c_ctdi,c_diff{i},'ok')
	axis([30 65 min(c_diff{i}) max(c_diff{i})]); grid on; hold on
	xlabel('ctd cond (mS/cm)'); ylabel('Cond Diff, ctd-mc (mS/cm)')
	title(['Microcat serial number ', num2str(mc_cfg{i}.temp_sn),...
		' : ', 'CTD Station ',sbe9_filename],'interpreter','none')
	plot(c_ctdi,nanmean(c_diff{i}),'r-');

	subplot(4,1,2),plot(p_ctdi,c_diff{i},'ok')
	axis([0 max(p_mci{i}) min(c_diff{i}) max(c_diff{i})]); grid on;   hold on
	xlabel('ctd press (dbar)'); ylabel('Cond Diff, ctd-mc (mS/cm)')
	plot(p_ctdi,nanmean(c_diff{i}),'r-')

	subplot(4,1,3),plot(p_mci{i},p_diff{i},'ok')
	axis([0 max(p_mci{i}) min(p_diff{i}) max(p_diff{i})]); grid on;   hold on
	xlabel('MC press (dbar)'); ylabel('Press diff, ctd-mc (dbar)')
	plot(p_mc_pf{i},p_fit{i},'r-')

	subplot(4,1,4),plot(t_ctdi,t_diff{i},'ok')
	axis([nanmin(t_ctdi) nanmax(t_ctdi) min(t_diff{i}) max(t_diff{i})]); grid on;   hold on
	xlabel('ctd temp ({^o}C)'); ylabel('Temp diff, ctd-mc (degC)')
	plot(t_ctd_pf{i},t_fit{i},'r')
	eval(['print -dpsc ',caldip_id,'Fig2a_sn',num2str(mc_cfg{i}.temp_sn),'.ps'])

	%%%%%%%%%%% FIGURE 2 %%%%%%%%%%%%%%%
	
	%%%%%%%%%%% FIGURE 3 %%%%%%%%%%%%%%%
	% Plot corrected data to check the calibrations
	figure(3); orient tall;clf
	subplot(4,1,1),plot(c_ctdi,K_corr{i},'ok')
	axis([32 55 min(K_corr{i}) max(K_corr{i})]);
	grid on; xlabel('ctd cond (mS/cm)'); ylabel('Cond diff, ctd-mc (mS/cm)')
	title(['Corrected : Microcat serial number ', num2str(mc_cfg{i}.temp_sn),...
		' : ', 'CTD Station ',sbe9_filename],'interpreter','none')

	subplot(4,1,2),plot(p_ctdi,K_corr{i},'ok')
	axis([0 max(p_mci{i}) min(K_corr{i}) max(K_corr{i})])
	grid on; xlabel('ctd press (dbar)'); ylabel('Cond diff, ctd-mc (mS/cm)')

	subplot(4,1,3),plot(p_corr{i},p_diff_corr{i},'ok')
	axis([0 max(p_mci{i}) min(p_diff_corr{i}) max(p_diff_corr{i})])
	grid on; xlabel('mc press (dbar)'); ylabel('Press diff, ctd-mc (dbar)')

	subplot(4,1,4),plot(t_ctdi,t_diff_corr{i},'ok')
	axis([nanmin(t_ctdi) nanmax(t_ctdi) min(t_diff_corr{i}) max(t_diff_corr{i})])
	grid on; xlabel('ctd press (dbar)'); ylabel('Temp diff, ctd-mc (degC)')


	%   bit more investigation as a function of time at the bottle stops
	figure(4);orient tall;clf
	subplot(3,1,1),plot(TI,p_diff_corr{i},'ko-')
	grid on 
	title([' Microcat serial number ', num2str(mc_cfg{i}.temp_sn),':  pdiff (ctd-mc) dbar'])
	subplot(3,1,2),plot(TI,t_diff_corr{i},'ro-')
	grid on
	title('tdiff (ctd-mc) ¡C')
	subplot(3,1,3),plot(TI,c_diff_corr{i},'go-')
	grid on
	title('cdiff (ctd-mc) mS/cm')
	
	%cd(sbe37_pathname)
	%for j=1:4
	%	eval(['print -dpsc -append ',caldip_id,'_sn',num2str(mc_cfg{i}.temp_sn),'.ps'])
	%end
end

diary off
