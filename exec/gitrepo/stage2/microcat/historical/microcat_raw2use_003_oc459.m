% MICROCAT_RAW2USE_003 is a script that performs stage2 processing% on microcat data:%      1. eliminates lauch and recovery periods%      2. saves data to rodb file%      3. creates data overview figures%% It calls rodbload.m, rodbsave.m, timeaxis.m, auto_filt.m,% julian.m, ddspike.m % 11.01.01 Kanzow % 13.08.02 Kanzow : debugged% 18.10.09 P Wright: changed for D344% 22.03.2010 ZB Szuts: modified for Oceanus 459%      clear allclose all% -----------------------------------------------------------------% --- This is the information that needs to be modified for -------% --- different users, directory trees, and moorings --------------% -----------------------------------------------------------------% the location where the processing is done (a cruise name, NOCS, etc)location = 'oc459';operator = 'zszuts';%moor     = 'wb1_6_200906';moor     = 'wbh2_3_200912';[out,host] = unix('hostname');if strfind(host,'dhcp108.science.oceanus.whoi.edu');  basedir  = '/Users/surman/rpdmoc/rapid/data/';else  basedir  = '/Volumes/surman/rpdmoc/rapid/data/';end% the start and end times of the time axis for plottingplot_interval = [2009 4 01 0;		 2010 4 15 0];% -----------------------------------------------------------------% --- set paths for data input and output ---inpath   = [basedir 'moor/proc/' moor '/microcat/'];outpath  = [basedir 'moor/proc/' moor '/microcat/'];infofile = [basedir 'moor/proc/' moor '/' moor 'info.dat'];mc_id    = [333 337] ;             % microcat id numbers% --- get mooring information from infofile ---[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr] = ...    rodbload(infofile,['instrument:serialnumber:z:Start_Time:Start_Date:'...                    'End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring']);ii = find(id >= mc_id(1) & id<=mc_id(2));sn = sn(ii);z  = z(ii);id = id(ii);%sn = sn(1:2); z = z(1:2); id = id(1:2); % --- sort instruments by their depth ---[z,zx] = sort(z);sn     = sn(zx);id     = id(zx);% --- write header info to log file ---fid_stat = fopen([outpath,'stage2_log'],'w');fprintf(fid_stat,'Processing steps taken by microcat_rdb2use.m:\n');fprintf(fid_stat,'  1. eliminate lauching and recovery period\n');fprintf(fid_stat,'  2. save data to rdb file\n');fprintf(fid_stat,'\n Operated by:%s on %s at %s\n',...        operator,location,datestr(clock)); fprintf(fid_stat,['        MicroCAT in Mooring ',moor,'\n\n\n']);fprintf(fid_stat,'     ID    Depth   Start         End      Cycles  Spikes  Gaps   Mean     STD     Max     Min\n');          % --- despike paramenters ---dodespike = 0; % whether to despike the instrumentsif dodespike  T_range = [-15 +15];  C_range = [-30 +30];    P_range = [-100 2000];   dT_range = 18;  % accepted standard deviation envelope of adjecent T-values  dC_range = 18;  % accepted standard deviation envelope of adjecent C-values   dP_range = 18;  % accepted standard deviation envelope of adjecent P-values   nloop    = 3;enddummy    = -9999;% ------------------------------------------% --- preprocessing loop -------------------% ------------------------------------------jd_s  = julian(s_d(1),s_d(2),s_d(3),s_t(1)+s_t(2)/60);  % start timejd_e  = julian(e_d(1),e_d(2),e_d(3),e_t(1)+e_t(2)/60);  % end timefor proc = 1 : length(sn),  infile  = [inpath,moor,'_',sprintf('%4.4d',sn(proc)),'.raw'];  if exist(infile)     rodbfile= [moor,'_',sprintf('%4.4d',sn(proc)),'.use'];    outfile = [outpath,rodbfile];    % --- load data ---    [YY,MM,DD,HH,C,T,P] = rodbload(infile,'YY:MM:DD:HH:C:T:P');    % -----------------------------------------------    % ----- cut off launching and recovery period ---    % -----------------------------------------------    disp('cut off launching and recovery period')     jd = julian(YY,MM,DD,HH);    ii = find(jd <= jd_e & jd >= jd_s );    YY = YY(ii);   MM = MM(ii);   DD = DD(ii);    HH = HH(ii);   c = C(ii);     t = T(ii);    jd = jd(ii);     if length(P) > 1,  p = P(ii); end    cycles     = length(ii);    Start_Date = [YY(1) MM(1) DD(1)];    Start_Time = HH(1);    End_Date = [YY(cycles) MM(cycles) DD(cycles)];    End_Time = HH(cycles);         % ------------------------------------------    % --- despike ------------------------------    % ------------------------------------------    if dodespike      disp('ddspike')         [t,tdx,tndx] = ddspike(T,T_range,dT_range,nloop,'y',dummy);       [c,cdx,cndx] = ddspike(C,C_range,dC_range,nloop,'y',dummy);       if length(P) > 1         [p,pdx,pndx] = ddspike(P,P_range,dP_range,nloop,'y',dummy);       end      % -----------------------------------------      % ---  basic statistics -------------------      % -----------------------------------------      tstat = t(find(t ~= dummy));      cstat = c(find(c ~= dummy));        end    tm = meannan(t);    cm = meannan(c);      tsd= stdnan(t);    csd= stdnan(c);    tmx = max(t);    cmx = max(c);    tmn = min(t);    cmn = min(c);    if length(P) > 1      % pstat = find(p ~= dummy);      pm  = meannan(p);      psd = stdnan(p);      pmx = max(p);      pmn = min(p);    end          % ------------------------------------------    % ---- fill time gaps  with dummy    % ------------------------------------------    disp(' fill time gaps  with dummy')    djd = diff(jd);           % time step      sr  = median(djd);        % sampling interval    ii  = find(djd > 1.5*sr);  % find gaps    gap = round(djd(ii)/sr)-1;    addt= [];     for i = 1 : length(gap),       addt = [addt; [[1:gap(i)]*sr + jd(ii(i))]'];                             end     [jd,xx] = sort([jd; addt]);   % add time    ngap    = length(addt);       % number of time gaps             gt      = gregorian(jd);    YY = gt(:,1);   MM = gt(:,2);   DD = gt(:,3);     if size(gt,2) == 6       HH=hms2h(gt(:,4),gt(:,5),gt(:,6));     else        HH= gt(:,4);    end               t = [t; dummy*ones(ngap,1)]; t = t(xx);    c = [c; dummy*ones(ngap,1)]; c = c(xx);     if length(P) > 1       p = [p; dummy*ones(ngap,1)]; p = p(xx);     end        % -----------------------------------------------------    % --- write output to logfile -------------------------    % -----------------------------------------------------    disp(' write output to logfile')    fprintf(fid_stat,'T   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d   %d         %d   %5.2f   %5.2f   %5.2f   %5.2f \n',...               sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,tm,tsd,tmx,tmn');     fprintf(fid_stat,'C   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d   %d         %d   %5.2f   %5.2f   %5.2f   %5.2f \n',...               sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,cm,csd,cmx,cmn');     if length(P) > 1      fprintf(fid_stat,'P   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d  %d        %d    %5.1f   %5.2f   %5.2f   %5.2f \n',...               sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,pm,psd,pmx,pmn');      end    fprintf(fid_stat,'\n');        %-----------------------------------      %--- write data to rodb format -----    %-----------------------------------    disp(['writing data to ',outfile])         rodboutvars = ['Latitude:Longitude:Columns:Start_Date:Start_Time:'...                   'SerialNumber:Mooring:WaterDepth:Instrdepth:'...                   'End_Date:End_Time'];    if length(P) <= 1      sub  = 2;      fort = '%4.4d   %2.2d   %2.2d   %8.5f   %6.4f   %6.4f';      cols = 'YY:MM:DD:HH:T:C';      rodbsave(outfile,rodboutvars,fort,...               lat,lon,cols,Start_Date,Start_Time,sn(proc),mr,wd,...               z(proc),End_Date,End_Time,[YY MM DD HH t c]);    else      sub  = 3;        fort = '%4.4d   %2.2d   %2.2d   %8.5f   %6.4f   %6.4f  %5.1f';      cols = 'YY:MM:DD:HH:T:C:P';      rodbsave(outfile,rodboutvars,fort,...               lat,lon,cols,Start_Date,Start_Time,sn(proc),mr,wd,...               z(proc),End_Date,End_Time,[YY MM DD HH t c p]);    end    %-----------------------------------      %--- make graphics -----------------    %-----------------------------------    jd1 = julian(plot_interval(1,:));    jd2 = julian(plot_interval(2,:));     % --- plot raw data ---    figure(1); clf, orient tall    subplot(sub,1,1); ii = find(~isnan(t) & t>dummy);        plot(jd(ii)-jd1,t(ii))    title(['MicroCAT s/n: ',num2str(sn(proc)),'; Target Depth: ',num2str(z(proc))])    ylabel('Temperature [deg C]')    grid on    xlim([0 jd2-jd1])    timeaxis(plot_interval(1,1:3));           subplot(sub,1,2);    ii = find(~isnan(c) & c>dummy);    plot(jd(ii)-jd1,c(ii))    ylabel('Conductivity [mS/cm]')    grid on    xlim([0 jd2-jd1])    timeaxis(plot_interval(1,1:3));             if sub == 3       subplot(sub,1,3); ii = find(~isnan(p)&p>dummy);      plot(jd(ii)-jd1,p(ii))      ylabel('Pressure [dbar]')      grid on       xlim([0 jd2-jd1])      timeaxis(plot_interval(1,1:3));       end    eval(['print -dps ',outfile,'.ps'])         % --- filter raw data ---    % set up empty matrices, in case there are dummy values    tf  = repmat(nan,size(t));    cf  = repmat(nan,size(t));    pf  = repmat(nan,size(t));    sal = repmat(nan,size(t));    sf  = repmat(nan,size(t));        sampling_rate = 1/median(diff(jd));    % remove dummy values for filtering, ZBS    % --> this only works if there are isolated spikes!    ii = find(isfinite(t) & t~=dummy);    tf(ii) = auto_filt(t(ii), sampling_rate, 1/2,'low',4);    ii = find(isfinite(c) & c~=dummy);    cf(ii) = auto_filt(c(ii), sampling_rate, 1/2,'low',4);    ii = find(isfinite(p) & p~=dummy);    pf(ii) = auto_filt(p(ii), sampling_rate, 1/2,'low',4);    ii = find(isfinite(c) & c~=dummy);    sal(ii) = sw_salt(c(ii)/sw_c3515,t(ii),p(ii));    sf(ii) = auto_filt(sal(ii), sampling_rate, 1/2,'low',4);    % give a warning if there are large gaps that are filtered over    if any(diff(ii)>5)      warning(['microcat_raw2use: exist gaps (of dummy values) '...               'that are longer than 5 data points'])      fprintf(fid_stat,...              ['microcat_raw2use: exist gaps (of dummy values) '...               'that are longer than 5 data points']);    end        % --- plot filtered data ---    figure(2); clf, orient tall    subplot(sub,1,1); ii = find(~isnan(t)&t>dummy);    plot(jd-jd1,tf)    title(['2-day low-pass; MicroCAT s/n: ',num2str(sn(proc)),...           '; Target Depth: ',num2str(z(proc))])    ylabel('Temperature [deg C]')    grid on      xlim([0 jd2-jd1])    timeaxis(plot_interval(1,1:3));       if 1==0 % plot conductivity only      subplot(sub,1,2); ii = find(~isnan(c)&c>dummy);      plot(jd-jd1,cf)      ylabel('Conductivity [mS/cm]')      grid on      xlim([0 jd2-jd1])      timeaxis(plot_interval(1,1:3));               elseif 1==1 % plot conductivity and salinity together      subplot(sub,1,2);      ii = find(isfinite(cf+sf) & c>dummy);      [ha,hc,hs] = plotyy(jd(ii)-jd1,cf(ii),jd(ii)-jd1,sf(ii));      set(ha(1),'ycolor','k') % cond in black      set(hc,'color','k')      set(hs,'linestyle','-','color',[1 1 1]*0.5) % sal in grey      set(ha(2),'ycolor',[1 1 1]*0.5) % sal in grey      set(ha(1),'ylim',prctile(c(ii),[1 99])+[-1 1]*0.02,'ytickmode','auto')      set(ha(2),'ylim',prctile(sal(ii),[1 99]+[-1 1]*0.02),'ytickmode','auto')      ylabel(ha(1),'Conductivity (-) [mS/cm]')      ylabel(ha(2),'Salinity (grey) [psu]')      grid on       xlim(ha(1),[0 jd2-jd1])      xlim(ha(2),[0 jd2-jd1])      timeaxis(ha(1),plot_interval(1,1:3));         set(ha(2),'xtick',[]), linkaxes(ha,'x')    end    if sub == 3       subplot(sub,1,3)      plot(jd-jd1,pf)      ylabel('Pressure [dbar]')      grid on       xlim([0 jd2-jd1])      timeaxis(plot_interval(1,1:3));       end    eval(['print -dpsc ',outfile,'_lowpass.ps'])     %pause      end % if exist(infile)  end % for proc = 1 : length(sn),comment = input('Enter additional comment to be save in Log file: ','s'); if ~isempty(comment)  fprintf(fid_stat,'\n COMMENT:\n %s',comment);endfclose(fid_stat);