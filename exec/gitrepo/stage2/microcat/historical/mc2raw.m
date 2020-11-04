% basic preprocessing for microcat data%   % features%      1. eliminate lauching and recovery period%      2. despike%      3. save data to rodb file%      4. create data overview sheet%% uses ddspike.m (Kanzow)% 11.01.01 Kanzow % 13.08.02 Kanzow : debugged      % --- get moring information from infofile mc_id    = [333 337] ;         % microcat idmoor     = 'eb1_1_200409';    % Mooring nameoperator = 'Kanzow';% -- set path for data input and outputinpath  = ['/data/rapid/cd170/moorings/',moor,'/microcat/'];outpath = ['/data/rapid/cd170/moorings/',moor,'/microcat/'];infofile =['/data/rapid/cd170/moorings/',moor,'/',moor,'info.dat'];[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');ii = find(id >= mc_id(1) & id<=mc_id(2));sn = sn(ii);z  = z(ii);id = id(ii);%sn = sn(1:2); z = z(1:2); id = id(1:2); [z,zx] = sort(z);  % sort instruments by their depthsn     = sn(zx);id     = id(zx);fid_stat= fopen([outpath,'stage2_log'],'w');fprintf(fid_stat,'Processing steps taken by mc2raw.m:\n');fprintf(fid_stat,'  1. eliminate lauching and recovery period\n');  fprintf(fid_stat,'  2. despike using ddspike.m\n');  fprintf(fid_stat,'  3. save data to rodb file\n');fprintf(fid_stat,'\n Operated by:%s on %s\n',operator,datestr(clock)); fprintf(fid_stat,['        MicroCAT in Mooring ',moor,'\n\n\n']);fprintf(fid_stat,'     ID    Depth   Start         End      Cycles  Spikes  Gaps   Mean     STD     Max     Min\n');          % ---- despike paramentersT_range = [-15 +15];C_range = [-30 +30];  P_range = [-100 2000]; dT_range = 18;  % accepted standard deviation envelope of adjecent T-valuesdC_range = 18;  % accepted standard deviation envelope of adjecent C-values dP_range = 18;  % accepted standard deviation envelope of adjecent P-values nloop    = 3;dummy    = -9999;%-----------------------------------------% --- preprocessing loop -------------------% ----------------------------------------inst = 1;jd_s  = julian(s_d(1),s_d(2),s_d(3),s_t(1)+s_t(2)/60);  % start timejd_e  = julian(e_d(1),e_d(2),e_d(3),e_t(1)+e_t(2)/60);  % end timefor proc = 1 : length(sn),  infile  = [inpath,moor,'_',sprintf('%4.4d',sn(proc)),'.raw'];  if exist(infile)   > 0      rodbfile= [moor,'_',sprintf('%3.3d',inst),'.use']; %    infile  = [inpath,moor,'_',sprintf('%4.4d',sn(proc)),'.mc'];    outfile = [outpath,rodbfile];    inst = inst +1;    [YY,MM,DD,HH,C,T,P] = rodbload(infile,'YY:MM:DD:HH:C:T:P');    %------------------------------------------     %----- cut off launching and recovery period    %------------------------------------------    disp('cut off launching and recovery period')     jd               = julian(YY,MM,DD,HH);    ii               = find(jd <= jd_e & jd >= jd_s );    YY=YY(ii);MM=MM(ii);DD=DD(ii);HH=HH(ii);C=C(ii);T=T(ii);    jd  = jd(ii);     if length(P) > 1   P = P(ii); end    cycles     = length(ii);    Start_Date = [YY(1) MM(1) DD(1)];    Start_Time = HH(1);    End_Date = [YY(cycles) MM(cycles) DD(cycles)];    End_Time = HH(cycles);         %------------------------------------------    %--- despike ------------------------------    %------------------------------------------    disp('ddspike')       [t,tdx,tndx] = ddspike(T,T_range,dT_range,nloop,'y',dummy);     [c,cdx,cndx] = ddspike(C,C_range,dC_range,nloop,'y',dummy);     if length(P) > 1       [p,pdx,pndx] = ddspike(P,P_range,dP_range,nloop,'y',dummy);     end    % -----------------------------------------    % ---  basic statistics -------------------    % -----------------------------------------    tstat = find(t ~= dummy);    cstat = find(c ~= dummy);        tstat = t(tstat);    cstat = c(cstat);    tm = meannan(tstat);    cm = meannan(cstat);      tsd= stdmiss(tstat);    csd= stdmiss(cstat);    tmx = max(tstat);    cmx = max(cstat);    tmn = min(tstat);    cmn = min(cstat);    if length(P) > 1      pstat = find(p ~= dummy);      pm  = meannan(pstat);      psd = stdmiss(pstat);      pmx = max(pstat);      pmn = min(pstat);    end          %------------------------------------------    %---- fill time gaps  with dummy    %------------------------------------------    disp(' fill time gaps  with dummy')    djd = diff(jd);           % time step      sr  = median(djd);        % sampling interval    ii  = find(djd > 1.5*sr);  % find gaps    gap = round(djd(ii)/sr)-1;    addt= [];     for i = 1 : length(gap),       addt = [addt; [[1:gap(i)]*sr + jd(ii(i))]'];                             end     [jd,xx] = sort([jd; addt]);   % add time    ngap    = length(addt)       % number of time gaps             gt      = gregorian(jd);    YY=gt(:,1); MM=gt(:,2); DD=gt(:,3);     if size(gt,2) == 6       HH=hms2h(gt(:,4),gt(:,5),gt(:,6));     else        HH= gt(:,4);    end                  t = [t;dummy*ones(ngap,1)]; t = t(xx);    c = [c;dummy*ones(ngap,1)]; c = c(xx);     if length(P) > 1       p = [p;dummy*ones(ngap,1)]; p = p(xx);     end    %-----------------------------------------------------    %  write output to logfile ---------------------------    %-----------------------------------------------------    disp(' write output to logfile')    fprintf(fid_stat,'T   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d   %d     %d      %d   %5.2f   %5.2f   %5.2f   %5.2f \n',...               sn(proc),z(proc),Start_Date,End_Date,cycles,tndx,ngap,tm,tsd,tmx,tmn');     fprintf(fid_stat,'C   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d   %d     %d      %d   %5.2f   %5.2f   %5.2f   %5.2f \n',...               sn(proc),z(proc),Start_Date,End_Date,cycles,cndx,ngap,cm,csd,cmx,cmn');     if length(P) > 1      fprintf(fid_stat,'P   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d  %d     %d      %d    %5.1f   %5.2f   %5.2f   %5.2f \n',...               sn(proc),z(proc),Start_Date,End_Date,cycles,pndx,ngap,pm,psd,pmx,pmn');      end    fprintf(fid_stat,'\n');    %-----------------------------------      %--- write data to rodb format -----    %-----------------------------------    disp(['writing data to ',outfile])              if length(P) <= 1       sub =2;       fort = '%4.4d   %2.2d   %2.2d   %8.5f   %6.4f   %6.4f';       cols = 'YY:MM:DD:HH:T:C';       rodbsave(outfile,...         'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...          fort,...          lat,lon,cols,Start_Date,Start_Time,sn(proc),mr,wd,z(proc),End_Date,End_Time,...          [YY MM DD HH t c]);    else      sub  = 3;        fort = '%4.4d   %2.2d   %2.2d   %8.5f   %6.4f   %6.4f  %5.1f';       cols = 'YY:MM:DD:HH:T:C:P';      rodbsave(outfile,...          'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...           fort,...          lat,lat,cols,Start_Date,Start_Time,sn(proc),mr,wd,z(proc),End_Date,End_Time,...          [YY MM DD HH t c p]);    end  %%%%%%%%%% Graphics %%%%%%%%%%%%%%%%   figure(1);clf     subplot(sub,1,1); ii = find(~isnan(t)&t>dummy);       plot(jd(ii)-jd(1),t(ii))       title(['Instrument: ',num2str(sn(proc)),'; Target Depth: ',num2str(z(proc))])       ylabel('Temperature [deg C]')       grid on     subplot(sub,1,2); ii = find(~isnan(c)&c>dummy);       plot(jd(ii)-jd(1),c(ii))	 ylabel('Conductivity [mS/cm]')      grid on     if sub < 3        xlabel(['Time [days] since ',datestr(gregorian(jd(1)))])     else        subplot(sub,1,3); ii = find(~isnan(p)&p>dummy);        plot(jd(ii)-jd(1),p(ii))	ylabel('Pressure [dbar]')        xlabel(['Time [days] since ',datestr(gregorian(jd(1)))])        grid on      end       eval(['print -dps ',outfile,'.ps'])   sampling_rate = 1/median(diff(jd));  tf            = auto_filt(t, sampling_rate, 1/2,'low',4);  cf            = auto_filt(c, sampling_rate, 1/2,'low',4);  pf            = auto_filt(p, sampling_rate, 1/2,'low',4);  figure(2);clf     subplot(sub,1,1); ii = find(~isnan(t)&t>dummy);       plot(jd-jd(1),tf)       title(['2-day low-pass; Instrument: ',num2str(sn(proc)),'; Target Depth: ',num2str(z(proc))])       ylabel('Temperature [deg C]')       grid on           subplot(sub,1,2); ii = find(~isnan(c)&c>dummy);       plot(jd-jd(1),cf)	 ylabel('Conductivity [mS/cm]')      grid on       if sub <3          xlabel(['Time [days] since ',datestr(gregorian(jd(1)))])           else                subplot(sub,1,3)       plot(jd-jd(1),pf)	ylabel('Pressure [dbar]')        xlabel(['Time [days] since ',datestr(gregorian(jd(1)))])        grid on      end      eval(['print -dps ',outfile,'_lowpass.ps'])      end end         