% basic preprocessing for microcat data
%   
% features
%      1. eliminate launching and recovery period
%      2. save data to rodb file
%      3. create data overview sheet
%
% uses timeaxis.m, auto_filt.m, julian.m 

% 11.01.01 Kanzow 
% 13.08.02 Kanzow : debugged
%
% 25.04.05          debugged      
  
% --- get moring information from infofile 

% jym 21 april 2005: Start with clean slate:
close all
clear all

moor     = 'eb2_6_200655';    % Mooring name
operator      = 'mpc';

plot_interval = [2006 12 02 12;   % start time of time axis on plot
		 2007 10 17 00];  % end time of time axis on plot

mc_id    = [333 337] ;         % microcat id

% -- set path for data input and output

%% inpath  = ['/local/users/pstar/data/moor/proc/',moor,'/microcat/'];
%% outpath = ['/local/users/pstar/data/moor/proc/',moor,'/microcat/'];
%% infofile =['/local/users/pstar/data/moor/proc/',moor,'/',moor,'info.dat'];


% inpath  = ['/noc/ooc/rpdmoc/rapid/data/moor/proc/',moor,'/microcat/'];
% outpath = ['/noc/ooc/rpdmoc/rapid/data/moor/proc/',moor,'/microcat/'];
% infofile =['/noc/ooc/rpdmoc/rapid/data/moor/proc/',moor,'/',moor,'info.dat'];

% mpc change paths for cruise d324:
inpath   = ['/data32/d324/rapid/data/moor/raw/d324/microcat/'];  
outpath  = ['/data32/d324/rapid/data/moor/proc/',moor,'/microcat/'];
infofile = ['/data32/d324/rapid/data/moor/proc/',moor,'/',moor,'info.dat'];

[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

ii = find(id >= mc_id(1) & id<=mc_id(2));
sn = sn(ii);
z  = z(ii);
id = id(ii);

%sn = sn(1:2); z = z(1:2); id = id(1:2);

 
[z,zx] = sort(z);  % sort instruments by their depth
sn     = sn(zx);
id     = id(zx);



fid_stat= fopen([outpath,'stage2_log'],'w');
% fprintf(fid_stat,'Processing steps taken by microcat_rdb2use.m:\n');
% jym 21 April 2005: Make info consistent with filename used.
fprintf(fid_stat,'Processing steps taken by microcat_raw2use.m:\n');
fprintf(fid_stat,'  1. eliminate launching and recovery period\n');
  fprintf(fid_stat,'  2. save data to rdb file\n');
fprintf(fid_stat,'\n Operated by:%s on %s\n',operator,datestr(clock)); 
fprintf(fid_stat,['        MicroCAT in Mooring ',moor,'\n\n\n']);
fprintf(fid_stat,'     ID    Depth   Start         End      Cycles  Spikes  Gaps   Mean     STD     Max     Min\n');          

% ---- despike paramenters

%T_range = [-15 +15];
%C_range = [-30 +30];  
%P_range = [-100 2000]; 

%dT_range = 18;  % accepted standard deviation envelope of adjecent T-values
%dC_range = 18;  % accepted standard deviation envelope of adjecent C-values 
%dP_range = 18;  % accepted standard deviation envelope of adjecent P-values 
%nloop    = 3;

dummy    = -9999;


%-----------------------------------------
% --- preprocessing loop -------------------
% ----------------------------------------

inst = 1;

jd_s  = julian(s_d(1),s_d(2),s_d(3),s_t(1)+s_t(2)/60);  % start time
jd_e  = julian(e_d(1),e_d(2),e_d(3),e_t(1)+e_t(2)/60);  % end time

for proc = 1 : length(sn),

  infile  = [inpath,moor,'_',sprintf('%4.4d',sn(proc)),'.raw'];
  if exist(infile)   > 0 
 
    rodbfile= [moor,'_',sprintf('%4.4d',sn(proc)),'.use']; 
    outfile = [outpath,rodbfile];

    inst = inst +1;
    [YY,MM,DD,HH,C,T,P] = rodbload(infile,'YY:MM:DD:HH:C:T:P');

    %------------------------------------------ 
    %----- cut off launching and recovery period
    %------------------------------------------
    disp('cut off launching and recovery period')
 
    jd               = julian(YY,MM,DD,HH);
    ii               = find(jd <= jd_e & jd >= jd_s );
    YY=YY(ii);MM=MM(ii);DD=DD(ii);HH=HH(ii);c=C(ii);t=T(ii);
    jd  = jd(ii); 
    if length(P) > 1   p = P(ii); end

    cycles     = length(ii);
    Start_Date = [YY(1) MM(1) DD(1)];
    Start_Time = HH(1);
    End_Date = [YY(cycles) MM(cycles) DD(cycles)];
    End_Time = HH(cycles);     

    %------------------------------------------
    %--- despike ------------------------------
    %------------------------------------------
 %   disp('ddspike')   

 %   [t,tdx,tndx] = ddspike(T,T_range,dT_range,nloop,'y',dummy); 
 %   [c,cdx,cndx] = ddspike(C,C_range,dC_range,nloop,'y',dummy); 
 %   if length(P) > 1 
 %     [p,pdx,pndx] = ddspike(P,P_range,dP_range,nloop,'y',dummy); 
 %   end

    % -----------------------------------------
    % ---  basic statistics -------------------
    % -----------------------------------------
 %   tstat = find(t ~= dummy);
 %   cstat = find(c ~= dummy);    
 %   tstat = t(tstat);
 %   cstat = c(cstat);
 
 % JYM SCU 20/4/05 change to std matlab
 %   tm = meannan(t);
 %   cm = meannan(c);
  
 %  tsd= stdmiss(t);
 %  csd= stdmiss(c);
     tm = nanmean(t);
     cm = nanmean(c);
     tsd = nanstd(t);
     csd = nanstd(c);
 
 %   tmx = max(t);
 %   cmx = max(c);
 %   tmn = min(t);
 %   cmn = min(c);
     tmx = nanmax(t);
     cmx = nanmax(c);
     tmn = nanmin(t);
     cmn = nanmin(c);

 
    if length(P) > 1
 %     pstat = find(p ~= dummy);
 %     pm  = meannan(p);
 %     psd = stdmiss(p);
 %     pmx = max(p);
 %     pmn = min(p);
       pm = nanmean(p);
       psd = nanstd(p);
       pmx = nanmax(p);
       pmn = nanmin(p);
  % end JYM SCU changes
    end 
     
    %------------------------------------------
    %---- fill time gaps  with dummy
    %------------------------------------------

    disp(' fill time gaps  with dummy')

    djd = diff(jd);           % time step  
    sr  = median(djd);        % sampling interval
    ii  = find(djd > 1.5*sr);  % find gaps
    gap = round(djd(ii)/sr)-1;
    addt= []; 

    for i = 1 : length(gap), 
      addt = [addt; [[1:gap(i)]*sr + jd(ii(i))]'];
                         
    end 

    [jd,xx] = sort([jd; addt]);   % add time
    ngap    = length(addt);       % number of time gaps         
    gt      = gregorian(jd);
    YY=gt(:,1); MM=gt(:,2); DD=gt(:,3); 
    if size(gt,2) == 6
       HH=hms2h(gt(:,4),gt(:,5),gt(:,6)); 
    else 
       HH= gt(:,4);
    end    
       
   
    t = [t;dummy*ones(ngap,1)]; t = t(xx);
    c = [c;dummy*ones(ngap,1)]; c = c(xx); 
    if length(P) > 1
       p = [p;dummy*ones(ngap,1)]; p = p(xx); 
    end
    %-----------------------------------------------------
    %  write output to logfile ---------------------------
    %-----------------------------------------------------

    disp(' write output to logfile')


    fprintf(fid_stat,'T   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d   %d         %d   %5.2f   %5.2f   %5.2f   %5.2f \n',...
               sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,tm,tsd,tmx,tmn'); 

    fprintf(fid_stat,'C   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d   %d         %d   %5.2f   %5.2f   %5.2f   %5.2f \n',...
               sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,cm,csd,cmx,cmn'); 

    if length(P) > 1
      fprintf(fid_stat,'P   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d  %d        %d    %5.1f   %5.2f   %5.2f   %5.2f \n',...
               sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,pm,psd,pmx,pmn');  
    end
    fprintf(fid_stat,'\n');

    %-----------------------------------  
    %--- write data to rodb format -----
    %-----------------------------------

    disp(['writing data to ',outfile]) 
         
    if length(P) <= 1
       sub =2;
       fort = '%4.4d   %2.2d   %2.2d   %8.5f   %6.4f   %6.4f';
       cols = 'YY:MM:DD:HH:T:C';
       rodbsave(outfile,...
         'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
          fort,...
	  lat,lon,cols,s_d,s_t,sn(proc),mr,wd,z(proc),e_d,e_t,...
          [YY MM DD HH t c]);
    else
      sub  = 3;  
      fort = '%4.4d   %2.2d   %2.2d   %8.5f   %6.4f   %6.4f  %5.1f';
       cols = 'YY:MM:DD:HH:T:C:P';
      rodbsave(outfile,...
          'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
           fort,...
	  lat,lon,cols,s_d,s_t,sn(proc),mr,wd,z(proc),e_d,e_t,...
          [YY MM DD HH t c p]);
    end

  %%%%%%%%%% Graphics %%%%%%%%%%%%%%%%

       jd1 = julian(plot_interval(1,:));
       jd2 = julian(plot_interval(2,:)); 

   figure(1);clf
     subplot(sub,1,1); ii = find(~isnan(t)&t>dummy);

       plot(jd(ii)-jd1,t(ii))
         title(['MicroCAT s/n: ',num2str(sn(proc)), ...
                '; Target Depth: ',num2str(z(proc))])
         ylabel('Temperature [deg C]')
         grid on
         xlim([0 jd2-jd1])
	 timeaxis(plot_interval(1,1:3));   
 
    subplot(sub,1,2); ii = find(~isnan(c)&c>dummy);

       plot(jd(ii)-jd1,c(ii))
	 ylabel('Conductivity [mS/cm]')
         grid on
         xlim([0 jd2-jd1])
	 timeaxis(plot_interval(1,1:3));   

     if sub == 3 

      subplot(sub,1,3); ii = find(~isnan(p)&p>dummy);

       plot(jd(ii)-jd1,p(ii))

         ylabel('Pressure [dbar]')
         grid on 
         xlim([0 jd2-jd1])
	 timeaxis(plot_interval(1,1:3));   

     end
       eval(['print -dps ',outfile,'.ps']) 

  sampling_rate = 1/median(diff(jd));
  innan1        = find(~isnan(t));             
  tf            = auto_filt(t(innan1), sampling_rate, 1/2,'low',4);
  innan2        = find(~isnan(c));             
  cf            = auto_filt(c(innan2), sampling_rate, 1/2,'low',4);
  innan3        = find(~isnan(p));             
  pf            = auto_filt(p(innan3), sampling_rate, 1/2,'low',4);

  figure(2);clf
     subplot(sub,1,1); 

      plot(jd(innan1)-jd1,tf)
         title(['2-day low-pass; MicroCAT s/n: ', ...
               num2str(sn(proc)),'; Target Depth: ',num2str(z(proc))])
         ylabel('Temperature [deg C]')
         grid on  
         xlim([0 jd2-jd1])
	 timeaxis(plot_interval(1,1:3));   
    
     subplot(sub,1,2); ii = find(~isnan(c)&c>dummy);
 
       plot(jd(innan2)-jd1,cf)
	 ylabel('Conductivity [mS/cm]')
         grid on
         xlim([0 jd2-jd1])
	 timeaxis(plot_interval(1,1:3));   

    if sub == 3 
           
       subplot(sub,1,3)
 
         plot(jd(innan3)-jd1,pf)
	   ylabel('Pressure [dbar]')
           grid on 
           xlim([0 jd2-jd1])
           timeaxis(plot_interval(1,1:3));   

    end
      eval(['print -dps ',outfile,'_lowpass.ps']) 
 
  

  end

 


end
       
  
