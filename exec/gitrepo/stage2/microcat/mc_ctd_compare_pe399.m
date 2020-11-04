% mc_ctd_comvpare_kn221.m
clear

cruise = 'pe399';
cast = '4'; % NB CHANGE BACK TO 1db DATA ON LINE 23-24
ctdnum = sprintf('%03d',str2num(cast));
doctd = 0 ; % whether to load and plot CTD data 

% setup paths 
basedir = ['/home/mstar/osnap/data/']; % KN221
ctddir = '/home/mstar/pe399/data/ctd/'; %PE399

% --- set paths for data input and output ---
inpath   = [basedir 'moor/proc_calib/' cruise '/cal_dip/microcat/cast' cast '/'];
infofile  = [basedir 'moor/proc_calib/' cruise '/cal_dip/cast',cast,'info.dat'];


% load ctd file

ctdinfile = [ctddir 'ctd_',cruise,'_',ctdnum,'_1hz']; % Fix for problem...
% % in 1db 005 ctd file
% % load ctdinfile ;
% eval(['load ',ctdinfile]);
% ctd=mload(ctdinfile,'/');
%   Time origin for jday
jd0 = julian(2015,1,0,0);
 [ctd h]=mload(ctdinfile,'time','press','temp1','cond1',' ','q');
 HH = h.data_time_origin(4)+h.data_time_origin(5)/60+h.data_time_origin(6)/3600;
 jtime=h.data_time_origin(1:3);
 jtime=[jtime HH];
 ctd.timeJ=julian(jtime)+ctd.time/86400-jd0;
 ctdtimetx = datestr(h.data_time_origin);
%  ptittxt = sprintf('Cast %s start %s CTD sensor set %s',cast,ctdtimetx,ctdsen);

% load infofile
[id,sn]= rodbload(infofile,'instrument:serialnumber');

% --- vector of serial numbers ---
ii = find(id >= 332 & id <= 337);
vec = sn(ii);

% load microcats
for i=1:length(vec);
% for i=1:2;
    disp(['Loading : ',num2str(sn(i))]);
    
    infile = ([inpath, 'cast', cast ,'_',sprintf('%4.4d',vec(i)),'.raw']) ;
    
    [yy,mm,dd,hh,c,t,p] = rodbload(infile,'yy:mm:dd:hh:c:t:p');
    
    jd = julian(yy,mm,dd,hh)-julian(2015,1,0,0); % jd relative to start of year
        
    % put CTD and MC on same time base
    TIstart = max(ctd.timeJ(1),jd(1)) ; TIend = min(ctd.timeJ(end),jd(end)) ;
    TI=TIstart:(10/86400):TIend ;
    
    ctd.pi = interp1(ctd.timeJ,ctd.press,TI) ;
    ctd.ti = interp1(ctd.timeJ,ctd.temp1,TI) ;
    ctd.ci = interp1(ctd.timeJ,ctd.cond1,TI) ;
%    ctd.alti = interp1(ctd.timeJ,ctd.altM,TI) ;
    
    pi = interp1(jd,p,TI) ;
    ti = interp1(jd,t,TI) ;
    ci = interp1(jd,c,TI) ;
    
    % compute differences of p, c, t [CTD - mc]
    dp = ctd.pi - pi ;
    dt = ctd.ti - ti ;
    dc = ctd.ci - ci ;
    
    
    
    % plot only values where CTD is at bottle stop
    ti = find(max(ctd.pi) == ctd.pi) ; % Roughly bottom of the down cast
    ctd.dpi = diff(ctd.pi) ;
    k = find(ctd.dpi < 1 & ctd.dpi > -1 & TI(2:end) > TI(ti)) ;
    
    figure(sn(i)*100); clf ;
    plot(TI(k),pi(k),'ko');
    grid on ; hold on ;
    plot(TI(k),ctd.pi(k),'r+')
    datetick('x')
    title(['Cast : ',num2str(ctdnum),' : MC s/n : ',num2str(vec(i))]) ;
eval(['print -djpeg ','f_',num2str(sn(i)*100),'.jpg'])
    
    %%
    figure(sn(i)); clf ;
    subplot(3,1,1),plot(TI(k),dp(k),'ko'); grid on ; hold on; datetick('x') ;
    plot([TI(1),TI(end)],[0,0],'k','linewidth',1)
    ylim([-10 10]);
    title(['Cast : ',num2str(ctdnum),' : MC s/n : ',num2str(vec(i))]) ; ylabel('dpress (db)')
    subplot(3,1,2),plot(TI(k),dt(k),'ko'); grid on ; hold on; datetick('x') ;
    plot([TI(1),TI(end)],[0,0],'k','linewidth',1);
    ylim([-0.01 0.01])
    ylabel('dtemp (�C)');
    subplot(3,1,3),plot(TI(k),dc(k),'ko'); grid on ; hold on ; datetick('x') ;
    plot([TI(1),TI(end)],[0,0],'k','linewidth',1);
    ylim([-0.01 0.01])
    ylabel('dcond (�C)')
eval(['print -djpeg ','f_',num2str(sn(i)),'.jpg'])

    
    %%
    figure(sn(i)*10); clf ;
    subplot(3,1,1),plot(ctd.pi(k),dp(k),'ko'); grid on ; hold on;
    plot([0,max(ctd.pi)],[0,0],'k','linewidth',1)
    ylim([-10 10]);
    title(['Cast : ',num2str(ctdnum),' : MC s/n : ',num2str(vec(i))]) ;ylabel('dpress (db)')
    subplot(3,1,2),plot(ctd.pi(k),dt(k),'ko'); grid on ; hold on;
    plot([0,max(ctd.pi)],[0,0],'k','linewidth',1);
    ylim([-0.01 0.01])
    ylabel('dtemp (�C)');
    subplot(3,1,3),plot(ctd.pi(k),dc(k),'ko'); grid on ; hold on ;
    plot([0,max(ctd.pi)],[0,0],'k','linewidth',1);
    ylim([-0.01 0.01])
    ylabel('dcond (mS/cm)')

eval(['print -djpeg ','f_',num2str(sn(i)*10),'.jpg'])

    
    
end





