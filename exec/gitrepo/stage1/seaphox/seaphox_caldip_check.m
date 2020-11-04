% seaphox_caldip_check reads seaphox caldip data from raw rodb file
% and compares data with the lowered CTD data
%
% LOH May 2017 - adapted from mc_caldip_check_dy078


clearvars -except MEXEC MEXEC_A MEXEC_G;
close all;

pathosnap = '/home/mstar/osnap';

cruise_ctd = 'dy078'; % used for ctd data dir
cruise = cruise_ctd; % used for microcat data
cruise2=cruise_ctd; % used for ctd data name
%   Time origin for jday
jd0 = julian(2018,1,0,0);

basedir = [pathosnap filesep 'data' filesep];
ctddir = [pathosnap filesep 'data' filesep cruise_ctd filesep]; 

% -----------------------------------------------------------------
% User supplied information
% Select cast number and choose which set of instruments from CTD
cast = input('Which cast number? ','s');
ctdsen = input('Which CTD sensors (1 or 2?) ','s');
ctdnum = sprintf('%03d',str2num(cast));

% --- set paths for data input and output ---
outpath   = [basedir 'moor/proc_calib/' cruise '/cal_dip/seaphox/cast' cast '/'];
infofile  = [basedir 'moor/proc_calib/' cruise '/cal_dip/cast',cast,'info.dat'];
%ctdinfile = [ctddir  'ctd_' cruise2 '_',ctdnum,'_raw.nc'];
ctdinfile = [ctddir  'ctd_' cruise2 '_',ctdnum,'_psal.nc'];



% ----------------- load CTD DATA   ----------------------------------
doctd = 1;
%  [d h]=mload(ctdinfile,'timeJ','prDM','t090C','c0S_slash_m',' ');
 [d h]=mload(ctdinfile,'time','press','temp1','cond1','temp2','temp','cond2','oxygen1','oxygen2','psal1','psal2','psal',' ','q');
 HH = h.data_time_origin(4)+h.data_time_origin(5)/60+h.data_time_origin(6)/3600;
 jtime=h.data_time_origin(1:3);
 jtime=[jtime HH];
 d.timeJ=julian(jtime)+d.time/86400-jd0;
 ctdtimetx = datestr(h.data_time_origin);
 ptittxt = sprintf('Cast %s start %s CTD sensor set %s',cast,ctdtimetx,ctdsen);

% Selecttime period that we will analyse
mxpi = max(d.press);
imp = d.press > mxpi-10;
tm1p = min(d.timeJ(imp))+0.25*(max(d.timeJ(imp))- min(d.timeJ(imp)));
tm2p = max(d.timeJ(imp))-0.25*(max(d.timeJ(imp))- min(d.timeJ(imp)));
imp2 = d.timeJ > tm1p & d.timeJ < tm2p;
meanctdpr = nanmean(d.press(imp2));

% CTD stats during period
meanctdpr = nanmean(d.press(imp2));
stdctdpr = nanstd(d.press(imp2));
ctd1_sal_mn = nanmean(d.cond1(imp2));
ctd1_sal_st = nanstd(d.cond1(imp2));
ctd2_sal_mn = nanmean(d.cond2(imp2));
ctd2_sal_st = nanstd(d.cond2(imp2));
ctd1_temp_mn = nanmean(d.temp1(imp2));
ctd1_temp_st = nanstd(d.temp1(imp2));
ctd2_temp_mn = nanmean(d.temp2(imp2));
ctd2_temp_st = nanstd(d.temp2(imp2));
ctd1_oxy_mn = nanmean(d.oxygen1(imp2));
ctd1_oxy_st = nanstd(d.oxygen1(imp2));
ctd2_oxy_mn = nanmean(d.oxygen2(imp2));
ctd2_oxy_st = nanstd(d.oxygen2(imp2));

% --- get mooring information from infofile ---
'here'
infofile
[zins,id,sn]= rodbload(infofile,'z:instrument:serialnumber');
'done'

[stdt,stti,endt,enti]= rodbload(infofile,'StartDate:StartTime:EndDate:EndTime');
jdstt = julian([stdt;stti(1)+stti(2)/60]')-jd0;
jdend = julian([endt;enti(1)+enti(2)/60]')-jd0;

% --- vector of serial numbers ---
ii = find(id == 375);
vec = sn(ii);
sn = vec;
id2=id(ii);
nvec = length(ii);
zmic = zins(ii);
mx_length= 0; 



% --- read data loop --
for i = 1:nvec
% display( [,num2str(vec(i)),])
   outfile = [outpath, 'cast', cast ,'_',sprintf('%0d',vec(i)),'.raw'];
% --- load rodb data ---
   [yy,mm,dd,hh,ph,ph_v,t,p,o2,s,batt1,batt2] = rodbload(outfile,'YY:MM:DD:HH:PH:PH_V:T:P:O:S:BATT1:BATT2');

%  if (i > 6 & i<=12)  lstr='--'; elseif i>12  lstr ='-.'; else lstr = '-'; end 
   if (i > 7 & i<=14)  lstr='--'; elseif i>14  lstr ='-.'; else lstr = '-'; end 
%   disp(['Checking ',num2str(i),': s/n:',num2str(sn(i))])
%    pause
% Time variable
   jd = julian(yy,mm,dd,hh)-jd0;
   if length(jd)>mx_length;mx_length=length(jd);end;

% Conversion now done in seaphox2rodb_O1.m
%    % conversion from ml/L to micromol/kg
%    % Molar Volume of O2 at STP = 22.391 L
%    % 1 micromol O2 = 0.022391 mL
%    % 1 ml/L = 44.661 micromol/L
%    % 1 micromol/kg = 1000/(sw_dens(s,t,p) micromol/L
%    o2_2 = o2.*44.661./(sw_dens(s,t,p)./1000);
%    
    eval(['jd' sprintf('%d',sn(i)) '= jd;']);
    eval(['s' sprintf('%d',sn(i)) '= s;']);
    eval(['t' sprintf('%d',sn(i)) '= t;']);
    eval(['p' sprintf('%d',sn(i)) '= p;']);
    eval(['o2' sprintf('%d',sn(i)) '= o2_2;']);
    eval(['ph' sprintf('%d',sn(i)) '= ph;']);    
    
   
end % for 


jd = ones(length(sn),mx_length)+nan;
s = ones(length(sn),mx_length)+nan;
t = ones(length(sn),mx_length)+nan;
p = ones(length(sn),mx_length)+nan;
o2 = ones(length(sn),mx_length)+nan;
ph = ones(length(sn),mx_length)+nan;

ctds = ones(length(sn),mx_length)+nan;
ctdt = ones(length(sn),mx_length)+nan;
ctdp = ones(length(sn),mx_length)+nan;
ctdo2 = ones(length(sn),mx_length)+nan;


for i = 1:length(vec)
    ll = eval(['length(jd' sprintf('%d',sn(i)) ')']);
    jd(i,1:ll)=eval(['jd' sprintf('%d',sn(i)) ';']);
    t(i,1:ll)=eval(['t' sprintf('%d',sn(i)) ';']);
    s(i,1:ll)=eval(['s' sprintf('%d',sn(i)) ';']);
    p(i,1:ll)=eval(['p' sprintf('%d',sn(i)) ';']);
    o2(i,1:ll)=eval(['o2' sprintf('%d',sn(i)) ';']);
    ph(i,1:ll)=eval(['ph' sprintf('%d',sn(i)) ';']);    
    % interp here if doctd
    if doctd
        ni = ~isnan(s(i,:)); ctds(i,ni)=interp1(d.timeJ, d.psal, jd(i,ni));
        ni = ~isnan(t(i,:)); ctdt(i,ni)=interp1(d.timeJ, d.temp, jd(i,ni));
        ni = ~isnan(p(i,:)); ctdp(i,ni)=interp1(d.timeJ, d.press, jd(i,ni));
        ni = ~isnan(o2(i,:)); ctdo2(i,ni)=interp1(d.timeJ, d.oxygen2, jd(i,ni));
    end
    
    
end

pdiff = ['mean p diff = ' num2str(nanmedian(p'-ctdp')) ...
    ' --- ' num2str(vec(i))];
sdiff = ['mean s diff = ' num2str(nanmedian(s'-ctds')) ...
    ' --- ' num2str(vec(i))];
tdiff = ['mean t diff = ' num2str(nanmedian(t'-ctdt')) ...
    ' --- ' num2str(vec(i))];
odiff = ['mean o2 diff = ' num2str(nanmedian(o2'-ctdo2')) ...
        ' --- ' num2str(vec(i))];

disp(pdiff)
disp(sdiff)
disp(tdiff)
disp(odiff)

outfig = [outpath, 'cast', cast ,'_all'];
%===========================
% figure with sampling (add also pH)

    figure
    sp1=subplot(5,1,1);
    hold on
    hl = plot(d.timeJ,d.press); set(hl(1),'color',[0 0 0])
    plot(jd,p,'-b')
    xlim([jdstt jdend]);
    ylabel('Pressure')
    
    title(['CAST ' cast ' Calibration Dip'])    
    
    sp2=subplot(5,1,2)   ;
    hl = plot(d.timeJ,d.temp1,d.timeJ,d.temp2); set(hl(1),'color',[0 0 0]);set(hl(2),'color',[.6 .6 .6]);
    hold on
    p1=plot(jd,t,'b');    
    legend([hl'],{'temp1','temp2'},'location','best')
    xlim([jdstt jdend]);
    ylabel('Temp.')
    
    sp3=subplot(5,1,3)  ;  
    hl = plot(d.timeJ,d.psal1,d.timeJ,d.psal2); set(hl(1),'color',[0 0 0]);set(hl(2),'color',[.6 .6 .6]);
    hold on
    p1=plot(jd,s,'b');    
    legend([hl'],{'sal1','sal2'},'location','best')
    xlim([jdstt jdend]);
    ylabel('Sal.')
     
    sp4=subplot(5,1,4);
    hl = plot(d.timeJ,d.oxygen1,d.timeJ,d.oxygen2); set(hl(1),'color',[0 0 0]);set(hl(2),'color',[.6 .6 .6]);    
    hold on
    plot(jd,o2,'b')
    legend([hl'],{'oxy1','oxy2'},'location','best')    
    xlim([jdstt jdend]);
    ylabel('Oxy.') 
    
    sp5=subplot(5,1,5);
    plot(jd,ph,'b')
    xlim([jdstt jdend]);
    ylabel('pH')     
    
    orient tall
    print(gcf,'-depsc',[outfig '_all.ps'])   
linkaxes([sp1 sp2 sp3 sp4 sp5],'x')

saveas(gcf,[outfig '_all.fig'],'fig')
