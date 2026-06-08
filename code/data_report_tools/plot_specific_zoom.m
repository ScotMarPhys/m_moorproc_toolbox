
close all
clearvars  -except pathosnap 
currentdir = pwd;


%-------------------------------------
% Definition of the different path
basedir  = [pathosnap '/data/moor/'];%'/home/sa02lh/Data/Dropbox/testmoorOSNAP/'; %
procpath = [basedir 'proc/'];
reportdir = [pathosnap '/Documents/datareports/'];%'/home/sa02lh/Data/Dropbox/testmoorOSNAP/proc/datareports/'; % /media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap/Documents/datareports/
outpathstats  = [reportdir 'stats/'];		
outpathfigs = [reportdir 'figs/'];
outpathfigszoom = '/Users/locupe/Dropbox/Work/Dataproc/OSNAP_mooring/moorOSNAP_quick_analyses/proc_plot/datareports/figszoom/'; 
%-------------------------------------
% Selection of the deployment year
% depyear ='01_2014';
% depyear ='02_2015';
depyear ='03_2016';
%-------------------------------------
% Selection of the mooring to process
moorlist ={'nocm1'};
% moorlist ={'nocm1','nocm2','nocm3','nocm4','nocm5'};
% moorlist = {'rteb1'}; %{'rteb1','rtwb1','rtwb2'};
%moorlist = {'rtadcp1'};%
 
proclvl = '3';%  If proclvl = '2': processing is using the .use files; proclvl = '3' : processing is using the .microcat and .edt files

%-------------------------------------------------------
% Specific period (e.g. SCVs)
% intervalstr = '[2014 08 05 00; 2014 09 05 00]'; % SCV?
%intervalstr = '[2018 01 01 00; 2018 05 01 00]'
intervalstr = [[2018 02 07 00]; [2018 02 09 00]];

for ijk= 1:length(moorlist)
	timelimit=intervalstr;
	strtime = [datestr(datenum(timelimit(1,1:3)),'yyyymmdd') '_' datestr(datenum(timelimit(2,1:3)),'yyyymmdd')];
	eval(['cd ' outpathfigszoom]);	  
	close all
	moor=[moorlist{ijk} '_' depyear];
	if exist([outpathfigszoom strtime filesep moor])~=7
	   eval(['!/bin/mkdir -p ' strtime filesep moor]);
	   cd([strtime filesep moor]);
	else
	   cd([strtime filesep moor]);
    end
    
   %--------------------------------------- 
   % i=1; % top microcat
   % infile = [procpath,moor,'/microcat/',moor,'_',sprintf('%0.3d',i),'.microcat'];
   i ='3218'
   infile = [procpath,moor,'/microcat/',moor,'_',i,'.use'];
   % read data into vectors and then into structure array

   [yy,mm,dd,hh,t,c,pmc] = ...
       rodbload(infile,'yy:mm:dd:hh:t:c:p');
   jd=julian(yy,mm,dd,hh);

   bad_data=find(c==-9999); c(bad_data)=NaN;
   bad_data=find(t==-9999); t(bad_data)=NaN;

   %compute salinity from standard seawater routines
   c3515=sw_c3515;
   c_ratio=c/c3515;
   %t_68=t68(t);
   %s=sw_salt(c_ratio,t_68,p);
   s=sw_salt(c_ratio,t,pmc); % version 3.2 of the Sea water routines used T-90 rather than T-68
   pden1000=sw_pden(s,t,pmc,0);
   ptmp1000=sw_ptmp(s,t,pmc,0);
 
   %--------------------------------------- 
   i='11067'; % top nortek  
   infile = [procpath,moor,'/nor/',moor,'_',i,'.edt'];    

       [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd');
       jd2=julian(yy,mm,dd,hh);

       bad_data=find(t==-9999); t(bad_data)=NaN;
       bad_data=find(p==-9999); p(bad_data)=NaN;
       bad_data=find(u==-9999); u(bad_data)=NaN;
       bad_data=find(v==-9999); v(bad_data)=NaN;
       bad_data=find(w==-9999); w(bad_data)=NaN;
       bad_data=find(hdg==-9999); hdg(bad_data)=NaN;
       bad_data=find(pit==-9999); pit(bad_data)=NaN;
       bad_data=find(rol==-9999); rol(bad_data)=NaN;
       bad_data=find(uss==-9999); uss(bad_data)=NaN;
       bad_data=find(vss==-9999); vss(bad_data)=NaN;
       bad_data=find(wss==-9999); wss(bad_data)=NaN;
       bad_data=find(ipow==-9999); ipow(bad_data)=NaN;
       bad_data=find(cs==-9999); cs(bad_data)=NaN;
       bad_data=find(cd==-9999); cd(bad_data)=NaN;   
   
       
  jdxticks =julian(intervalstr(1,:)):(julian(intervalstr(2,:))-julian(intervalstr(1,:)))/5:julian(intervalstr(2,:));
  gxticks = gregorian(jdxticks);
  xticks = gxticks(:,1:4);
  xticklabels = datestr(gxticks,'dd/mm HH:00');     
  
%set figure size on screen for better viewing
layout='portrait';
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

%set print area of figure
    pos1  = [1/8*scnsize(3),8*bdwidth,1/2*scnsize(3),(scnsize(4) - 30*bdwidth)];
   fig_plot=figure('Position',pos1);
    set(fig_plot,'PaperUnits','centimeters');
    set(fig_plot, 'PaperType', 'A4');
    set(fig_plot, 'PaperOrientation',layout);
    papersize = get(fig_plot,'PaperSize');
    width=17; height=26; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
    figuresize = [left, bottom, width, height];
    set(fig_plot, 'PaperPosition', figuresize);

 
  ax(1)=subplot(5,1,1);
  plot(jd,pmc)
  ylabel('Pres. (dbar)')
  set(gca,'XLim',[min(jdxticks) max(jdxticks)])
  set(gca,'XTick',jdxticks,'XTickLabel',xticklabels)  
  
  ax(2)=subplot(5,1,2);
  plot(jd,ptmp1000)
  ylabel('Ptemp (degC)')
  set(gca,'XLim',[min(jdxticks) max(jdxticks)])
  set(gca,'XTick',jdxticks,'XTickLabel',xticklabels)
  
  ax(3)=subplot(5,1,3);
  plot(jd,s)
  ylabel('Salinity')
  set(gca,'XLim',[min(jdxticks) max(jdxticks)])
  set(gca,'XTick',jdxticks,'XTickLabel',xticklabels)
  
  ax(4)=subplot(5,1,4);
  plot(jd,pden1000)  
  ylabel('Pden (kg/m^3)')
  set(gca,'XLim',[min(jdxticks) max(jdxticks)])
  set(gca,'XTick',jdxticks,'XTickLabel',xticklabels)
  
  ax(5)=subplot(5,1,5);
  plot(jd2,w)
  ylabel('W (cm/s)')
  set(gca,'XLim',[min(jdxticks) max(jdxticks)])
  set(gca,'XTick',jdxticks,'XTickLabel',xticklabels)

  linkaxes(ax,'x')
  
  print('-dpng',[moor '_mc_nor' i ])
  savefig([moor '_mc_nor' i ])
       
end