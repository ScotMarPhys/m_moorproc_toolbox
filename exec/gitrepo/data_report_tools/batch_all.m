
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
%outpathfigszoom = '/Users/locupe/Dropbox/Work/Dataproc/OSNAP_mooring/moorOSNAP_quick_analyses/proc_plot/datareports/figszoom/'; 
outpathfigszoom = [reportdir 'figs/'];
%-------------------------------------
% Selection of the deployment year
% depyear ='01_2014';
% depyear ='02_2015';
%depyear ='03_2016';
%depyear ='04_2017';

%depyear ='05_2018';

depyear ='01_2018';

%-------------------------------------
% Selection of the mooring to process
% moorlist ={'nocm1'};
% moorlist ={'nocm1','nocm2','nocm3','nocm4','nocm5'};
% moorlist = {'rteb1','rtwb1','rtwb2'};
% moorlist = {'rtwb1','rtwb2'};
%moorlist = {'rtadcp1'};%
%moorlist = {'ib4','ib3','ib5'};
moorlist = {'ib5'};
 
proclvl = '2';%  If proclvl = '2': processing is using the .use files; proclvl = '3' : processing is using the .microcat and .edt files
%  
% %----------------------------------------------------------------------------------------------------
% Stats table on .use files
for ijk= 1:length(moorlist)
 	stats_table([moorlist{ijk} '_' depyear],'procpath',procpath,'outpath',outpathstats)					
end

%----------------------------------------------------------------------------------------------------		
% PLots												
if 1
for ijk= 1:length(moorlist)
  eval(['cd ' outpathfigs])	  
    close all
    moor=[moorlist{ijk} '_' depyear];
	if exist([outpathfigs moor])~=7
        eval(['mkdir ' moor]);
        cd(moor);
    else
        cd(moor);
    end
    
     batch_data_report(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
%  %   ADCP_allbin_plot(moor,'procpath',procpath,'proclvl',proclvl) % SJ


% Script to look into correlation between pressure of the instrument and
% horizontal current on IB3:
%vel_vs_sal_pres_corr(moor,'procpath',procpath,'proclvl','3','unfiltered')

eval(['cd ' currentdir])	   
end

else
%-------------------------------------------------------
% Specific period (e.g. SCVs)
% intervalstr = '[2014 08 05 00; 2014 09 05 00]'; % SCV?
%intervalstr = '[2018 01 01 00; 2018 05 01 00]'
intervalstr = '[2018 02 05 00; 2018 02 15 00]'
for ijk= 1:length(moorlist)
	timelimit=str2num(intervalstr);
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
    batch_data_report(moor,'procpath',procpath,'proclvl','3','plot_interval',intervalstr,'unfiltered')
end
end %fi 0
