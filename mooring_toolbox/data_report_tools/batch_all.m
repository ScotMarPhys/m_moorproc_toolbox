clear all
close all

currentdir = pwd;

%-------------------------------------
% Definition of the different path
basedir  = '/media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap/data/moor/';%'/home/sa02lh/Data/Dropbox/testmoorOSNAP/'; %
procpath = [basedir 'proc/'];
reportdir = '/media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap/Documents/datareports/';%'/home/sa02lh/Data/Dropbox/testmoorOSNAP/proc/datareports/'; % /media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap/Documents/datareports/
outpathstats  = [reportdir 'stats/'];		
outpathfigs = [reportdir 'figs/'];
outpathfigszoom = '/home/sa02lh/Data/Dropbox/testmoorOSNAP/proc/datareports/figszoom/'; 
%-------------------------------------
% Selection of the deployment year
depyear ='01_2014';
%depyear ='02_2015';
%-------------------------------------
% Selection of the mooring to process
moorlist ={'nocm1','nocm2','nocm3','nocm4','nocm5'};
%moorlist = {'rteb1','rtwb1','rtwb2'};
%moorlist = {'rtadcp1'};%
 
% for ijk= 1:length(moorlist)
% 	stats_table([moorlist{ijk} '_' depyear],'procpath',procpath,'outpath',outpathstats)					
% end		
% 												
if 1
for ijk= 5 %1:length(moorlist)
  eval(['cd ' outpathfigs])	  
    close all
    moor=[moorlist{ijk} '_' depyear];
	if exist([outpathfigs moor])~=7
        eval(['mkdir ' moor]);
        cd(moor);
    else
        cd(moor);
    end
    
   % batch_data_report(moor,'procpath',procpath,'proclvl','3','unfiltered')
    ADCP_allbin_plot(moor,'procpath',procpath,'proclvl','3')
    ADCP_allbin_plot(moor,'procpath',procpath,'proclvl','2')   
 eval(['cd ' currentdir])	   
end

else
%-------------------------------------------------------
% Specific period (e.g. SCVs)
intervalstr = '[2014 08 05 00; 2014 09 05 00]'; % SCV?
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
