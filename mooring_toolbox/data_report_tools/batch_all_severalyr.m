clear all
close all

currentdir = pwd;

basedir  = '/media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap/data/moor/'; %'/home/sa02lh/Data/Dropbox/testmoorOSNAP/';%
procpath = [basedir 'proc/'];
outpathstats  = '/media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap/Documents/datareports/stats/'; %'/home/sa02lh/Data/Dropbox/testmoorOSNAP/proc/datareports/stats/';%		
outpathfigs = '/media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap/Documents/datareports/figs/';  %'/home/sa02lh/Data/Dropbox/testmoorOSNAP/proc/datareports/figs/';%'

outpathfigszoom = '/media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap/Documents/datareports/figs/'; %'/home/sa02lh/Data/Dropbox/testmoorOSNAP/proc/datareports/figszoom/';%

proclvlselec = '3'
depyear ={'01_2014';'02_2015'};
%moorlist = {'rteb1','rtwb1','rtwb2'};
moorlist = {'nocm1','nocm2','nocm3','nocm4','nocm5'};
% 
% for ijk= 1:length(moorlist)
% 	stats_table([moorlist{ijk} '_' depyear],'procpath',procpath,'outpath',outpathstats)					
% end		
% 												

for ijk= 1:length(moorlist)
    moorselect  ={};
    depyrselect =[];
    eval(['cd ' outpathfigs])	  
    close all
    for ill=1:length(depyear)
        moorselect={ moorselect{:}; [moorlist{ijk} '_' depyear{ill}]};
        depyrselect = [depyrselect '_' depyear{ill}];
    end
    moorrep=[moorlist{ijk} '_' depyrselect];
	if exist([outpathfigs moorrep])~=7
        eval(['mkdir ' moorrep]);
        cd(moorrep);
    else
        cd(moorrep);
    end

    batch_data_report_severalyr(moorselect,'procpath',procpath,'proclvl',proclvlselec);%,'unfiltered')
 eval(['cd ' currentdir])	   
end

