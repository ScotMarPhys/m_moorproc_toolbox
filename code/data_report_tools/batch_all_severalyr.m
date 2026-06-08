% clear all
close all
clearvars  -except pathosnap 

currentdir = pwd;

basedir  = [pathosnap '/data/moor/']; %'/home/sa02lh/Data/Dropbox/testmoorOSNAP/';%
procpath = [basedir 'proc/'];
reportdir = [pathosnap '/Documents/datareports/'];%'/home/sa02lh/Data/Dropbox/testmoorOSNAP/proc/datareports/'; % /media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap/Documents/datareports/
outpathstats  = [reportdir 'stats/'];		
outpathfigs = [reportdir 'figs/'];
outpathfigszoom = [pathosnap '/Documents/datareports/figs/']; %'/home/sa02lh/Data/Dropbox/testmoorOSNAP/proc/datareports/figszoom/';%

proclvlselec = '2'; % '2': processing is using the .use files; '3' : processing is using the .microcat and .edt files
%depyear ={'01_2014';'02_2015';'03_2016'};
%depyear ={'02_2015';'03_2016';'04_2017'};
depyear ={'04_2017';'05_2018'};
moorlist = {'rteb1'}; %,'rtwb1','rtwb2'
%moorlist = {'nocm1'};
% moorlist = {'nocm1','nocm2','nocm3','nocm4','nocm5'};

										

for ijk= 1:length(moorlist)
    moorselect  ={};
    depyrselect =[];
    eval(['cd ' outpathfigs])	  
    close all
    for ill=1:length(depyear)
        moorselect{ill}= [moorlist{ijk} '_' depyear{ill}];
        depyrselect = [depyrselect '_' depyear{ill}];a
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

