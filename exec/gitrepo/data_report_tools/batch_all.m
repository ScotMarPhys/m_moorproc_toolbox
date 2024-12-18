% script to run data report scripts and output 
global MOORPROC_G %in case cleared

moorlist = {'wb1_17_2023','wb2_18_2023'};%,'wb1l1_1_2020','wb2l15_15_2020','wbal9_9_2020','wbadcp_16_2023'};

%  If proclvl = '2': processing is using the .use files; proclvl = '3' : processing is using the .microcat and .edt files
proclvl = '2';

% Stats table on .use files
for ijk= 1:length(moorlist)
 	stats_table(moorlist{ijk}) %paths set in moor_inoutpaths
end

%----------------------------------------------------------------------------------------------------		
% Plots	
for ijk= 1:length(moorlist)
    pd = moor_inoutpaths('reports');
    outpathfigs = fullfile(pd.figsdir,moor);
    if ~exist(outpathfigs,'dir')
        mkdir(outpathfigs)
    end
    batch_data_report(moor,'procpath',procpath,'outpathfigs',outpathfigs,'proclvl',proclvl,'unfiltered')

    %  %   ADCP_allbin_plot(moor,'procpath',procpath,'proclvl',proclvl) % SJ
% Script to look into correlation between pressure of the instrument and
% horizontal current on IB3:
%vel_vs_sal_pres_corr(moor,'procpath',procpath,'proclvl','3','unfiltered')
end

% else
% %-------------------------------------------------------
% % Specific period (e.g. SCVs)
% % intervalstr = '[2014 08 05 00; 2014 09 05 00]'; % SCV?
% %intervalstr = '[2018 01 01 00; 2018 05 01 00]'
% intervalstr = '[2018 02 05 00; 2018 02 15 00]'
% for ijk= 1:length(moorlist)
% 	timelimit=str2num(intervalstr);
% 	strtime = [datestr(datenum(timelimit(1,1:3)),'yyyymmdd') '_' datestr(datenum(timelimit(2,1:3)),'yyyymmdd')];
% 	eval(['cd ' outpathfigszoom]);	  
% 	close all
% 	moor=[moorlist{ijk} '_' depyear];
% 	if exist([outpathfigszoom strtime filesep moor])~=7
% 	   eval(['!/bin/mkdir -p ' strtime filesep moor]);
% 	   cd([strtime filesep moor]);
% 	else
% 	   cd([strtime filesep moor]);
% 	end
%     batch_data_report(moor,'procpath',procpath,'proclvl','3','plot_interval',intervalstr,'unfiltered')
% end
% end %fi 0
