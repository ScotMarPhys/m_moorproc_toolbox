
%======================== SET DIRECTORIES FOR MERGE =======================
basedir                 = [pathosnap filesep];
hydrodir                = [basedir 'data/moor/proc/hydro_grid_LD_TEST/'];
grdatdir                = [basedir 'data/moor/proc/hydro_grid_merged/'];
boundarydir             = [basedir 'exec/ar30/stage3/gridding/MCAT/']; 

%======================== MERGE INITIALISATION =============================
gridding                = 1  ;  % 1: linear, 2: using climatological profiles
bathy                   = false ;  % turns on/off the bathy charts. off = flase
mc_check_plot           = true; %false ;  % turns on/off the microcat check plots. off =false

jg_start                = datenum(2014,6,01,00,00,00);
jg_end                  = datenum(2018,7,10,00,00,00);
lastyeardata            = '_2017'; % for datafilename

JG                      = jg_start: 0.5: jg_end; % full time series using 2 samples per day
pgg                     = 0:20:2000; % depths in 20dbar bins
depthminforhoriz_interp = 40; % in case no data are available at a specific time 
% (e.g.: mooring turn around, knock-down of the mooring head) don't interpolate on a time basis for level above 40 m 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================== MOORING FILE NAMES =============================
% ------------------ UPDATE NEW MOORING FILES HERE ----------------------- 
% path of the mooring data is defined in the startup file under osnap/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
moorings = {'rteb1_01_2014';...
            'rtwb1_01_2014';...
            'rtwb2_01_2014';...
            'rtwb1_02_2015';...
            'rtwb2_02_2015';...
            'rteb1_02_2015';...
            'rtwb1_03_2016';...
            'rtwb2_03_2016';...
            'rteb1_03_2016';...
            'rtwb1_04_2017';...
            'rtwb2_04_2017';...
            'rteb1_04_2017'};

microcat_order=(   {[1:4 6:8],...
                    [1:8],...
                    [1:2],... 
                    [1:8],...
                    [1:2],... 
                    [1:8],...
                    [1:8],...
                    [1:2],...
                    [1 3:8],...
                    [1:9],...
                    [1:2],...
                    [2:6 9:11 13]});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[indx,tf] = listdlg('PromptString',{'Select moorings for merge.',...
    'Make sure only western or eastern arrays are selected',''},...
    'ListString',moorings);
if tf==1
    moor=moorings(indx);       
else
    disp('No moorings selected');
end

% change file naming convention for .dat files by inclduing 'osnap' in
% file name and removing number after rteb
for ii=1:numel(moor)
   moor{ii}=[moor{ii}(1:4) '_osnap' moor{ii}(6:end)];
end

moorstr=strjoin(moor,', ');
answer = questdlg(['Merge moorings : ' moorstr '?'], ...
    'Merge',...
    'Continue','Cancel','Wrong moorings','Wrong moorings');
% Handle response
switch answer
    case 'Continue'
        disp([answer ' Merging.'])

        % if data is from western or eastern array, merge accordingly,
        % the difference is the way the data is indexed in the code,
        % and could probably do with an overhaul
        if ismember(moor{1}(3:4),'eb')
            run merge_osnap_data_east.m
        else
            run merge_osnap_data_west.m
        end

    case 'Cancel'
    case 'Wrong moorings'
        % re-select moorings
        [indx,tf] = listdlg('PromptString',{'Select moorings for merge.',...
            'Make sure only western or eastern arrays are selected',''},...
            'ListString',moorings);
        if tf==1
            moor=moorings(indx);       
        else
            disp('No moorings selected');
        end

        % change file naming convention for .dat files by inclduing 'osnap' in
        % file name
        for ii=1:numel(moor)
           moor{ii}=[moor{ii}(1:5) '_osnap' moor{ii}(6:end)];
        end

        disp([answer ' Merging.'])

        if ismember(moor{1}(3:4),'eb')
            run merge_osnap_data_east.m
        else
            run merge_osnap_data_west.m
        end
end