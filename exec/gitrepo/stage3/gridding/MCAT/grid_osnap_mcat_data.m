function grid_osnap_mcat_data(varargin)
% grid_osnap_mcat_data()
% grid_osnap_mcat_data(sec_seg)
% grid_osnap_mcat_data(sec_seg, parameter, value)
%
% sec_seg (optional) is 'rt' or 'ib'

%% GRID AND MERGE DATA for OSNAP mooring
% _*Original author: Lewis Drysdale, 2020*_

global MOORPROC_G
close('all')

if nargin==0
    sec_seg = input('which set of moorings (e.g. rt or ib)?  ','s');
else
    sec_seg = varargin{1}; varargin(1) = [];
end

%defaults and input arguments
p_hydrogrid.uncal = 0; %use calibrated data

% 1.SET GRID AND QA PARAMETERS
p_hydrogrid.dum         = -9999.0000;
p_hydrogrid.c1535       = 42.914;
p_hydrogrid.mcat        = [332:337];    % instrument model numbers
p_hydrogrid.int_step    = 10;           % vertical interpolation step
p_hydrogrid.preverse    = 4000;         % 4000 pressure level below whch deep temperature reversion may occur

% 2.REPAIR AND DESPIKE SETTINGS
p_hydrogrid.gap_max      = 10;          % allow for a maximum of gap [days] in data to interpolate across
p_hydrogrid.y_tol        = [-10 10];    % deviation in PSU allowed by depsike routine
p_hydrogrid.stddy_tol    = 4;           % tolerance range of differences between adjacent values of y
p_hydrogrid.nloop        = 5;           % despike loop number
p_hydrogrid.graphics     = 'n';

% 3.GRID INITIALISATION
p_hydrogrid.p_gridsize   = 20;
p_hydrogrid.max_depth               = 1780; %set max depth of valid data
p_hydrogrid.pgg = 0:p_hydrogrid.p_gridsize:p_hydrogrid.max_depth;
p_hydrogrid.v_interp     = 'linear';
p_hydrogrid.co           = 1/2; % filter cut off frequency [1/days]
p_hydrogrid.iss          = 12; % initial sub-sampling frequency [1/days]
p_hydrogrid.fss          = 2;  % final sub-sampling frequency [1/days]

% 4.SET DIRECTORIES FOR MERGE
pd = moor_inoutpaths('mcgrid',[]);
fn = fieldnames(pd);
for no = 1:length(fn)
    if ~exist(pd.(fn{no}),'dir')
        mkdir(pd.(fn{no}))
    end
    p_hydrogrid.(fn{no}) = pd.(fn{no});
end

% 5.MERGE INITIALISATION
% gridding                = 1  ;  % 1: linear, 2: using climatological profiles
% bathy                   = false ;  % turns on/off the bathy charts. off = flase
% mc_check_plot           = true; %false ;  % turns on/off the microcat check plots. off =false

jg_start                = datenum(2014,07,01,00,00,00);
jg_end                  = datenum(MOORPROC_G.YEAR,07,31,00,00,00); %***
% jg_start_str = datestr(jg_start,'YYYYmm');
% jg_end_str = datestr(jg_end,'YYYYmm');
% lastyeardata            = num2str(MOORPROC_G.YEAR-2);
% 
% JG                      = jg_start: 0.5: jg_end; % full time series using 2 samples per day
% depthminforhoriz_interp = 0; % in case no data are available at a specific time, % multiples of 20
% % (e.g.: mooring turn around, knock-down of the mooring head) don't interpolate on a time basis for level above 40 m
% 
% data_version = 'v0';

% MOORING FILE NAMES
% UPDATE NEW MOORING FILES HERE
% path of the mooring data is defined in the startup file under osnap/

% overwrite defaults with other inputs
for no = 1:2:length(varargin)
    eval([varargin{no} ' = varargin{no+1};']);
end

if ~exist('moorings','var')
switch sec_seg
    case 'rt'
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
            'rteb1_04_2017';...
            'rtwb1_05_2018';...
            'rtwb2_05_2018';...
            'rteb1_05_2018';...
            'rtwb1_06_2020';...
            'rtwb2_06_2020';...
            'rteb1_06_2020'};
    case 'ib'
        moorings = {'ib3_01_2018';...
            'ib4_01_2018';...
            'ib5_01_2018';...
            'ib3_02_2020';...
            'ib4_02_2020';...
            'ib5_02_2020';...
            'ib3_03_2022';...
            'ib4_03_2022';...
            'ib5_03_2022'};
end
end

if p_hydrogrid.uncal
% order of microCATs (serial numbers 335:337) in info.dat file, but exclude double deployed ODO
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
    [2:6 9:11 13],...
    [1:9],...
    [1:3],...
    [1:5 8:12] ...
    [1:9] ...
    [1:3] ...
    [1:5 8:12]});
instrument_order=(  {[1:4 6:8],...
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
    [2:6 9:11 13],...
    [9140 14355 3481 11329 10559 14368 11337 14353 14354],...
    [9377 3218 9372],...
    [11290 11289 11288 11287 10577 3276 10560 10562 9113 9375],...
    [11343,11342,14364,7923,7924,11137,11139,10576,11465],...
    [10575,11341,13020],...
    [9141,11322,9396,21560,11327,11330,11334,11335,9390,11338]});
end


% CHECK FOR FOLDERS

for ii=1:numel(moorings)
    if ~exist(fullfile(pd.figdir,moorings{ii}),'dir')
        mkdir(fullfile(pd.figdir,moorings{ii}))
    end
end

%% GRID AND MERGE

if ~exist('gridans','var')
gridans = questdlg(['Grid or merge moorings?'],...
    'Grid and merge',...
    'Grid','Merge','Merge');
end
% Handle response
switch gridans

    case 'Grid'

        if exist('indx','var')
            tf = 1;
        else
        [indx,tf] = listdlg('PromptString',{'Select the moorings for interpolation.',...
            '',''},...
            'ListString',moorings);
        end
        if tf==1
            mlist=moorings(indx);
            if p_hydrogrid.uncal
            ilist=microcat_order(indx);
            tlist=instrument_order(indx);
            end
        else
            disp('No mooring selected');
        end

        for ii=1:numel(indx) % for each mooring to process

            moor=mlist{ii}; % select the first moooring to process

            if p_hydrogrid.uncal
            p_hydrogrid.mc_ind=cell2mat(ilist(ii)); % get instrument order
            end

            %%
            % _*if there is a requirement to process uncalibrtaed data make next line uncommented,
            % and update serial numbers above*_
            %if p_hydrogrid.uncal
            % p_hydrogrid.mc_int=cell2mat(tlist(ii)); %get instrument list
            %else
            % p_hydrogrid.mc_int=[]; % get instrument list
            %end

            % SET PATHS

            p_hydrogrid.moor         = moor;
            p_hydrogrid.mooringpath  = pd.mooringpath;
            p_hydrogrid.outname      = [moor '_grid.mat'];

            % INTERPOLATE
            hydro_grid_osnap_linear_interp(p_hydrogrid);
        end
        %
        %  OPTIONS TO MERGE

        answer = questdlg(['Do you want to merge data?'], ...
            'Merge',...
            'Yes','No','No');
        % Handle response
        switch answer
            case 'Yes'
                [indx,tf] = listdlg('PromptString',{'Select moorings for merge.',...
                    'Make sure only western or eastern arrays are selected',''},...
                    'ListString',moorings);
                if tf==1
                    moor=moorings(indx);
                else
                    disp('No moorings selected');
                end

                % change file naming convention for .dat files by inclduing 'osnap' in
                % file name and removing deployment number
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
                        if ismember(moor{1}(3:4),'eb')
                            merge_osnap_data_east
                        else
                            merge_osnap_data_west
                        end
                    case 'Cancel'
                    case 'Wrong moorings'
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
                            moor{ii}=[moor{ii}(1:4) '_osnap' moor{ii}(6:end)];
                        end

                        disp([answer ' Merging.'])

                        if ismember(moor{1}(3:4),'eb')
                            merge_osnap_data_east
                        else
                            merge_osnap_data_west
                        end

                end
            case 'No'
                return

        end

        % if merge only then prompt for files to merge.
    case 'Merge'

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
end