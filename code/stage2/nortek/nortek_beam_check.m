% function nortek_beam_check(moor)
% function nortek_beam_check(moor, pd)
%
% Function to load .dia files from Norteks to generate figures for assessing
% for beam problems. Up till now this hasn't been part of normal moored
% data processing.
%
% Function can be called at any point once .dia files have been transferred
% to data directory on workstation, but if running after .dat files have
% been converted to rodb format and processed by stage 2 to .use files then
% can easily plot the actual timeseries data with the noise floor as 
% detected by this script. So storing this script in the stage 2 folder.
%
% Diagnostic filenames will all have .dia extension but prefix varies depending
% on setup. Prefix should match that for the .dat files processed by
% nortek2rodb_01.m
% If Nortek filenames do not have one of a few conventional formats
% (serialnum.dia, serialnum_data.dia, serialnum01.dia), they can be read
% from a separate text file with serial number in first column and filename
% in second, either space or tab delimited.  
% Naming of text file is as eb3_1_200407_filenames.txt, and it (if
% required) should be in same directory as raw files. 
%
% required inputs: moor - mooring name e.g. 'wb2_3_200606'
%
% optional inputs: pd - structure containing paths to raw data files,
%   info.dat file, and output files.  Otherwise these will be set by
%   calling moor_inoutpaths.m
%
% functions called: moor_inoutpaths.m (if single input argument)
%                   hms2h.m
%                   rodbload.m
%
% written by Darren Rayner in July 2026 on DY214


function nortek_beam_check(moor, varargin)

global MOORPROC_G

if nargin==0
    help nortek_beam_check
    return
end

if nargin==1
    pd = moor_inoutpaths('nor',moor);
else
    pd = varargin{1};
end

operator = MOORPROC_G.operator;

% ----- read infofile / open logfile  ------------------------------------

infovar = 'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth'; 
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd]  =  rodbload(pd.infofile,infovar);

bg = datenum([s_d(:)' s_t(:)' 0]); %start
ed = datenum([e_d(:)' e_t(:)' 0]); %end

mnor = (id==368 | id==370); % Nortek instrument code
serial_nums=sn(mnor)
z = z(mnor);
id = id(mnor);

%first try to find files with names in a few standard forms (adapting code
%from nortek2rodb_01.m
filenames = cell(length(serial_nums),1);
pats = {'%d_data.dia';'%d.dia';'%d01.dia'};
for sno = 1:length(serial_nums)
    n = 1;
    while isempty(filenames{sno}) && n<=length(pats)
        fn = sprintf(pats{n},serial_nums(sno));
        if exist(fullfile(pd.rawpath,fn),'file')
            filenames{sno} = fn;
        else
            n = n+1;
        end
    end
end
%now, if necessary, try looking in text file
fnmiss = find(cellfun('isempty',filenames));
if ~isempty(fnmiss)
    % load text file containing serial numbers with related filenames. NB: No
    % header in this file, but column 1 should be serial number and column 2
    % should be filename. This is necessary as there is currently no standard
    % naming convention for Aquadopp files. Do not include path though.
    try
       %load and add to list
       [file_ser_nums, filenames_list]=textread(pd.listfile,'%d %s');
       [~,iinfo,ilist] = intersect(serial_nums(fnmiss),file_ser_nums,'stable');
       filenames(fnmiss(iinfo)) = filenames_list(ilist);
       list_not_used = setdiff(filenames_list,filenames);
       if ~isempty(list_not_used)
           warning('some files in list not used:')
           fprintf(1,'%d\n',list_not_used(:))
       end
    catch
    end
end
%check and warn about any still not found
fnmiss = find(cellfun('isempty',filenames));
if ~isempty(fnmiss)
    warning('these S/Ns from info file do not have standard names for .dat/.dia files and were not found in listfile:')
    fprintf(1,'%d\n',serial_nums(fnmiss))
    warning('Create list of S/Ns and data file names for these in %s and try again',pd.listfile)
    c = input('or ''skip'' to continue processing skipping these S/Ns (e.g. if data not downloaded)  ','s');
    if strcmp(c,'skip')
        serial_nums(fnmiss) = [];
        z(fnmiss) = [];
        filenames(fnmiss) = [];
    else
        error('missing Nortek files')
    end
end

% -------- load data --------------
for proc = 1:length(filenames)
    infile=fullfile(pd.rawpath,filenames{proc});
    if exist(infile,'file')
        fprintf(1,'Processing diagnostic file for %d\n', serial_nums(proc))
    else
       checkans = input(['Diagnostic file for serial number ' serial_nums(proc) 'not found. Do you want to continue to next? y/n [y]:'],'s');
       if strcmpi(checkans,'n')
           error('stopping on S/N %d file not found',serial_nums(proc))
       else
           %on to next
           continue
       end
    end
    
    all_data=load(infile);
    month=all_data(:,1); day=all_data(:,2); year=all_data(:,3); 
    hour=all_data(:,4); minute=all_data(:,5); second=all_data(:,6);
    Amp{1}=all_data(:,12); Amp{2}=all_data(:,13); Amp{3}=all_data(:,14);
        
    dnum      = datenum([year month day hour minute second]);
    valI      = find(dnum<ed & dnum>bg);
    dnum=dnum(valI);
    Amp{1}=Amp{1}(valI);
    Amp{2}=Amp{2}(valI);
    Amp{3}=Amp{3}(valI);

    dnum_diff=diff(dnum);

    % detection of first sample of diagnsotic burst assumes diagnostics are
    % always at 1Hz - this will likely need updating for Mk2 Aquadopp when
    % get round to including them
    FirstSampleIDX=[0;find(round(dnum_diff*60*60*24,1)~=1)]+1; % use round as may not be exactly 1s difference in the wonderful land of Matlab

    % setup figure for data from diagnostic file
    figure
    max_y=0; % set values to be replaced in first run through loop below
    min_y=100;
    NoiseFloor=cell(1,3); % pre-define size
    for sub=1:3
        subplot(1,3,sub)
        hold on
        plot(dnum,Amp{sub}','b.');
        plot(dnum(FirstSampleIDX),Amp{sub}(FirstSampleIDX),'r.');
        ylabel('Counts')
        legend('All samples','First sample')
        title(['Beam ' num2str(sub) ' diagnostic signal strength'])

        % create variable for noise floor from 1st sample so can use with
        % non-diagnostic data later
        NoiseFloor{sub}=Amp{sub}(FirstSampleIDX);
        % ensure all 3 panels have same ylimits - find here and apply in next loop
        ylims=get(gca,'ylim');
        if ylims(2)>max_y
            max_y=ylims(2);
        end
        if ylims(1)<min_y
            min_y=ylims(1);
        end
    end
    for sub=1:3
        subplot(1,3,sub)
        set(gca,'ylim',[min_y max_y])
        datetick('x','mmm-yy','keeplimits')
    end
    sgtitle(['Beam checking from diagnostic signal strength from file: ' regexprep(filenames{proc},'_','\\_')])

    % check if output folder exists for saving figure to and if not create it
    if ~exist(pd.stage2figpath,'dir')
        mkdir(pd.stage2figpath)
    end
    figname = fullfile(pd.stage2figpath,sprintf(pd.stage2form,serial_nums(proc)));
    % trim .use from file name - probably supposed to update this in
    % wherever defines pd structure but will get to that another time
    figname=figname(1:end-4);

    orient landscape
    print(gcf,'-dpng','-r300',[figname,'_diagnostic_beam_strength.png']);
    savefig([figname '_diagnostic_beam_strength.fig']);
    % ----- not going to save data to rodb for diagnostic data as further processing only really using .dat files -----------------

    % Now load actual moored timeseries data (from .dat files) to plot
    % signal strength from them with overlayed noise floor found above
    % but first check if stage2 script has been run for the mooring in
    % question - if not then put warning message but continue loop for
    % processing other instruments on the mooring
    stage2file=fullfile(pd.stage2path,sprintf(pd.stage2form,serial_nums(proc)));
    if ~exist(stage2file,'file')
        disp(['Warning: stage 2 processing not been completed for ' num2str(serial_nums(proc)) '. Unable to show normal (non-diagnostic) signal strength plot'])
        continue
    else
        % load data from stage 2 processing
        [YY,MM,DD,HH,stage2Amp{1},stage2Amp{2},stage2Amp{3}]=rodbload(stage2file,'YY:MM:DD:HH:USS:VSS:WSS');
        stage2_dnum=datenum(YY,MM,DD,HH,0,0);
        
        % setup figure for data from .dat file
        figure
        max_y=0; % set values to be replaced in first run through loop below
        min_y=100;
        for sub=1:3
            subplot(1,3,sub)
            hold on
            plot(stage2_dnum,stage2Amp{sub}','b.');
            plot(dnum(FirstSampleIDX),Amp{sub}(FirstSampleIDX),'r.');
            ylabel('Counts')
            legend('Main timeseries data','Noise floor from diagnostics')
            title(['Beam ' num2str(sub) ' signal strength'])
    
            % ensure all 3 panels have same ylimits - find here and apply in next loop
            ylims=get(gca,'ylim');
            if ylims(2)>max_y
                max_y=ylims(2);
            end
            if ylims(1)<min_y
                min_y=ylims(1);
            end
        end
        for sub=1:3
            subplot(1,3,sub)
            set(gca,'ylim',[min_y max_y])
            datetick('x','mmm-yy','keeplimits')
        end
        sgtitle(['Beam checking from main data signal strength from file: ' regexprep([moor '_' num2str(serial_nums(proc)) '.use'],'_','\\_')])
        figname = fullfile(pd.stage2figpath,sprintf(pd.stage2form,serial_nums(proc)));

        orient landscape
        print(gcf,'-dpng','-r300',[figname,'_main_data_beam_strength.png']);
        savefig([figname '_main_data_beam_strength.fig']);
    end

end
disp(' ')
disp('Bad beams can manifest as drops in signal strength to near the noise floor level.')
disp('Each beam should be fairly similar for signal strength, so drops in one but not the others could also indicate beam failure.')
disp('Noise floor is found from the first sample of a diagnostic burst (saved in .dia raw file).')
disp('Timeseries may be questionable if signal levels are near noise floor throughout,')
disp('but this isn''t necessarily beam failure and could be due to low backscatters or too large a blanking distance.')
