function varargout=stats_table(moor,varargin)
% function varargout=stats_table(moor,'procpath','layout','non-verbose')
%
% function to generate a simple stats table for data reports
%
% required inputs:-
%   moor: complete mooring name as string. e.g. 'wb1_1_200420'
%
% optional inputs:-
%   layout: orientation of figure portrait/lanscape (default = portrait)
%           input of 'landscape' or 'portrait'
%           e.g. pressure_overlay('wb1_1_200420','layout','landscape')
%   procpath: specify exact procpath
%           e.g. pressure_overlay('wb1_1_200420','procpath','/Volumes/noc/mpoc/hydro/rpdmoc/rapid/data/moor/proc/')
%   outpath: specify exact output path
%   non-verbose: a mode to output the stats for other routines without
%           saving an ascii file. If called in this mode there will be an
%           output to the function consisiting of the Microcat serial
%           numbers and statistics
%
%   output: an ascii text file containing the summary statistics for the
%   mooring - naming convention is moor_stats.asc where moor is the mooring
%   name as input to the function
%
% functions called:-
%   rodbload, julian, auto_filt, gregorian
%   from .../exec/moor/tools and .../exec/moor/rodb paths
%
% Routine written by Darren Rayner February 2008.
%
% Doesn't do anything with MMP or ADCP data.
%
% Loic Houpert, 4/10/16, add output path
%
%***************************
%***************************
% MAJOR PROBLEM
%***************************
%***************************
% .use files do not have bad data changed to -9999 so statistics will be
% badly skewed by bad data such as when a microcat has a pressure sensor
% failure. Need to check with Torsten if there is a better input file to
% use.

global MOORPROC_G

if nargin <1
    help stats_table
    return
end

%defaults
layout = 'portrait';
procpath = fullfile(MOORPROC_G.moordatadir,'proc');
outpath = fullfile(MOORPROC_G.reportdir,'stats');
non_verbose = 0;
%and optional inputs overwrite them
n = 2; 
while n<=nargin
    if ischar(varargin{n}) 
        if strcmp(varargin{n},'non-verbose')
            non_verbose = 1;
            n = n+1;
        else
            eval([varargin{n} ' = varargin{n+1};'])
            n = n+2;
        end
    end
end

infofile = fullfile(procpath,moor,[moor 'info.dat']);

%-----------------------------------------------------
% Load vectors of mooring information
% id instrument id, sn serial number, z nominal depth of each instrument
% s_t, e_t, s_d, e_d start and end times and dates
% lat lon mooring position, wd corrected water depth (m)
% mr mooring name
[id,sn,z,s_t,s_d,e_t,e_d]  =  rodbload(infofile,...
    'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date');


% JULIAN Convert Gregorian date to Julian day.
% JD = JULIAN(YY,MM,DD,HH) or JD = JULIAN([YY,MM,DD,HH]) returns the Julian
% day number of the calendar date specified by year, month, day, and decimal
% hour.
% JD = JULIAN(YY,MM,DD) or JD = JULIAN([YY,MM,DD]) if decimal hour is absent,
% it is assumed to be zero.
% Although the formal definition holds that Julian days start and end at
% noon, here Julian days start and end at midnight. In this convention,
% Julian day 2440000 began at 00:00 hours, May 23, 1968.
jd_start = julian([s_d' hms2h([s_t;0]')]);
jd_end   = julian([e_d' hms2h([e_t;0]')]);

if non_verbose==0
    disp('z : instrument id : serial number')
    for iid = 1:length(id)
        disp([z(iid),id(iid),sn(iid)])
    end
end

%get information about all the instruments including prefixes and filenames/paths, in a table we can loop through
id_z_sn = all_inst_table(id, z, sn);
id_z_sn.data_loaded = false(length(id),1);

% -----------------------------------
% START OF READING IN INSTRUMENT DATA
% -----------------------------------

for iid=1:length(id_z_sn.id)
    clear data
    if ~isempty(id_z_sn.dirs{iid})
        if ~non_verbose
            disp('*************************************************************')
            disp(['Reading ' id_z_sn.inst{iid} ' - ',num2str(id_z_sn.sn(iid))])
            disp('*************************************************************')
        end
        infile = fullfile(procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%0.4d%s.use',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
        iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));

        %read data into structure array
        fileopen=fopen(infile,'r');
        if fileopen>0
            varstr = ['yy:mm:dd:hh:' id_z_sn.vars{iid}];
            clear a; a(1).a = [];
            [yy,mm,dd,hh,a(1).a,a(2).a,a(3).a,a(4).a,a(5).a,a(6).a,a(7).a,a(8).a,a(9).a] = rodbload(infile,varstr);
            vars = split(id_z_sn.vars{iid},':');
            for no = 1:length(vars)
                data.(vars{no}) = a(no).a;
            end
            data.jd=julian(yy,mm,dd,hh);
        else
            disp('File does not exist!')
            disp(['infile = ' infile])
            data = [];
        end

        if ~isempty(data)
            %calculate stats
            %actually calculate the stats, which are the same for all fields other than
            %jd no matter the instrument***
            fn = fieldnames(data);
            fn = setdiff(fn,'jd');
            for no = 1:length(fn)
                vnam = fn{no};
                data.(vnam)(data.(vnam)==-9999) = NaN;
                data.([vnam 'mean']) = nanmean(data.(vnam));
                data.([vnam 'std']) = nanstd(data.(vnam));
                data.([vnam 'max']) = max(data.(vnam));
                data.([vnam 'min']) = min(data.(vnam));
            end
            data.samples = length(data.jd);
            data.start_date = gregorian(data.jd(1));
            data.end_date = gregorian(data.jd(end));

            if isfield(data,'u')
                % calculate speed and direction
                data.spd = sqrt(data.u.^2 + data.v.^2);
                data.dir=atan2(data.v, data.u)*180/pi;
                % calculate spd and dir means and stds
                data.dirmean=atan2(data.vmean,data.umean)*180/pi;
                data.spdmean=sqrt(data.umean.^2 + data.vmean.^2);
                data.spdstd = nanstd(data.spd);
                data.spdmin = min(data.spd);
                data.spdmax = max(data.spd);
                % for dir STD convert directions to values around mean direction
                % i.e. mean direction becomes 0, and all other directions are relative to
                % that.
                dir_new = data.dir-data.dirmean;
                dir_new(dir_new>180) = dir_new(dir_new>180)-360;
                dir_new(dir_new<-180) = dir_new(dir_new<-180)+360;
                data.dirstd = nanstd(dir_new);
                data.dirmin = min(dir_new)+data.dirmean;
                data.dirmax = max(dir_new)+data.dirmean;
            end
            alldata.(iname) = data;

            id_z_sn.data_loaded(iid) = true;
        end
    end
end

%only keep the rows where we loaded data
id_z_sn = id_z_sn(id_z_sn.data_loaded,:);

% ============================================
% START OF OUPUTTING SUMMARY DATA TO TEXT FILE
% ============================================
if ~exist(outpath,'dir')
    mkdir(outpath)
end
outfile=fullfile(outpath,[moor '_stats.asc']);
check=0;
while check==0
    if exist(outfile,'file')
        disp(['Outfile ' outfile ' already exists!'])
        overwrite=input('Do you wish to overwrite the file? y/n (default=n):- ','s');
        if overwrite~='y'
            outfile=input('Please enter alternative outfile name:- ','s');
        else
            check=1;
        end
    else
        check=1;
    end
end
fid=fopen(outfile,'w');

fprintf(fid,'%s Moorig Array. \n%s %s \n',MOORPROC_G.project,'Simple Statisctics for Mooring:- ',moor);
fprintf(fid,'%s %0.2i/%0.2i/%i %0.2i:%0.2i\n','Mooring deployment - start: ',s_d(3),s_d(2),s_d(1),s_t(1),s_t(2));
fprintf(fid,'%s %0.2i/%0.2i/%i %0.2i:%0.2i\n\n','                       end: ',e_d(3),e_d(2),e_d(1),e_t(1),e_t(2));
fprintf(fid,'%s\n','-----------------------------------------------------------------------------------------');
fprintf(fid,'%s\n','       SN  var       first            last        valid    mean   stdev     min     max');
fprintf(fid,'%s\n','                    record          record      records');
fprintf(fid,'%s\n','-----------------------------------------------------------------------------------------');

% remove double occurence of Sonteks with SMPs
% and remove MMPs as this routine is too simple for them
% also remove releases as some info.dat files have them in

% NOTE this routine assumes that all current meters have pressure and
% temperature sensors, and that the RCM11 and S4 also have conductivity
% this is likely to change in future as new CMs are unlikely to have cond
% and possibly not temp.

for iid=1:size(id_z_sn,1)
    variables=[' p ';' t ';' c ';'o2 ';'ph ';'phv';' u ';' v ';'spd';'dir']; % pad single character variables with spaces
    if sum(strcmp({'MC' 'RBR' 'IDR'},id_z_sn.inst{iid}))
        stats_rows=3;
        variables=variables(1:3,:);
    elseif sum(strcmp({'ODOMC'},id_z_sn.inst{iid}))
        stats_rows=4;
        variables=variables(1:4,:);
    elseif sum(strcmp({'RCM11' 'S4'},id_z_sn.inst{iid}))
        stats_rows=7;
        variables=variables([1:3 7:10],:); %no o or ph
    elseif sum(strcmp({'ARG' 'NOR'},id_z_sn.inst{iid})) % i.e. is a cm with t and p but not c (NOR or ARG)
        stats_rows=6;
        variables=variables([1:2 7:10],:); % remove c, o, and ph variables
    elseif sum(strcmp({'BPR'},id_z_sn.inst{iid})) % i.e. is a SBE BPR
        stats_rows=2;
        variables=variables(1:2,:);
    elseif sum(strcmp({'DVS'},id_z_sn.inst{iid})) % current meter with t but not c or p
        stats_rows=5;
        variables=variables([2 7:10],:); % only t and currents
    end

    SN = repmat({' '},1,stats_rows);
    SN{1}=num2str(id_z_sn.sn(iid));
    iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
    data = alldata.(iname);

    % Use automatic detection of first and last useable record
    % NB: this works on the assumption that bad data has been replaced with
    % -9999 values - which are then replaced with NaNs during loading
    j2=size(variables,1);
    for j=1:j2
        norexist = isfield(alldata,iname);
        if norexist ==0
            validRecs = [];
        else
            vname = strtrim(variables(j,:));
            validRecs = find(~isnan(data.(vname)));
        end
        if ~isempty(validRecs)
            firstRec(j)=validRecs(1);
            lastRec(j)=validRecs(end);
            numRecs(j)=length(validRecs);
            firstRecJD = data.jd(firstRec(j));
            lastRecJD = data.jd(lastRec(j));
            firstRecGREG(j,:)=gregorian(firstRecJD);
            lastRecGREG(j,:)=gregorian(lastRecJD);

            % combine mean, max, min and std into columns for all variables per
            % instrument
            meanRec(j) = data.([vname 'mean']);
            stdRec(j) = data.([vname 'std']);
            minRec(j) = data.([vname 'min']);
            maxRec(j) = data.([vname 'max']);

            short_year1=num2str(firstRecGREG(j,1));
            short_year2=num2str(lastRecGREG(j,1));
            short_year1=short_year1(3:4);
            short_year2=short_year2(3:4);
            % write data to file
            fprintf(fid,' %8s  %s  %02.0f/%02.0f/%s %02.0f:%02.0f  %02.0f/%02.0f/%s %02.0f:%02.0f  %7.0f  %6.1f  %6.1f  %6.1f  %6.1f\n',...
                char(SN{j}),variables(j,:),...
                firstRecGREG(j,3),firstRecGREG(j,2), ...
                short_year1,firstRecGREG(j,4),firstRecGREG(j,5),...
                lastRecGREG(j,3),lastRecGREG(j,2),short_year2,lastRecGREG(j,4),lastRecGREG(j,5),...
                numRecs(j),meanRec(j),stdRec(j),minRec(j),maxRec(j));
        else
            fprintf(fid,' %8s  %s  No valid data \n',char(SN{j}),variables(j,:));
        end
    end

    % print dividing line between instruments
    fprintf(fid,'%s\n','-----------------------------------------------------------------------------------------');
end
fclose(fid);

if non_verbose==1
    fn = fieldnames(alldata);
    ii = find(strcmp(fn,'MC') | strcmp(fn,'ODOMC'));
    deletefields={'p','t','c','jd'};
    for no = 1:length(ii)
        d = alldata.(fn{ii(no)});
        d = rmfield(d, deletefields);
        dataout(no) = d;
        dataout(no).sn = id_z_sn.sn(ii(no));
    end
end

