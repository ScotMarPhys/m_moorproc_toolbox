% function nortek2rodb_01(moor)
% function nortek2rodb_01(moor, pd)
% Function to convert Nortek Aquadopp Single Point CM data to rodb format.
%
% Nortek filenames will all have .dat extension but prefix varies depending
% on setup.
% If Nortek filenames do not have one of a few conventional formats
% (serialnum.dat, serialnum_data.dat, serialnum01.dat), they can be read
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

% Rayner, 15 January 2007 - converted from argocat2rodb_004
% DR, 14/11/12 - updated so can set basedir

function nortek2rodb_01(moor, varargin)

global MOORPROC_G

if nargin==0
    help nortek2rodb_01
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

if ~exist(pd.stage1log,'file')
    fp = fileparts(pd.stage1log);
    if ~exist(fp,'dir')
        mkdir(fp)
    end
end
fidlog   = fopen(pd.stage1log,'a');
fprintf(fidlog,'Transformation of Nortek ascii data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));
fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n\n\n',lon);

bg = datenum([s_d(:)' s_t(:)' 0]); %start
ed = datenum([e_d(:)' e_t(:)' 0]); %end

mnor = (id==368 | id==370); % Nortek instrument code
serial_nums=sn(mnor)
z = z(mnor);
id = id(mnor);

%first try to find files with names in a few standard forms
filenames = cell(length(serial_nums),1);
pats = {'%d_data.dat';'%d.dat';'%d01.dat'};
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
    % naming convention for Argonaut files. Do not include path though.
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
    warning('these S/Ns from info file do not have standard names and were not found in listfile:')
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
for i = 1:length(filenames)
    infile=fullfile(pd.rawpath,filenames{i});
    if exist(infile,'file')
        fprintf(1,'Processing file %d\n', serial_nums(i))
    else
       checkans = input(['File for serial number ' serial_nums(i) 'not found. Do you want to continue to next? y/n [y]:'],'s');
       if strcmpi(checkans,'n')
           error('stopping on S/N %d file not found',serial_nums(i))
       else
           %on to next
           continue
       end
    end
    outfile=fullfile(pd.stage1path,sprintf(pd.stage1form,serial_nums(i)));

    fprintf(fidlog,'infile : %s\n',infile);
    fprintf(fidlog,'outfile: %s\n',outfile);
    fprintf(fidlog,'Nortek serial number  : %d\n',serial_nums(i));
    
    all_data=load(infile);
    month=all_data(:,1); day=all_data(:,2); year=all_data(:,3); 
    hour=all_data(:,4); minute=all_data(:,5); second=all_data(:,6);
    %errcode=all_data(:,7); statcode=all_data(:,8);
    u=all_data(:,9); v=all_data(:,10); w=all_data(:,11);
    Amp1=all_data(:,12); Amp2=all_data(:,13); Amp3=all_data(:,14);
    bat=all_data(:,15); 
    
    % updated for differences in file format depending on which version of the software used to download and convert - dr400
    % cruise kn200-4
    a=size(all_data);
    if a(2)==27 % I THINK THIS IS RIGHT FOR USING "SOUNDSPEED USED" rather than just SOUNDSPEED for newer Norteks. Need to check really.
        disp(' I THINK THIS IS RIGHT FOR USING "SOUNDSPEED USED" rather than just SOUNDSPEED for newer Norteks. Need to check really.')
        %soundspeed=all_data(:,17); 
        heading=all_data(:,18); pitch=all_data(:,19); roll=all_data(:,20);
        pressure=all_data(:,21); temperature=all_data(:,23);
        %analogue1=all_data(:,24); analogue2=all_data(:,25);
        speed=all_data(:,26); direction=all_data(:,27);
    elseif a(2)==25 % converted in older software
        %soundspeed=all_data(:,16);  
        heading=all_data(:,17); pitch=all_data(:,18); roll=all_data(:,19);
        pressure=all_data(:,20); temperature=all_data(:,21);
        %analogue1=all_data(:,22); analogue2=all_data(:,23);
        speed=all_data(:,24); direction=all_data(:,25);
    end

    dat       = [year month day hms2h(hour,minute,second)];
    dnum      = datenum([year month day hour minute second]);
    valI      = find(dnum<ed & dnum>bg);
    
    % ----- save data to rodb -----------------
    

    columns = 'YY:MM:DD:HH:T:P:U:V:W:HDG:PIT:ROL:USS:VSS:WSS:IPOW:CS:CD';
    % errcode, statcode analogue1, analogue2 and soundspeed not saved in this version.
    data = [dat temperature pressure u.*100 v.*100 w.*100 heading pitch roll Amp1 Amp2 Amp3 bat speed.*100 direction];
    fort ='%4.4d %2.2d %2.2d  %6.4f  %4.2f %7.3f  %4.1f %4.1f %4.1f  %4.1f %4.1f %4.1f  %2d %2d %2d  %2.1f  %4.1f %5.2f';
    infovar = ['Mooring:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:' ...
        'Columns:SerialNumber:InstrDepth'];
    rodbsave(outfile,infovar,fort,moor,s_t,s_d,e_t,e_d,lat,lon,wd,columns,...
        serial_nums(i),z(i),data);
    fprintf(1,'Data written to %s',outfile)



    % -------- generate logfile entries --------------

    sz   =   size(data);

    fprintf(fidlog,'Instrument Target Depth[m]  : %d\n',z(i));
    fprintf(fidlog,'Start date and time         : %s \n',datestr(dnum(1)));
    fprintf(fidlog,'End date and time           : %s \n',datestr(dnum(end)));
    sampling_rate = round(1./median(diff(dnum)));
    ex_samples = round((dnum(end)-dnum(1))*sampling_rate+1);
    fprintf(fidlog,'Sampling Frequency [per day]: %d \n',sampling_rate);
    fprintf(fidlog,'Number of samples           : %d; expected: %d \n',sz(1),ex_samples);

    m_u = median(u(valI,:));
    m_v = median(v(valI,:));
    m_w = median(w(valI,:));
    m_uss = median(Amp1(valI,:));
    m_vss = median(Amp2(valI,:));
    m_wss = median(Amp3(valI,:));
    m_hdg = median(heading(valI,:));
    m_pit = median(pitch(valI,:));
    m_rol = median(roll(valI,:));
    m_volt = median(bat(valI));
    m_vel = median(speed(valI));
    m_dir = median(direction(valI));
    m_t = median([temperature(valI)]);
    m_p = median([pressure(valI)]);
    
    fprintf(fidlog,'Median temperature Nortek [deg C]                   : %4.2f \n',m_t);
    fprintf(fidlog,'Median pressure Nortek [dbar]                       : %7.3f \n',m_p);
    fprintf(fidlog,'Median velocity u / v / w [cm/s]                    : %4.1f  %4.1f  %4.1f\n',m_u, m_v, m_w);
    fprintf(fidlog,'Median heading / pitch / roll [deg]                 : %4.1f  %4.1f  %4.1f\n',m_hdg, m_pit, m_rol);
    fprintf(fidlog,'Median velocity signal amplitude u / v / w [count]  : %2d  %2d  %2d\n',m_uss, m_vss, m_wss);
    fprintf(fidlog,'Median speed [cm/s]                                 : %4.1f\n',m_vel);
    fprintf(fidlog,'Median direction [cm/s]                             : %5.2f\n',m_dir);
    fprintf(fidlog,'Median power input level [V]                        : %4.1f\n\n\n',m_volt);

end  % loop of serial numbers
fclose(fidlog);
