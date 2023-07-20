% function nortek2rodb_01(moor,'inpath','outpath','procpath')
% Function to convert Nortek Aquadopp Single Point CM data to rodb format.
%
% Requires a seperate text file containing the Nortek filenames as this
% is non conventional. Format of required text file is serial number in
% first column and filename in second, either space or tab delimited.
% Naming of text file is as eb3_1_200407_filenames.txt
% Filenames will all have .dat extension but prefix varies depending on
% setup.
% filenames.txt file should be in same directory as raw files.
%
% required inputs: moor - mooring name e.g. 'wb2_3_200606'
%                  inpath - path to raw data files (currently no standard paths)
%
% optional inputs: procpath - if not using standard paths for info.dat file
%                  outpath - if not using standard paths 
%
% functions called: hms2h.m
%                   julian.m
%                   rodbload.m
%

% Rayner, 15 January 2007 - converted from argocat2rodb_004
% DR, 14/11/12 - updated so can set basedir

function nortek2rodb_01(moor, varargin)

global MOORPROC_G

basedir = MOORPROC_G.moordatadir;
cruise= MOORPROC_G.cruise;

if nargin==0
    help nortek2rodb_01
    return
end

% check for optional arguments
a=strmatch('procpath',varargin,'exact');
if a>0
    procpath=char(varargin(a+1));
else
    procpath = fullfile(MOORPROC_G.moordatadir,'proc');
end

a=strmatch('inpath',varargin,'exact');
if a>0
    inpath=char(varargin(a+1));
else
    inpath = fullfile(MOORPROC_G.moordadtadir,'raw',cruise,'nor');end

a=strmatch('outpath',varargin,'exact');
if a>0
    outpath=char(varargin(a+1));
else
    outpath = fullfile(procpath,moor,'nor');
end


% --- get moring information from infofile 
infofile = fullfile(procpath,moor,[moor 'info.dat']);

operator = MOORPROC_G.operator;

% ----- read infofile / open logfile  ------------------------------------

infovar = 'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth'; 
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd]  =  rodbload(infofile,infovar);

fidlog   = fopen([outpath,moor,'_Nortek_stage1.log'],'a');
fprintf(fidlog,'Transformation of Nortek ascii data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));
fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n\n\n',lon);

bg = julian([s_d(:)' hms2h([s_t(:)' 0])]); %start
ed = julian([e_d(:)' hms2h([e_t(:)' 0])]); %end

vec=find((id==368|id==370)); % Nortek instrument code

serial_nums=sn(vec)

% load text file containing serial numbers with related filenames. NB: No
% header in this file, but column 1 should be serial number and column 2
% should be filename. This is necessary as there is currently no standard
% naming convention for Argonaut files. Do not include path though.

textfile=fullfile(inpath,[moor '_filenames.txt']);
%check if file exists
if ~exist(textfile)
    disp(textfile)
    error('File containing filenames does not exist - stopping routine.')

end
[input_ser_nums filenames]=textread(textfile,'%d %s');

% -------- load data --------------
for i = 1:length(filenames)
    infile=fullfile(inpath,filenames{i});
    outfile=[moor '_' num2str(input_ser_nums(i)) '.raw'];
    indep=sn==input_ser_nums(i);
    indep=z(indep);
    if length(indep)>1
        indep=indep(1);
    end
    fprintf(fidlog,'infile : %s\n',infile);
    fprintf(fidlog,'outfile: %s\n',outfile);
    fprintf(fidlog,'Nortek serial number  : %d\n',input_ser_nums(i));

    if ~exist(infile,'file')
        checkans = input(['File for serial number ' serial_nums(i) 'not found. Do you want to continue? y/n [y]:'],'s');
       if isempty(checkans)
          reply = 'Y';
       end
       continue
    else 
        disp(['Processing file ' serial_nums(i)])
    end
    
    all_data=load(infile);
    month=all_data(:,1); day=all_data(:,2); year=all_data(:,3); 
    hour=all_data(:,4); minute=all_data(:,5); second=all_data(:,6);
    errcode=all_data(:,7); statcode=all_data(:,8);
    u=all_data(:,9); v=all_data(:,10); w=all_data(:,11);
    Amp1=all_data(:,12); Amp2=all_data(:,13); Amp3=all_data(:,14);
    bat=all_data(:,15); 
    
    % updated for differences in file format depending on which version of the software used to download and convert - dr400
    % cruise kn200-4
    a=size(all_data);
    if a(2)==27 % I THINK THIS IS RIGHT FOR USING "SOUNDSPEED USED" rather than just SOUNDSPEED for newer Norteks. Need to check really.
        disp(' I THINK THIS IS RIGHT FOR USING "SOUNDSPEED USED" rather than just SOUNDSPEED for newer Norteks. Need to check really.')
        soundspeed=all_data(:,17); heading=all_data(:,18); pitch=all_data(:,19); roll=all_data(:,20);
        pressure=all_data(:,21); temperature=all_data(:,23);
        analogue1=all_data(:,24); analogue2=all_data(:,25);
        speed=all_data(:,26); direction=all_data(:,27);
    elseif a(2)==25 % converted in older software
        soundspeed=all_data(:,16);  heading=all_data(:,17); pitch=all_data(:,18); roll=all_data(:,19);
        pressure=all_data(:,20); temperature=all_data(:,21);
        analogue1=all_data(:,22); analogue2=all_data(:,23);
        speed=all_data(:,24); direction=all_data(:,25);
    end

    dat       = [year month day hms2h(hour,minute,second)];
    jd        = julian(dat);
    valI      = find(jd<ed & jd>bg);
    
    % ----- save data to rodb -----------------
    

    if find(intersect(input_ser_nums(i),serial_nums))
        columns = 'YY:MM:DD:HH:T:P:U:V:W:HDG:PIT:ROL:USS:VSS:WSS:IPOW:CS:CD';
        % errcode, statcode analogue1, analogue2 and soundspeed not saved in this version.
        data = [dat temperature pressure u.*100 v.*100 w.*100 heading pitch roll Amp1 Amp2 Amp3 bat speed.*100 direction]; 
        fort =['%4.4d %2.2d %2.2d  %6.4f  %4.2f %7.3f  %4.1f %4.1f %4.1f  %4.1f %4.1f %4.1f  %2d %2d %2d  %2.1f  %4.1f %5.2f'];
        infovar = ['Mooring:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:' ...
                   'Columns:SerialNumber:InstrDepth']; 
        rodbsave(fullfile(outpath,outfile),infovar,fort,moor,s_t,s_d,e_t,e_d,lat,lon,wd,columns,...
                 input_ser_nums(i),indep,data);
        eval(['disp(''Data written to ' fullfile(outpath,outfile) ''')']);

    else
        disp(['Serial number does not match those in info.dat file - sn: ' num2str(input_ser_nums(i))])
        disp('Stopping routine')
        return
    end

    

     
    % -------- generate logfile entries --------------

    sz   =   size(data);

    fprintf(fidlog,'Instrument Target Depth[m]  : %d\n',indep);
    fprintf(fidlog,'Start date and time         : %s \n',datestr(gregorian(jd(1))));
    fprintf(fidlog,'End date and time           : %s \n',datestr(gregorian(jd(end))));
    sampling_rate = round(1./median(diff(jd)));
    ex_samples = round((jd(end)-jd(1))*sampling_rate+1);
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
