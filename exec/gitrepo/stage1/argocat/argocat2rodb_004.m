% function argocat2rodb_004(moor,'inpath','outpath','procpath')
% Function to convert Sontek Argonaut data to rodb format.
% Works for both stand alone instruments and those coupled with Microcats
%
% Requires a seperate text file containing the Argonaut filenames as this
% is non conventional. Format of required text file is serial number in
% first column and filename in second, either space or tab delimited.
% Naming of text file is as eb3_1_200407_filenames.txt
% Filenames will all have .dat extension but prefix varies depending on
% setup.
% filenames.txt file should be in same directory as raw files.
%
% required inputs: moor - mooring name e.g. 'wb2_3_200606'
%
% optional inputs: procpath - if not using standard paths for info.dat file
%                  inpath - if not using standard paths
%                  outpath - if not using standard paths
%
% functions called: hms2h.m
%                   julian.m
%                   rodbload.m
%

% Kanzow, 13 April 2005
% Changes - 11/9/06  - Included noise, signal 2 noise ratio and percent good
%                      pings into output format. DR.
%         - 11/12/06 - adapted to a function with paths options. DR.
%         - 12/12/06 - adapted so deals with both Argonauts paired with
%                      Microcats and those that are not. DR.
%         - 11/01/07 - fixed glitch where wasn't passing outpath to
%                      rodbsave. DR.
%         - 21/04/11 - fixed bug where data format changed slightly for
%                      files produced in ViewArgonaut version 3.71. New
%                      format has 2 less columns in it (CellBegin and
%                      CellEnd are missing)
% 
%         - 12/09/11 - edited paths and added cruise variable - gdm jc064

function argocat2rodb_004(moor, varargin)

cruise = 'd382';

if nargin==0
    help argocat2rodb_04
    return
end

% check for optional arguments
a=strmatch('procpath',varargin,'exact');
if a>0
    procpath=char(varargin(a+1));
else
    procpath='/noc/users/pstar/rpdmoc/rapid/data/moor/proc/';
end

a=strmatch('inpath',varargin,'exact');
if a>0
    inpath=char(varargin(a+1));
else
    inpath=['/noc/users/pstar/rpdmoc/rapid/data/moor/raw/' cruise '/arg/'];
end

a=strmatch('outpath',varargin,'exact');
if a>0
    outpath=char(varargin(a+1));
else
    outpath = [procpath moor '/arg/'];
end


% --- get moring information from infofile 
infofile = [procpath '/' moor '/' moor 'info.dat'];


[gash, operator]=system('whoami'); % This line will not work if run from a PC. May need to edit it out.
% % gash is not used, but a second variable needs to be specified for the system command


% ----- read infofile / open logfile  ------------------------------------

infovar = 'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth'; 
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd]  =  rodbload(infofile,infovar);

fidlog   = fopen([outpath,moor,'_Argonaut_stage1_log'],'w');
fprintf(fidlog,'Transformation of ascii data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));
fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n\n\n',lon);

bg = julian([s_d(:)' hms2h([s_t(:)' 0])]); %start
ed = julian([e_d(:)' hms2h([e_t(:)' 0])]); %end

combo=find(id==366337); % Argonauts coupled with Microcats 
individual=find(id==366); % uncoupled Argonauts

individual_serial_nums=sn(individual);

for i=1:length(individual)
    for j=1:length(combo)
        if individual_serial_nums(i)==sn(combo(j))
            individual_serial_nums(i)=NaN;
        end
    end
end

% load text file containing serial numbers with related filenames. NB: No
% header in this file, but column 1 should be serial number and column 2
% should be filename. This is necessary as there is currently no standard
% naming convention for Argonaut files. Do not include path though.

textfile=[inpath moor '_filenames.txt'];
%check if file exists
if ~exist(textfile)
    disp(textfile)
    disp('File containing filenames does not exist - stopping routine.')
end
[input_ser_nums filenames]=textread(textfile,'%d %s');

% -------- load data --------------
for i = 1:length(filenames)
    infile=[inpath char(filenames(i))]
    outfile=[moor '_' num2str(input_ser_nums(i)) '.raw'];
    indep=find(sn==input_ser_nums(i));
    indep=z(indep);
    if length(indep)>1
        indep=indep(1);
    end
    fprintf(fidlog,'infile : %s\n',infile);
    fprintf(fidlog,'outfile: %s\n',outfile);
    fprintf(fidlog,'Argonaut serial number  : %d\n',input_ser_nums(i));
    
    if find(intersect(input_ser_nums(i),sn(combo)));
        a=find(sn==input_ser_nums(i));
        b=z(a);
        if length(b)>1
            b=b(1);
        end
        c=find(z==b & sn~=input_ser_nums(i));
        mc_ser_num=sn(c);
        fprintf(fidlog,'Argonaut is coupled with Microcat serial number: %d\n',mc_ser_num);
    end
    
    fid = fopen(infile,'r');

    readline = fscanf(fid,'%c');
    fclose(fid);

    ret = sprintf('\n');
    retx = findstr(readline,ret);  % car. return indices

    % convert string to numbers

    data      = str2num(readline(retx(1)+1:end));
    dat       = [data(:,[1:3]) hms2h(data(:,4:6))];
    jd        = julian(dat);
    valI      = find(jd<ed & jd>bg);
    uvw       = data(:,7:9);
    uvw_sd    = data(:,10:12); % velocity standard deviation
    uvw_snr   =  data(:,13:15); %signal to noise ratio
    uvw_ss    = data(:,16:18); % signal strength
    uvw_noise = data(:,19:21); % signal noise

    pgp     =  data(:,22); % percent good pings
    hpr     = data(:,23:25); %heading pitch roll
    hpr_sd  = data(:,26:28); %heading pitch roll stand. dev.

    t     = data(:,29);
    p     = data(:,30);
    p_sd  = data(:,31);
    volt  = data(:,32);
    
    if strfind(readline,'CellBegin') % processed using ViewArgonut software version <3.71
        cell  = data(:,33:34);
        vel   = data(:,35);
        dir   = data(:,36);
        if find(intersect(input_ser_nums(i),sn(combo)))
            t_mc  = data(:,37);
            c_mc  = data(:,38);
            p_mc  = data(:,39);
        end
    else
        %cell  = data(:,33:34); % no cell begin and end values in data if
        %processed in ViewArgonaut software 3.71
        vel   = data(:,33);
        dir   = data(:,34);
        if find(intersect(input_ser_nums(i),sn(combo)))
            t_mc  = data(:,35);
            c_mc  = data(:,36);
            p_mc  = data(:,37);
        end
    end
    % ----- save data to rdb -----------------
    if find(intersect(input_ser_nums(i),sn(combo)))
        % Microcat attached so include in data output
        columns = 'YY:MM:DD:HH:T:TCAT:P:PCAT:C:U:V:W:HDG:PIT:ROL:USD:VSD:WSD:USS:VSS:WSS:HDGSD:PITSD:ROLSD:IPOW:UNOISE:VNOISE:WNOISE:USNR:VSNR:WSNR:PGP';
        data = [dat t t_mc p p_mc c_mc*10 uvw hpr uvw_sd  uvw_ss hpr_sd volt uvw_noise uvw_snr pgp]; 
        fort =['%4.4d %2.2d %2.2d  %6.4f    %7.4f %7.4f %6.2f %6.2f %7.4f   %6.2f %6.2f %6.2f   %4.1f %4.1f %4.1f   %6.2f %6.2f %6.2f   %3.3d %3.3d %3.3d   %4.1f %4.1f %4.1f   %4.2f   %3.3d %3.3d %3.3d   %3.1f %3.1f %3.1f   %3.3d'];
        infovar = ['Mooring:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:' ...
                   'Columns:SerialNumber:InstrDepth:MicrocatSN']; 
        rodbsave([outpath outfile],infovar,fort,moor,s_t,s_d,e_t,e_d,lat,lon,wd,columns,...
                 input_ser_nums(i),indep,mc_ser_num,data);


    elseif find(intersect(input_ser_nums(i),individual_serial_nums))
        % No Microcat attached
        columns = 'YY:MM:DD:HH:T:P:U:V:W:HDG:PIT:ROL:USD:VSD:WSD:USS:VSS:WSS:HDGSD:PITSD:ROLSD:IPOW:UNOISE:VNOISE:WNOISE:USNR:VSNR:WSNR:PGP';
        data = [dat t p uvw hpr uvw_sd  uvw_ss hpr_sd volt uvw_noise uvw_snr pgp]; 
        fort =['%4.4d %2.2d %2.2d  %6.4f    %7.4f %6.2f   %6.2f %6.2f %6.2f   %4.1f %4.1f %4.1f   %6.2f %6.2f %6.2f   %3.3d %3.3d %3.3d   %4.1f %4.1f %4.1f   %4.2f   %3.3d %3.3d %3.3d   %3.1f %3.1f %3.1f   %3.3d'];
        infovar = ['Mooring:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:' ...
                   'Columns:SerialNumber:InstrDepth']; 
        rodbsave([outpath outfile],infovar,fort,moor,s_t,s_d,e_t,e_d,lat,lon,wd,columns,...
                 input_ser_nums(i),indep,data);

    else
        disp(['Serial number does not match those in info.dat file - sn: ' num2str(input_ser_nums(i))])
        disp('Stopping routine')
        return
    end

    

     
    % -------- generate logfile entries --------------

    sz   =   size(data);

    fprintf(fidlog,'Instrument Target Depth[m]  : %d\n',indep);
    fprintf(fidlog,'Start date and time         : %s \n',datestr(gregorian(jd(1))));
    fprintf(fidlog,'End date and time           :   %s \n',datestr(gregorian(jd(end))));
    sampling_rate = round(1./median(diff(jd)));
    ex_samples = round((jd(end)-jd(1))*sampling_rate+1);
    fprintf(fidlog,'Sampling Frequency [per day]: %d \n',sampling_rate);
    fprintf(fidlog,'Number of samples           : %d; expected: %d \n',sz(1),ex_samples);


    m_uvw = median(uvw(valI,:));
    m_uvwsd = median(uvw_sd(valI,:));
    m_uvwss = median(uvw_ss(valI,:));
    m_uvwsnr = median(uvw_snr(valI,:));
    m_uvwnoise = median(uvw_noise(valI,:));

    m_pgp = median(pgp(valI));
    m_hpr = median(hpr(valI,:));
    m_hprsd = median(hpr_sd(valI,:));
    m_psd = median(p_sd(valI));
    m_volt = median(volt(valI));
    m_vel = median(vel(valI));
    m_dir = median(dir(valI));
    
    if find(intersect(input_ser_nums(i),sn(combo)))
        m_t = median([t(valI) t_mc(valI)]);
        m_p = median([p(valI) p_mc(valI)]);
        m_c = median(c_mc(valI)*10);
        fprintf(fidlog,'Median temperature Argonaut/Microcat [deg C]        : %5.2f / %5.2f\n',m_t);
        fprintf(fidlog,'Median MicroCAT conductivity [mS/cm]                : %5.2f\n',m_c);
        fprintf(fidlog,'Median pressure Argonaut/Microcat [dbar]            : %6.2f / %6.2f\n',m_p);
    else
        m_t = median([t(valI)]);
        m_p = median([p(valI)]);
        fprintf(fidlog,'Median temperature Argonaut [deg C]                 : %5.2f \n',m_t);
        fprintf(fidlog,'Median pressure Argonaut [dbar]                     : %6.2f \n',m_p);
    end
    
    fprintf(fidlog,'Median velocity u / v / w [cm/s]                    : %4.1f  %4.1f  %4.1f\n',m_uvw);
    fprintf(fidlog,'Median heading / pitch / roll [deg]                 : %4.1f  %4.1f  %4.1f\n',m_hpr);

    fprintf(fidlog,'Median velocity STD u / v / w [cm/s]                : %4.1f  %4.1f  %4.1f\n',m_uvwsd);
    fprintf(fidlog,'Median velocity signal strength u / v / w [count]   : %3.3d  %3.3d  %3.3d\n',m_uvwss);
    fprintf(fidlog,'Median velocity signal to noise ratio u / v / w [??]: %4.1f  %4.1f  %4.1f\n',m_uvwsnr);
    fprintf(fidlog,'Median velocity noise  u / v / w [??]               : %4.1f  %4.1f  %4.1f\n',m_uvwnoise);
    fprintf(fidlog,'Median percent good pings                           : %4.1f \n',m_pgp);
    fprintf(fidlog,'Median heading / pitch / roll STD [deg]             : %4.1f  %4.1f  %4.1f\n',m_hprsd);
    fprintf(fidlog,'Median power input level [V]                        : %4.1f\n\n\n',m_volt);




end  % loop of serial numbers
