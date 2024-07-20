% function seaphox2rodb_01(moor,'outpath',outpath,'procpath',procpath,'cruise',cruise,'p_coeffs',p_coeffs)
% Function to convert Satlantic Deep SeapHOx data to rodb format.
%
% required inputs: moor - mooring name e.g. 'wb2_3_200606' or ctd cast e.g.
% 'cast5' 
%
% Requires a seperate text file containing the SeapHOx filenames as this
% is based on the day the isntrument started logging
% Format of required text file is serial number in first column and 
% filename in second, either space or tab delimited.
% Naming of text file is as eb3_1_200407_filenames.txt
% Filenames will all have .CSV extension but prefix varies depending on
% start of sampling. Will also be associated .CTD file. Copying folders of
% data per instrument so as to not copy over things with same filenames.
% ?????_filenames.txt file should be in same directory as folders of raw 
% files - e.g.
% .../rpdmoc/rapid/raw/dy039/seaphox_cal_dip/cast5/?????_filenames.txt
%
% optional inputs: procpath - if not using standard paths for info.dat file
%                  outpath - if not using standard paths 
%                  cruise - if running on historical cruise, otherwise will
%                  take cruise name from mfilename path. Also used if
%                  running for CTD cal dip so knows where to loog for data
%                  files
%                  p_coeffs - recalculating the pH from the raw sensor
%                  voltages requires the pressure correction coefficients
%                  as supplied with the instrument calibration. These are
%                  not output in the standard data file from the version of
%                  the SeaFET software that shipped with the instruments in
%                  2015. If not supplied then these will be prompted for if
%                  choosing to recalculate the pH.
%                  p_coeffs in format [p1, p2, p3, p4, p5, p6] with one row 
%                  per instrument on mooring being processed - k0 and k2
%                  not needed as these are stored in the data file
%
% functions called: hms2h.m
%                   julian.m
%                   rodbload.m
%
% Rayner, 07 November 2015 - converted from nortek2rodb_01
% LOH, May 2017 - Add plot for each individual instrument
% LOH, July 2018 - Conversion of oxygen from ml/L to micromol/kg
% LAD, Oct 2020 inserted new code for reading output from seafets that have
% been updated with new software 
% LAD, 0ct 2020 udated to print .png 
% LAD Oct, 2020 updated Xlim and Ylim
% TSD Jul, 2024 - updated to use global MOORPROC_G 
%               - Added new pieces of code of D. Rayner to read Version 1
%                 and 2: (D. Rayner -> 14 March 2024 - updated to allow 
%                           conversion of V2 SeaFET data, which should now 
%                           be used but keeping capability in case was to 
%                           re-process legacy v1 data)

function seaphox2rodb_01(moor, varargin)

%global basedir
%
%if exist('basedir','var')==0
%    if strfind(mfilename('fullpath'),'/Volumes/rpdmoc/rapid/data/exec/')
%        % using Mac with mount to rpdmoc either on a cruise or at NOC
%        basedir = '/Volumes/rpdmoc/rapid/data/';
%    % else % using NOC network
%    %     basedir = '/home/mstar/osnap/data/';
%    else % on pstar machine
%        basedir ='/local/users/pstar/osnap/'
%    end
%else
%   
%end

global MOORPROC_G

basedir = MOORPROC_G.moordatadir;
cruise= MOORPROC_G.cruise;

if nargin==0
    help seaphox2rodb_01
    return
end

%if nargin==0
%    help(mfilename)
%    return
%end

% check for optional arguments
if strcmpi(moor(1:4),'cast')
    % indiciates is a CTD cal dip
    cal_dip=1;
else
    cal_dip=0;
end

a=strmatch('procpath',varargin,'exact');
if a>0
    procpath=char(varargin(a+1));
else
    procpath=[basedir '/proc'];
end

a=strmatch('cruise',varargin,'exact');
if a>0
    cruise=char(varargin(a+1));
else
    script=mfilename('fullpath');
    cruise=script(strfind(script,'/data/exec/')+11:strfind(script,'/stage1/seaphox/')-1);
end

a=strmatch('inpath',varargin,'exact');
if a>0
    inpath=char(varargin(a+1));
else
    if cal_dip==1
        inpath=[basedir '/raw/' cruise '/seaphox_cal_dip/' moor '/'];
    else
        inpath=[basedir '/raw/' cruise '/seaphox/'];
    end
end

a=strmatch('outpath',varargin,'exact');
if a>0
    outpath=char(varargin(a+1));
else
    if cal_dip==1
        outpath=[basedir '/proc_calib/' cruise '/cal_dip/seaphox/' moor '/'];
    else
        outpath=[procpath '/' moor '/seaphox/'];
    end
end

a=strmatch('p_coeffs',varargin,'exact');
if a>0
    p_coeffs=varargin{a+1};
else
    p_coeffs=[];
end

% --- get moring information from infofile 
if cal_dip==1
    infofile =[basedir '/proc_calib/' cruise '/cal_dip/' moor 'info.dat'];
else
    infofile =[procpath '/' moor '/' moor 'info.dat'];
end

[gash, operator]=system('whoami'); % This line will not work if run from a PC. May need to edit it out.
% % gash is not used, but a second variable needs to be specified for the system command


% ----- read infofile / open logfile  ------------------------------------

infovar = 'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth'; 
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd]  =  rodbload(infofile,infovar);

disp(['Writing to: ',outpath,moor,'_seaphox_stage1.log'])
fidlog   = fopen([outpath,moor,'_seaphox_stage1.log'],'a');
fprintf(fidlog,'Transformation of SeapHOx data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));
fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n\n\n',lon);

bg = (julian([s_d(:)' hms2h([s_t(:)' 0])])); %start
ed = julian([e_d(:)' hms2h([e_t(:)' 0])]); %end
% bgvec = ([s_d(:)' ([s_t(:)' 0])]); % ---- D. RAYNERs version do not use this
% bg = datenum(bgvec); %start
% edvec = ([e_d(:)' ([e_t(:)' 0])]); %end
% ed = datenum(edvec);

vec=find(id==375); % SeapHOx instrument code

serial_nums=sn(vec);

% load text file containing serial numbers with related filenames. NB: No
% header in this file, but column 1 should be serial number and column 2
% should be filename. This is necessary as there is currently no standard
% naming convention for Argonaut files. Do not include path though.
% If have multiple files include wildcard and extension and routine will
% search all matching files in that directory
% e.g. 104/DATA/*.CSV

textfile=[inpath moor '_filenames.txt'];
%check if file exists
if ~exist(textfile)
    disp(textfile)
    disp('File containing filenames does not exist - stopping routine.')
end
[input_ser_nums filenames]=textread(textfile,'%d %s');

% -------- load data --------------
for i = 1:length(filenames)
    rawfiles=eval(['dir(''' basedir  '/raw/' cruise '/seaphox/' filenames{i} ''')']);
    rawpath=[basedir  '/raw/' cruise '/seaphox/' filenames{i}(1:strfind(filenames{i},'*')-1)];
    % remove files that start ._????? as the ._ is a remnant from copying
    % from OSX and is a hidden file, not valid data.
    valid_file=[];
    for j=1:length(rawfiles)
        if isempty(strfind(rawfiles(j).name,'._'))
            valid_file=[valid_file j];
        end
    end
    
    for j=1:length(valid_file)
        infile=[rawpath char(rawfiles(valid_file(j)).name)];
        fprintf(fidlog,'infile : %s\n',infile);
        fid=fopen(infile,'r');
        % check whether V1 or V2 format file
        line1=fgetl(fid);
        if contains(line1,'SATFHR.APP-VERSION')
            version=1;
        elseif contains(line1,'FrameSync')
            version=2;
        else
           disp('File format unknown. Doesn''t seem to be V1 or V2 SeapHOx - stopping')
           dbstop if error
           return
        end
        % close an reopen file to reset pointer
        fclose(fid);
        fid=fopen(infile,'r');
        if version==1
            % run through header first
            for k=1:8
                eval(['header' num2str(k) '=textscan(fid,''%s %s %s'',1,''delimiter'','','');'])
            end
            part_data=textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %f','delimiter',',');
            if j==1
                all_data=part_data;
            else
                for k=1:length(all_data)
                    all_data{k}=[all_data{k};part_data{k}];
                end
            end
        elseif version==2
            header1=textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s',1,'delimiter',',');
            part_data=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');
            if j==1
                all_data=part_data;
            else
                for k=1:length(all_data)
                    all_data{k}=[all_data{k};part_data{k}];
                end
            end
        end
        fclose(fid);
    end
    
    outfile=[moor '_' sprintf('%03d',input_ser_nums(i)) '.raw'];
    indep=find(sn==input_ser_nums(i));
    indep=z(indep);
    if length(indep)>1
        indep=indep(1);
    end
    
    fprintf(fidlog,'outfile: %s\n',outfile);
    fprintf(fidlog,'SeapHOx serial number  : %d\n',input_ser_nums(i));
  

    if version==1;
        date=all_data{2};
        year=floor(date/1000); doy=date-year*1000; 
        hour=all_data{3};
        ph_ext=all_data{5}; % pH value for Deep SeaFETs is pHExt as there is no internal pH sensor as per shallow SeaFETs
                            % The data should be converted from the raw voltages to
                            % a pH value, but this was not happening sometimes on
                            % DY039 with the reason not clear from discussions with
                            % the manufacturer. However the raw sensor voltage is
                            % recorded so this is used to calculate the pH.
        mc_t=all_data{7};
        mc_s=all_data{8};
        mc_o_mlL=all_data{9};
        mc_p=all_data{10};
        v_ext=all_data{12}; % raw sensor voltage from external (ie. the only one for Deep Seafet) pH sensor 
        main_v=all_data{18}; % main battery voltage
        isol_v=all_data{20}; % isolated battery voltage

        % conversion from ml/L to micromol/kg
        % Molar Volume of O2 at STP = 22.391 L
        % 1 micromol O2 = 0.022391 mL
        % 1 ml/L = 44.661 micromol/L
        % 1 micromol/kg = 1000/(sw_dens(s,t,p) micromol/L
        mc_o = mc_o_mlL.*44.661./(sw_dens(mc_s,mc_t,mc_p)./1000);

        [DD MM YY]=DOY2date(doy,year);
    
    % 	% DO NOT USE THIS CALCULATION FROM THE SHALLOW SEAFET MANUAL AS THERE
    % 	IS NO INCLUSION OF THE PRESSURE CORRECTION - INSTEAD USE THE SEABIRD
    % 	SUPPLIED PHCALC.M ROUTINE
    %     calculate pH from sensor voltage as not always output by instrument
    %     % formulae from SeaFET manual
    %     MC_T=mc_t+273.15; % convert deg C to Kelvin for use in calculations
    %     I=19.924 .* mc_s ./ (1000 - 1.005 .* mc_s);
    %     ADH=0.00000343 * mc_t .^2 + 0.00067524 * mc_t + 0.49172143; % uses temperature in deg C, not Kelvin
    %     logLhcl=((-ADH .*I .^0.5) ./ (1 + 1.1394 * I .^0.5)) + (0.08885 - 0.000111 * mc_t) .* I; % uses temperature in deg C, not Kelvin
    %     Ks=(1-0.001005.*mc_s) .* exp((-4276.1./MC_T) + 141.328 - 23.093.*log(MC_T) + (-13856./MC_T + 324.57 - 47.986*log(MC_T)).* I.^0.5 + (35474./MC_T - 771.54 + 114.723.*log(MC_T)) .* I - (2698./MC_T).*I.^1.5 + (1776./MC_T).*I.^2);
    %     St=(0.1400/96.062)*(mc_s/1.80655);
    %     Clt=0.99889/35.453*mc_s/1.80655;
    %     R=8.314472;
    %     F=96485.3415;
    %     Snernst=(R.*(MC_T).*log(10))./F;
        if strcmpi(header5{2},'CAL_PHEXT_OFFSET_COEFF')
            k0e=str2double(header5{3});
        else
            disp('pH_Ext Offset coefficient not in expected position of 5th header line. Code needs improving to search header lines. STOPPING!');
            dbstop if error
            return
        end
        if strcmpi(header6{2},'CAL_PHEXT_SLOPE_COEFF')
            k2e=str2double(header6{3});
        else
            disp('pH_Ext Slope coefficient not in expected position of 5th header line. Code needs improving to search header lines. STOPPING!');
            dbstop if error
            return
        end
    
        % % line below as per Shallow SeaFET manual - gives erroneous pH readings
        % ph_recalc=(v_ext-k0e-k2e.*MC_T)./Snernst + log10(Clt) + 2* logLhcl - log10(1+St./Ks);
        % line below as per Bresnahan et al 2014 Methods in Oceanography paper.
        % But doesn't include a pressure correction so not actually correct.
        % % ph_recalc_without_p_correction=(v_ext-k0e-k2e.*(mc_t-25))./Snernst + log10(Clt) + 2* logLhcl - log10(1+St./Ks);
    
    
        if find(~isnan(ph_ext)) % implies pH was being output by the instrument when sometimes it is now - reason not clarified by SeaBird/Satlantic yet 
            response_valid=0;
            j=1;
            while response_valid~=1;
                if j>1
                    disp('Invalid response!')
                end
                use_ph_ext=input('Values detected for ph_ext in datafile. Do you want to use these for pH (u), or recalculate from raw voltages (r)? u/r?:- ','s');
                if strcmpi(use_ph_ext,'u')
                    recalc_ph=0;
                response_valid=1;
                elseif strcmpi(use_ph_ext,'r')
                    recalc_ph=1;
                    response_valid=1;
                end
                j=j+1;
            end
        else
            disp('No values detected for pH External - using values as recalculated from raw voltages')
            recalc_ph=1;
        end

        if recalc_ph==1
            % recalculate pH from raw sensor voltages using the phcalc.m routine supplied by seabird when on DY039
            % check if p_coeffs entered to function
            if ~isempty(p_coeffs)
                p_coeffs_inst=p_coeffs(i,:); % assign coeffs for correct serial number
            else
                disp('pressure coefficients not provided to function, please enter here from calibration certificate:')
                disp(['Serial number: ' num2str(input_ser_nums(i))])
                p1=input('P1 = ');
                p2=input('P2 = ');
                p3=input('P3 = ');
                p4=input('P4 = ');
                p5=input('P5 = ');
                p6=input('P6 = ');
                p_coeffs_inst=[p1 p2 p3 p4 p5 p6];
            end
            [ph_free,ph]=phcalc(v_ext, mc_p, mc_t, mc_s, [k0e k2e p_coeffs_inst]);
        else 
            ph=ph_ext;
        end
    elseif version ==2
        date_time=datenum(all_data{2});
        date_time_vec=datevec(date_time);
        YY=date_time_vec(:,1);
        MM=date_time_vec(:,2);
        DD=date_time_vec(:,3);
        hour=date_time_vec(:,4)+date_time_vec(:,5)/60+date_time_vec(:,6)/60/60;
        
        mc_t=all_data{5};
        ph_ext=all_data{6}; % pH value for Deep SeaFETs is pHExt as there is no internal pH sensor as per shallow SeaFETs
        v_ext=all_data{7}; % raw sensor voltage from external (ie. the only one for Deep Seafet) pH sensor 
        
        % assuming V2 has fixed bug where sometimes didn't output pH value
        % from raw voltages so not going to recalculate here. If need to
        % put it in in the future then see the code under the verrsion==1
        % check
        ph=ph_ext;
        
        mc_p=all_data{8};
        mc_s=all_data{9};
        mc_c=all_data{10};
        mc_o_mlL=all_data{11};

        % conversion from ml/L to micromol/kg
        % Molar Volume of O2 at STP = 22.391 L
        % 1 micromol O2 = 0.022391 mL
        % 1 ml/L = 44.661 micromol/L
        % 1 micromol/kg = 1000/(sw_dens(s,t,p) micromol/L
        mc_o = mc_o_mlL.*44.661./(sw_dens(mc_s,mc_t,mc_p)./1000);
    
    end
                
       
    dat       = [YY MM DD hour];
    jd        = julian(dat);
    valI      = find(jd<ed & jd>bg);

 
    % ----- save data to rodb -----------------
    %NB oxygen saved in ml/l (so "O" variable rather than "O2")

    if find(intersect(input_ser_nums(i),serial_nums))
        if version==1
            columns = 'YY:MM:DD:HH:PH:PH_V:T:P:O:S:BATT1:BATT2';
            data = [dat ph v_ext mc_t mc_p mc_o mc_s main_v isol_v]; 
            fort =['%4.4d %2.2d %2.2d  %6.4f  %7.5f %10.8f  %6.4f %7.3f %6.3f  %6.4f %5.3f %5.3f'];
        elseif version==2
            columns = 'YY:MM:DD:HH:PH:PH_V:T:C:P:O:S';
            data = [dat ph v_ext mc_t mc_c*10 mc_p mc_o mc_s]; 
            fort =['%4.4d %2.2d %2.2d  %6.4f  %6.4f %8.6f  %6.4f %6.4f %7.3f  %5.3f %6.4f'];
        end
        infovar = ['Mooring:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:' ...
                   'Columns:SerialNumber:InstrDepth']; 
        rodbsave([outpath outfile],infovar,fort,moor,s_t,s_d,e_t,e_d,lat,lon,wd,columns,...
                 input_ser_nums(i),indep,data);
        eval(['disp(''Data written to ' outpath outfile ''')']);

    else
        disp(['Serial number does not match those in info.dat file - sn: ' num2str(input_ser_nums(i))])
        disp('Stopping routine')
        return
    end
    
    % -------- generate logfile entries --------------
    
    sz   =   size(data);
    
    fprintf(fidlog,'Instrument Target Depth [m] : %d\n',indep);
    fprintf(fidlog,'Start date and time         : %s \n',datestr(datenum(gregorian(jd(1)))));
    fprintf(fidlog,'End date and time           : %s \n',datestr(datenum(gregorian(jd(end)))));
    sampling_rate = round(1./median(diff(jd)));
    ex_samples = round((jd(end)-jd(1))*sampling_rate+1);
    fprintf(fidlog,'Sampling Frequency [per day]: %d \n',sampling_rate);
    fprintf(fidlog,'Number of samples           : %d; expected: %d \n',sz(1),ex_samples);
    
    ph2=ph(valI,:);
    t2 =mc_t(valI,:);
    s2 =mc_s(valI,:);
    o2 =mc_o(valI,:);
    p2 =mc_p(valI,:);

    fprintf(fidlog,['Oxygen conversion: conversion from ml/L to micromol/kg \n ...' ...
        'Molar Volume of O2 at STP = 22.391 L \n ...' ...
        '1 micromol O2 = 0.022391 mL \n ...' ...
        '1 ml/L = 44.661 micromol/L \n ...' ...
        '1 micromol/kg = 1000/(sw_dens(s,t,p) micromol/L \n']);
    
    m_ph = median(ph2(~isnan(ph2)));
    m_t = median(t2(~isnan(t2)));
    m_s = median(s2(~isnan(s2)));
    m_o = median(o2(~isnan(o2)));
    m_p = median(p2(~isnan(p2)));
    
    fprintf(fidlog,'Median pH                             : %4.2f \n',m_ph);
    fprintf(fidlog,'Median MicroCAT temperature [deg C]   : %4.2f \n',m_t);
    fprintf(fidlog,'Median MicroCAT salinity              : %4.2f \n',m_s);
    fprintf(fidlog,'Median MicroCAT oxygen [ml/l]         : %5.2f \n',m_o);
    fprintf(fidlog,'Median MicroCAT pressure [dbar]       : %6.2f \n\n',m_p);
       
    %fclose(fid);
end  % loop of serial numbers
fclose(fidlog);

% TODO: Add CTD data on plots with all SeapHOx (if there is several on a cast) 
%end
