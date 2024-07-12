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
% .../rpdmoc/rapid/data/moor/raw/dy039/seaphox_cal_dip/cast5/?????_filenames.txt
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
% TSD Jul, 2024 modified to use global MOORPROC_G 

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
    help nortek2rodb_01
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
    procpath=[basedir '\data\moor/proc'];
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
        inpath=[basedir '/data/moor/raw/' cruise '/seaphox_cal_dip/' moor '/'];
    else
        inpath=[basedir '/data/moor/raw/' cruise '/seaphox/'];
    end
end

a=strmatch('outpath',varargin,'exact');
if a>0
    outpath=char(varargin(a+1));
else
    if cal_dip==1
        outpath=[basedir '/data/moor/proc_calib/' cruise '/cal_dip/seaphox/' moor '/'];
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
    infofile =[basedir '/data/moor/proc_calib/' cruise '/cal_dip/' moor 'info.dat'];
else
    infofile =[procpath '/' moor '/' moor 'info.dat'];
end

[gash, operator]=system('whoami'); % This line will not work if run from a PC. May need to edit it out.
% % gash is not used, but a second variable needs to be specified for the system command


% ----- read infofile / open logfile  ------------------------------------

infovar = 'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth'; 
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd]  =  rodbload(infofile,infovar);

fidlog   = fopen([outpath,moor,'_seaphox_stage1.log'],'a');
fprintf(fidlog,'Transformation of SeapHOx data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));
fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n\n\n',lon);

% bg = (julian([s_d(:)' hms2h([s_t(:)' 0])])); %start
% ed = julian([e_d(:)' hms2h([e_t(:)' 0])]); %end
bgvec = ([s_d(:)' ([s_t(:)' 0])]);
bg = datenum(bgvec); %start
edvec = ([e_d(:)' ([e_t(:)' 0])]); %end
ed = datenum(edvec);

vec=find(id==375); % SeapHOx instrument code

serial_nums=sn(vec);

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
    outfile=[moor '_' sprintf('%0d',input_ser_nums(i)) '.raw'];
    indep=find(sn==input_ser_nums(i));
    indep=z(indep);
    if length(indep)>1
        indep=indep(1);
    end
    fprintf(fidlog,'infile : %s\n',infile);
    fprintf(fidlog,'outfile: %s\n',outfile);
    fprintf(fidlog,'SeapHOx serial number  : %d\n',input_ser_nums(i));

    
    fid=fopen(infile,'r');
    
    % for older instrument used on DY120
    if serial_nums==4;
        header=textscan(fid,[repmat('%s ',1,12)  '%s'],1,'delimiter',',');
        all_data=textscan(fid,['%s %s %f' repmat(' %f',1,10)],'delimiter',',');
        date=datenum(all_data{2},'mm/dd/yyyy HH:MM:SS');
        [year,~,~,hour,MI,S]=datevec(date); 
        if strcmp(cruise,'dy120')==1
        doy=floor(date-datenum(year,1,1)) + 1; % leap year??
        else
            doy=floor(date-datenum(year,1,1));
        end
        
ph_ext=all_data{6}; % pH value for Deep SeaFETs is pHExt as there is no internal pH sensor as per shallow SeaFETs
            % The data should be converted from the raw voltages to
            % a pH value, but this was not happening sometimes on
            % DY039 with the reason not clear from discussions with
            % the manufacturer. However the raw sensor voltage is
            % recorded so this is used to calculate the pH.
mc_t=all_data{5};
mc_s=all_data{9};
mc_o_mlL=all_data{11}; 
mc_p=all_data{8};
v_ext=all_data{7}; % raw sensor voltage from external (ie. the only one for Deep Seafet) pH sensor 

% conversion from ml/L to micromol/kg
% Molar Volume of O2 at STP = 22.391 L
% 1 micromol O2 = 0.022391 mL
% 1 ml/L = 44.661 micromol/L
% 1 micromol/kg = 1000/(sw_dens(s,t,p) micromol/L
mc_o = mc_o_mlL.*44.661./(sw_dens(mc_s,mc_t,mc_p)./1000);

    %    [DD MM YY]=DOY2date(doy,year); 
        for idd=1:length(doy)
            [DD(idd,1) MM(idd,1) YY(idd,1)]=DOY2date(doy(idd),year(idd)); % created by LOH on DY078        
        end
    
    elseif serial_nums==117
        for k=1:8
            eval(['header' num2str(k) '=textscan(fid,''%s %s %s'',1,''delimiter'','','');'])
        end
all_data=textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %f','delimiter',',');
date=all_data{2};
year=floor(date/1000); doy=date-year*1000; 
hour=floor(all_data{3}); MI=floor((all_data{3}-hour)*60);
S=round(((all_data{3}-hour)*60)-MI)*60;
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
    
%    [DD MM YY]=DOY2date(doy,year); 
    for idd=1:length(doy)
        [DD(idd,1) MM(idd,1) YY(idd,1)]=DOY2date(doy(idd),year(idd)); % created by LOH on DY078        
    end

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
    
    end
     
if find(~isnan(ph_ext)) % implies pH was being output by the instrument when sometimes it is now - reason not clarified by SeaBird/Satlantic yet 
    response_valid=0;
    j=1;
    while response_valid~=1;
        if j>1;
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

%%%%%%%%%%%%%%%%%%%%% date and time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dat       = [YY MM DD hour];
jd        = datenum(YY,MM,DD,hour,MI,S);

% make sure the datetime array is in order and sort all daata to that order
[jd,ix]  =sort(jd,'ascend');
if serial_nums==117
    dat     =dat(ix,:);
    ph      =ph(ix);
    v_ext   =v_ext(ix);
    mc_t    =mc_t(ix);
    mc_p    =mc_p(ix);
    mc_o    =mc_o(ix);
    mc_s    =mc_s(ix);
    main_v  =main_v(ix);
    isol_v  =isol_v(ix);
else
    dat     =dat(ix,:);
    ph      =ph(ix);
    v_ext   =v_ext(ix);
    mc_t    =mc_t(ix);
    mc_p    =mc_p(ix);
    mc_o    =mc_o(ix);
    mc_s    =mc_s(ix);
end

valI      = find(jd<ed & jd>bg);

 
% ----- save data to rodb -----------------
if serial_nums ==4
    if find(intersect(input_ser_nums(i),serial_nums))
        columns = 'YY:MM:DD:HH:PH:PH_V:T:P:O:S:BATT1:BATT2';
        data = [dat ph v_ext mc_t mc_p mc_o mc_s ]; 
        fort =['%4.4d %2.2d %2.2d  %6.4f  %7.5f %10.8f  %6.4f %7.3f %6.3f  %6.4f'];
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

elseif serial_nums==117
    if find(intersect(input_ser_nums(i),serial_nums))
        columns = 'YY:MM:DD:HH:PH:PH_V:T:P:O:S:BATT1:BATT2';
        data = [dat ph v_ext mc_t mc_p mc_o mc_s main_v isol_v]; 
        fort =['%4.4d %2.2d %2.2d  %6.4f  %7.5f %10.8f  %6.4f %7.3f %6.3f  %6.4f %5.3f %5.3f'];
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
else
    disp('INSTRUMENT SERIAL NUMBER NEEDS UPDATES IN CODE HERE!!')
    return
end



    sz   =   size(data);
    
    % -------- generate logfile entries --------------
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
    
    fclose(fid);
    
    % -------- generate plot  --------------

    jd = datenum(data(:,1),data(:,2),data(:,3),data(:,4),MI,S);
    dl = length(jd);
    nplots = 5; %sz(2) - 4; % variables other than the four time vars

    ccplot=0;
    for k = [1 3 4 6 5]
        ccplot=ccplot+1;
        subplot(nplots,1,ccplot);
        if k==1
              title(['SeapHOx:  ',num2str(sn(vec)),'   Depth:  ',num2str(z(vec)),' m']) 
        end
        hold off;
        plot(jd(valI),data(valI,4+k));

        set(gca,'Ygrid','on')
        set(gca,'Xgrid','on')
        
        mdata=nanmedian(data(valI,4+k));
        if k == 1;
            ylabel('pH')
            fprintf(fidlog,'Median pH %5.2f\n',mdata);        
        elseif k == 3;
            ylabel('Temp.')
            fprintf(fidlog,'Median Temperature %5.2f\n',mdata);
        elseif k == 4;
            ylabel('Press.')
            fprintf(fidlog,'Median Pressure %5.2f\n',mdata);
         elseif k == 6;
            ylabel('Sal.')
            fprintf(fidlog,'Median Salinity %5.2f\n',mdata);           
        elseif k == 5;
            ylabel('Oxygen')
            fprintf(fidlog,'Median Oxygen %5.2f\n',mdata);
        end
    %     val = find(abs(data(:,4+k))<40);
    %     sdt = std(dt(val(100:end-100),1));
%         set(gca,'ylim',prctile(data(valI,4+k),[1 99])+[-1 1]*std(data(valI,4+k))/10) % SJ
%        set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
% set ylimits to be 3* the std deviation of the data
        if k==1
            set(gca,'ylim',[mdata-(1/2*(nanstd(data(valI,4+k)))) mdata+(1/2*(nanstd(data(valI,4+k))))])
        else
            set(gca,'ylim',[mdata-(3*(nanstd(data(valI,4+k)))) mdata+(3*(nanstd(data(valI,4+k))))])
        end
            title(['SeapHOx:  ',num2str(sn(vec)),'   Depth:  ',num2str(z(vec)),' m']) 

        set(gca,'xlim',[jd(valI(1)) jd(valI(end)) ])
        datetick('x','mmm-yy','Keeplimits');
        yl = get(gca,'ylim');
        xl = get(gca,'xlim');
        tx =text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['median: ',sprintf('%3.2f',mdata)]) ;
        set(tx,'FontWeight','bold','FontSize',11)
    end;

      orient tall
%       eval(['print -dpsc ',[outpath outfile],'.ps'])
      print(gcf, [outpath moor '_'  num2str(serial_nums(i))],'-dpng', '-r300');
%       eval(['print -dpsc ',[outpath outfile],'.ps'])
    
    
end  % loop of serial numbers
fclose(fidlog);

% TODO: Add CTD data on plots with all SeapHOx (if there is several on a cast) 
end