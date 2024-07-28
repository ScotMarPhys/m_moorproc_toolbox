function O2 = microcat_odo_recalculate_oxy(infile_raw)

% Function to recalculate SBE37-ODO oxygen concentration using SBE37
% temperature instead of SBE63. Gnereally not required, unless in case of
% failure of the SBE63 temperature sensor.

% Formula from SBE63 manual (rev 011):
% O2 (ml/L) = [((a0 + a1T + a2V^2) / (b0 + b1V) – 1) / Ksv] [SCorr] [PCorr]
% where:
%   O2 is oxygen concentration (ml/L)
%   T is temperature output from SBE 63’s thermistor in °C
%   V is raw measured phase delay in volts
%   a0, a1, a2, b0, b1 are calibration coefficients (Uchida et al, 2008)
%   Ksv is Stern-Volmer constant equation (Demas et al, 1999)
%       Ksv = c0 + c1T + c2T^2
%       where:
%           c0, c1, c2 are calibration coefficients
%           T is temperature output from SBE 63’s thermistor in °C
%   SCorr is salinity correction function (with calibration coefficients SolB0, SolB1, SolB2, SolB3, SolC0)
%       SCorr = exp [S * (SolB0 + SolB1 * Ts + SolB2 * Ts^2 + SolB3 * Ts^3) + SolC0 * S^2]
%       where:
%           SolBx = salinity correction coefficients (constants, Benson and Krause, 1984)
%               SolB0 = -6.24523e-3
%               SolB1 = -7.37614e-3
%               SolB2 = -1.03410e-2
%               SolB3 = -8.17083e-3
%               SolC0 = -4.88682e-7
%           Ts = ln [(298.15 – T) / (273.15 + T)]
%           S = salinity
%   PCorr is pressure correction function (with calibration coefficient E)
%       Pcorr = exp (E * P / K)
%       where:
%           E = pressure correction coefficient = 0.011
%           K = temperature in Kelvin = T + 273.15
%           P = pressure


%% Load ODO data

% Based on microcat2rodb.
% Can't use that function directly as 1) don't want to overwite the original
% rodb files; 2) don't want to change the output rodb format in case it
% crashes otherscripts somewhere else in the toolbo; and 3) need other
% variables from the header (SBE63 calibration coefficients)
% For now just use a copied version of the "read" part of microcat2rodb


%--------------------------------------------------------
% open input file, read data into string
%--------------------------------------------------------

SerialNumber=[]; Start_Time=[]; Start_Date=[];

fid1 = fopen(infile_raw,'r');
if fid1 == -1
    disp(['unable to open infile:  ',infile_raw])
    return
else
    disp(['loading ',infile_raw]);
end
zeile = fscanf(fid1,'%c');  %read data into string
fclose(fid1);

ret = sprintf('\n');
retx = findstr(zeile,ret);% car. return indices


%---------------------------------------------------------
% find serial number inheader
%--------------------------------------------------------
a=findstr(zeile,'Temperature SN'); % find serial number for .cnv files
if ~isempty(a)
    SerialNumber = str2num(zeile(a+17:a+22));
else
    disp('unable to find serial number\n\n')
    input('please enter valid serial number')
    SerialNumber = str2num(input('please enter valid serial number'));
end


%-------------------------------------------
% get data
%-------------------------------------------
% There is a bug here when runnning this code for older instrumnets (e.g.
% instruments with early serial numbers - 4608). The code below reads the
% number of lines in the .cnv file and estimaetes a start postion. The
% issue is that older microcats have significantly less header lines so the
% start position is well below the threshold number.I have added an elseif
% for just now as a sticking plaster - LAD 06/10/2020

% detect data column length
if length(retx)<700 %in case of short record from test (added by Loic H on DY078)
    if SerialNumber==4608
        [XXX,data_length] = max(hist(diff(retx(73:end)),1:300));% length of data columns
    else
        [XXX,data_length] = max(hist(diff(retx(304:end)),1:300));% length of data columns
    end
else
    [XXX,data_length] = max(hist(diff(retx),1:300));% length of data columns
end
% DR increased final term to 200 from 100 due to longer header in .cnv
% files   % changed to 300 PW

ii0 = find(diff(retx) == data_length); % data row index
ii1 = find(diff(ii0)>1);  %

bg = 2;
[XXX,sede] = find(abs(diff(retx(ii0(bg):ii0(length(ii0))))-data_length) >3);

if ~isempty(XXX)
    bg = bg+1;
    [XXX,sede] = find(abs(diff(retx(ii0(bg):ii0(length(ii0))))-data_length) >3);
end

if ~isempty(XXX)
    bg = bg +1;
    [XXX,sede] = find(abs(diff(retx(ii0(4):ii0(length(ii0))))-data_length) >3);
end


if ~isempty(XXX)
    disp('severe deviation from input format')
    % fprintf(fidlog,'conversion stopped - severe deviation from input format: %d\n',sede+ii0(1));
    msgbox(['MC',sprintf('%4.4d','SerialNumber'),'  ',sprintf('severe deviation from input format: %d\n',sede+ii0(1))],'conversion stopped')
    %    return
end

% check if there are deviations in data row length
data_begin = retx(ii0(bg))+1;
data_end   = retx(ii0(length(ii0)));
disp('data begin detected')

if ~isempty(ii1)
    disp('warning: deviation from format')
    %fprintf(fidlog,'warning: deviation from format \n');
end

% define data stream
dt = zeile(data_begin:data_end);
dret = find(dt == ret);
comx = findstr(dt(1:dret(1)),',');
retn = length(dret);% number of data columns

ii = findstr(dt,','); % replace comma by space
dt(ii) = ' ';
ii = findstr(dt,':'); % replace colon by space
dt(ii) = ' ';

% if cnv==0 %microcat_month2 only valid for .asc format files
%     dt = microcat_month2(dt);
% end

whos dt;
ii = findstr(dt,' .'); % temp. check for missing number
if ~isempty(ii)
    disp('insert dummies for missing elements')
    disp(['missing element element ',num2str(ii+1)])
    %fprintf(fidlog,'warning: missing element %s \n',num2str(ii+1));
    for i = 1 : length(ii);
        dt(ii(i):ii(i)+1) = '99';
    end
end

% --- convert data string
disp('Convert datastr into numbers')
%save tmp.mat dt ;

% This next line won't work if there're still ??? (carriage returns are ok) in dt
dt = str2num(dt);
sz = size(dt);

% catch and locate errors in raw microcat files where have spurious data
% scans for whatever reason. Typically these come at the end of the record,
% but on D382 we found some bad data in the middle of a record that was
% difficult to track down
if isempty(dt);
    load tmp.mat % reload dt data
    k2=1;
    for k=1:length(dt)/57
        string=dt(57*(k-1)+1:57*(k-1)+55);
        dt_convert=str2num(string);
        if isempty(dt_convert)
            bad_data_line(k2)=k;
            k2=k2+1;
        end
    end

    msgbox(['MC',sprintf('%4.4d',SerialNumber),': format deviation at sample ' num2str(bad_data_line(1))],'conversion stopped');
    %fprintf(fidlog,'conversion stopped - deviation from element number \n');

    return
end

% ---- check order of variables if cnv file ----

% --- short stub inserted by gdm, dy039 to ensure the indices match the cnv file ---
% --- aditional conditions added by bim dy039, to enable DR PC conversion files to be used.

namei = findstr(zeile,'name');
spani = findstr(zeile,'span');
toti = [namei spani(1)];

for k = 1:length(toti)-1
    varstr=zeile(toti(k):toti(k+1)-5);

    if ~isempty( strfind(varstr(10), 't') ) & ~isempty( strfind(varstr, '90') )
        tempi = k;
        continue;
    elseif ~isempty(strfind(varstr(10), 'c')) & ~isempty(strfind(varstr, '0S/m'))
        condi = k;
        continue;
    elseif ~isempty(strfind(varstr, 'pr')) |  ~isempty(strfind(varstr, 'prdM'))
        presi = k;
        continue;
    elseif ~isempty(strfind(varstr, 'sbeoxTC'))
        oxyti = k;
        continue;
    elseif ~isempty(strfind(varstr, 'sbeopoxMm/Kg'))
        oxyi = k;
        continue;
    elseif ~isempty(strfind(varstr, 'timeK'))|  ~isempty(strfind(varstr, 'timeS'))
        timeformat = 'sec';
        timei = k;
        continue;
    elseif ~isempty(strfind(varstr, 'timeJ'))
        timeformat = 'JD';
        timei = k;
        continue;
    elseif ~isempty(strfind(varstr(10), 'c')) & ~isempty(strfind(varstr, '0mS/cm'))
        condi = k;
        dt(:,condi)=dt(:,condi)/10;
        continue;
    elseif ~isempty(strfind(varstr, 'sbeopoxML/L'))
        oxymlli = k;
        continue;
    elseif ~isempty(strfind(varstr, 'sbeoxpdv'))
        oxyvi = k;
        continue;

    else
    end

end

if ~isempty(strfind(timeformat, 'sec'))
    if size(dt,2)==5 % SMP with pressure
        cnv_secs=dt(:,timei); % seconds since 1st Jan 2000.
        jd=cnv_secs/(60*60*24)+julian(2000,1,1,0) - toffset;
        disp(['toffset : ',num2str(toffset)]);
        gtime=gregorian(jd);
        HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6));
        Start_Date = gtime(1,[1 2 3]);
        End_Date = gtime(dtl,[1 2 3]);
    elseif size(dt,2)==7 % SMP-ODO with pressure
        cnv_secs=dt(:,timei); % seconds since 1st Jan 2000.
        jd=cnv_secs/(60*60*24)+julian(2000,1,1,0) - toffset;
        disp(['toffset : ',num2str(toffset)]);
        gtime=gregorian(jd);
        HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6));
        Start_Date = gtime(1,[1 2 3]);
        End_Date = gtime(dtl,[1 2 3]);
    end
elseif ~isempty(strfind(timeformat, 'JD'))
    if size(dt,2)==5 % SMP with pressure
        cnv_days=dt(:,4); % seconds since 1st Jan 2000.
        jd=cnv_days+julian(toffset,1,1,0) - 1;
        disp(['dateoffset : ',num2str(toffset)]);
        gtime=gregorian(jd);
        HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6));
        Start_Date = gtime(1,[1 2 3]);
        End_Date = gtime(dtl,[1 2 3]);
    elseif size(dt,2)==7 % SMP-ODO with pressure
        cnv_days=dt(:,4); % seconds since 1st Jan 2000.
        jd=cnv_days+julian(toffset,1,1,0) - 1;
        disp(['dateoffset : ',num2str(toffset)]);
        gtime=gregorian(jd);
        HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6));
        Start_Date = gtime(1,[1 2 3]);
        End_Date = gtime(dtl,[1 2 3]);                
    end      
end


% Get vars for oxygen calculations
time_s = dt(:,timei);
mtime = time_s/3600/24 + datenum('01-Jan-2000 00:00:00');
P = dt(:,presi);
T = dt(:,tempi);
C = dt(:,condi); if nanmedian(C)<10; C=C*10; end; % convert from S/m to mS/cm
S = sw_salt((C/sw_c3515),T,P);

if ~exist('oxyvi','var')

    % No raw sbe63 sensor voltage supplied, need it to recaculate oxygen
    % concentration
    disp('No raw SBE63 phase delay voltage found in data file, please re-export the data')
    disp('Oxygen concentration not recalculated')

else
    %% Recalculate oxygen concentration
    V_sbe63 = dt(:,oxyvi);

    %-----------------------------------------------
    % get SBE63 calibration coefficients from header
    %-----------------------------------------------
    oxcoeffi = findstr(zeile,'<OXA0>'); oxcoeffii = findstr(zeile,'</OXA0>'); oxA0 = str2double(zeile(oxcoeffi+6:oxcoeffii-1)); clear oxcoeff*;
    oxcoeffi = findstr(zeile,'<OXA1>'); oxcoeffii = findstr(zeile,'</OXA1>'); oxA1 = str2double(zeile(oxcoeffi+6:oxcoeffii-1)); clear oxcoeff*;
    oxcoeffi = findstr(zeile,'<OXA2>'); oxcoeffii = findstr(zeile,'</OXA2>'); oxA2 = str2double(zeile(oxcoeffi+6:oxcoeffii-1)); clear oxcoeff*;
    oxcoeffi = findstr(zeile,'<OXB0>'); oxcoeffii = findstr(zeile,'</OXB0>'); oxB0 = str2double(zeile(oxcoeffi+6:oxcoeffii-1)); clear oxcoeff*;
    oxcoeffi = findstr(zeile,'<OXB1>'); oxcoeffii = findstr(zeile,'</OXB1>'); oxB1 = str2double(zeile(oxcoeffi+6:oxcoeffii-1)); clear oxcoeff*;
    oxcoeffi = findstr(zeile,'<OXC0>'); oxcoeffii = findstr(zeile,'</OXC0>'); oxC0 = str2double(zeile(oxcoeffi+6:oxcoeffii-1)); clear oxcoeff*;
    oxcoeffi = findstr(zeile,'<OXC1>'); oxcoeffii = findstr(zeile,'</OXC1>'); oxC1 = str2double(zeile(oxcoeffi+6:oxcoeffii-1)); clear oxcoeff*;
    oxcoeffi = findstr(zeile,'<OXC2>'); oxcoeffii = findstr(zeile,'</OXC2>'); oxC2 = str2double(zeile(oxcoeffi+6:oxcoeffii-1)); clear oxcoeff*;

    % Recalculate oxygen

    % Calculate SCorr
    SolB0 = -6.24523e-3;
    SolB1 = -7.37614e-3;
    SolB2 = -1.03410e-2;
    SolB3 = -8.17083e-3;
    SolC0 = -4.88682e-7;
    Ts=nan(length(S),1);
    SCorr=nan(length(S),1);
    for j=1:length(S)
        Ts(j,1) = log((298.15 - T(j,1)) / (273.15 + T(j,1)));
        SCorr(j,1) = exp(S(j,1)* (SolB0 + Ts(j,1)*SolB1 + SolB2*Ts(j,1)^2 + SolB3*Ts(j,1)^3) + SolC0 * S(j,1)^2);
    end

    % Calculate Pcorr
    K = T + 273.15;
    PCorr(:,1) = exp(0.011*P./K.*1);

    % Calculate Ksv
    Ksv = oxC0 + T.*oxC1 + T.^2*oxC2;

    % Calculate oxygen concentration in ml/l
    for j=1:length(P)
        O_sbe37(j,1) = (((oxA0 + oxA1*T(j,1) + oxA2*V_sbe63(j,1)^2) / (oxB0 + oxB1*V_sbe63(j,1)) - 1) / Ksv(j,1)) * SCorr(j,1) * PCorr(j,1);
    end
    ind_bad = find(imag(O_sbe37)); O_sbe37(ind_bad) = NaN;


    %% Convert oxygen from ml/l to umol/kg
    %formula from seaphox2rodb_01.m)
    % Molar Volume of O2 at STP = 22.391 L
    % 1 micromol O2 = 0.022391 mL
    % 1 ml/L = 44.661 micromol/L
    % 1 micromol/kg = 1000/(sw_dens(s,t,p) micromol/L
    O2 = O_sbe37.*44.661./(sw_dens(S,T,P)./1000);

    disp('Oxygen concentration recalculated using sbe37 temperature')

end

