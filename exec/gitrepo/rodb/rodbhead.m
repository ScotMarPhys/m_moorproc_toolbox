function [vars,formats,lengths,names] = rodbhead(ident)
%RODBHEAD Decode RO database keywords.
%  [VARS,FORMATS,LENGTHS,NAMES] = RODBHEAD(IDENT)
%
%  input  :      ident           : header array of RODB file
%
%  output :      vars            : keyword numbers
%                formats         : keyword formats
%                lengths         : data lengths of keyword information
%
%  TYPE rodbhead  help   FOR MORE INFORMATIONS ON IDENTIFIERS

%  G.Krahmann, IfM Kiel, Oct 1995
%  number changed to EPIC conventions    G.Krahmann, Feb 1996, 0.1.0-->0.1.1
%  added identifier                      G.Krahmann, Feb 1996, 0.1.1-->0.1.2
%  added identifier and 'disp'           G.Krahmann, Feb 1996, 0.1.2-->0.1.3
%  about 50% faster                      C. Mertens, May 1997, 0.1.3-->0.1.4
%  version 0.1.4         last change 26.05.1997
%  added several header and data variables   d.kieke, Aug 2000
%  added several variables  T Kanzow, summer 2005
%  added O2 oxygen variable in umol/kg, D. Rayner 2014
%  added pH and pH voltage variables, L. Drysdale 2021
global allkeys allkeys1 nall


if isempty(allkeys)
%
% keyword array (first 20 characters) , its C-format (8 characters),
% the number of rows, which are necessary for the information and the
% keyword identifying number, negative values indicate header keywords
%
allkeys = [ ...
'Filename            %s       0   -1'; ...  % (raw data) file name
'Shipname            %s       0   -2'; ...  % ship name
'Cruise              %s       0   -3'; ...  % cruise identifier
'Mooring             %s       0   -4'; ...  % mooring identifer
'Model               %s       0   -5'; ...  % model identifer
'OldMooringName      %S       0   -6'; ...  % Old mooring name to cover lander name changes
'Instrument          %s       0  -10'; ...  % instrument identifier
'InstrType           %s       0  -10'; ...  % instrument identifier
'Serial_Number       %d       1  -11'; ...  % serial number
'SerialNumber        %d       1  -11'; ...  % serial number
'Station             %s       0  -12'; ...  % station number or name
'Profile             %s       0  -13'; ...  % profile identifier
'ArgonautSN          %d       1  -15'; ...  % Argonaut Serial number   
'MicrocatSN          %d       1  -16'; ...  % MicroCat Serial number
'CTD_File            %s       0  -20'; ...  % CTD file name
'Meteo_File          %s       0  -21'; ...  % meterological data file name
'WaterTemp           %g       1  -31'; ...  % water temp (deg C)
'WaterDepth          %d       1  -32'; ...  % water depth (m)
'Max_Depth           %d       1  -33'; ...  % maximum profile depth (m)
'MaxDepth            %d       1  -34'; ...  % maximum profile depth (m)
'Max_Press           %d       1  -41'; ...  % maximum profile pressure (dbar)
'MaxPress            %d       1  -41'; ...  % maximum profile pressure (dbar)
'InstrDepth          %d       1  -35'; ...  % instrument depth (m)
'Depth               %d       1  -36'; ...  % (instrument) depth (m)
'SensorSpacing       %d       1  -37'; ...  % sensor spacing (m) (thermistor strings)
'LayerThickness      %g       1  -38'; ...  % 
'Press_Offset        %g       1  -39'; ...  % pressure offset (dbar)
'PressOffset         %g       1  -39'; ...  % pressure offset (dbar)
'Mag_Deviation       %g       1  -40'; ...
'MagDeviation        %g       1  -40'; ...
'LatentHeat          %g       1  -50'; ...
'ThermalExpansion    %g       1  -51'; ...
'HalineContraction   %g       1  -52'; ...
'ReferenceTemperature%g       1  -53'; ...
'ReferenceSalinity   %g       1  -54'; ...
'ReferenceDensity    %g       1  -55'; ...
'Temp_Coef           %s       0  -60'; ...
'Cond_Coef           %s       0  -61'; ...
'Gravity             %g       1  -70'; ...
'Start_Date          %d/%d/%d 3  -80'; ...  % start date
'StartDate           %d/%d/%d 3  -80'; ...  % start date
'Date                %d/%d/%d 3  -80'; ...  % (start) date
'Start_Time          %d:%d    2  -81'; ...  % start time
'StartTime           %d:%d    2  -81'; ...  % start time
'Time                %d:%d    2  -81'; ...  % (start) time
'End_Date            %d/%d/%d 3  -82'; ...  % end date
'EndDate             %d/%d/%d 3  -82'; ...  % end date
'End_Time            %d:%d    2  -83'; ...  % end time
'EndTime             %d:%d    2  -83'; ...  % end time
'Time_Step           %g       1  -84'; ...
'TimeStep            %g       1  -84'; ...
'SamplingInterval    %g       1  -84'; ...
'Start_Lat           pos      1  -90'; ...  % start latitude
'StartLat            pos      1  -90'; ...  % start latitude
'Latitude            pos      1  -90'; ...  % (start) latitude
'Start_Lon           pos      1  -91'; ...  % start longitude
'StartLon            pos      1  -91'; ...  % start longitude
'Start_Long          pos      1  -91'; ...  % start longitude
'StartLong           pos      1  -91'; ...  % start longitude
'Longitude           pos      1  -91'; ...  % (start) longitude
'End_Lat             pos      1  -92'; ...  % end latitude
'EndLat              pos      1  -92'; ...  % end latitude
'End_Lon             pos      1  -93'; ...  % end longitude
'EndLon              pos      1  -93'; ...  % end longitude
'End_Long            pos      1  -93'; ...  % end longitude
'EndLong             pos      1  -93'; ...  % end longitude
'deployment          %d       1  -95'; ...  % number of deployment
'Columns             col      0 -100'; ...  % Columns in data file
'Columns cont.       col      0 -101'; ...  % Columns in data file
'ColumnUnits         %s       0 -102'; ...  % Physical Units of columns
'P                   dat      0    1'; ...  % pressure (dbar)
'D                   dat      0    2'; ...  % water depth (m)
'Z                   dat      0    3'; ...  % depth (m)
'PARG                dat      0    5'; ...  % press. of sontek argonaut (dbar)
'PCAT                dat      0    6'; ...  % press. of seabird MicroCAT (dbar)
'MLD                 dat      0    8'; ...  % mixed layer depth (m)
'T                   dat      0   20'; ...  % temperature (deg C)
'WT                  dat      0   19'; ...  % wet bulb temperature (deg C)
'AT                  dat      0   21'; ...  % air temperature (deg C)
'SST                 dat      0   25'; ...  % sea surface temperature (deg C)
'PT                  dat      0   30'; ...  % potential temp (deg C)
'TARG                dat      0   33'; ...  % temp. of sontek argonaut (deg C)
'TCAT                dat      0   34'; ...  % temp. of seabird MicroCAT (deg C)
'TC                  dat      0   35'; ... % temperatre taken from C sensor (deg C) (Seaguard)
'TP                  dat      0   36'; ... % temperatre taken from P sensor (deg C) (Seaguard)
'S                   dat      0   41'; ...  % salinity
'C                   dat      0   50'; ...  % conductivity (mS/cm)
'O                   dat      0   60'; ...  % oxygen (ml/l)
'BO                  dat      0   61'; ...  % bottle oxygen (ml/l)
'O2                  dat      0   62'; ...  % oxygen (umol/kg)
'DO                  dat      0   67'; ...  % oxygen difference (ml/l)
'PH_V                dat      0   68'; ...  % pH
'PH                  dat      0   69'; ...  % pH Voltage
'ST                  dat      0   70'; ...  % sigma-t (kg/m^3)
'STH                 dat      0   71'; ...  % sigma-theta (kg/m^3)
'SV                  dat      0   80'; ...  % sound velocity (m/s)
'BATT1               dat      0   81'; ...  % Main battery voltage for seaphox 
'BATT2               dat      0   82'; ...  % 2nd battery voltage for seaphox
'TT                  dat      0   85'; ...  % acoustic travel time (s)
'TT1                 dat      0   86'; ...  % acoustic travel time 1 (PIES) (s)
'TT2                 dat      0   87'; ...  % acoustic travel time 2 (PIES) (s)
'TT3                 dat      0   88'; ...  % acoustic travel time 3 (PIES) (s)
'TT4                 dat      0   89'; ...  % acoustic travel time 4 (PIES) (s)
'PFIT                dat      0   90'; ...  % empirical drift estimate for pressure record
'PTIDE               dat      0   92'; ...  % tidal estimate for pressure record
'PTIDEM              dat      0   93'; ...  % monthly and semi-monthly tidal estimate for pressure record
'RCMC1               dat      0   95'; ...  % lower limit of RCM11 conductivity range [mS/cm]
'RCMC2               dat      0   96'; ...  % upper limit of RCM11 conductivity range [mS/cm]
'REF                 dat      0   97'; ...  % reference number fo Aanderra instruments
'MSS                 dat      0   98'; ...  % Mean Signal strength
'SN                  dat      0  100'; ...  % scan number (or serian number)
'SCN                 dat      0  100'; ...  % scan number
'BTL                 dat      0  103'; ...  % niskin bottle number
'PRF                 dat      0  104'; ...  % profile/cast number
'STN                 dat      0  105'; ...  % station
'ID                  dat      0  105'; ...  % identifier
'OC                  dat      0  110'; ...  % oxygen current (micro A)
'OT                  dat      0  111'; ...  % oxygen temperature (deg C)
'DOC                 dat      0  112'; ...  % dOc/Dt (micro A/s)
'QS                  dat      0  133'; ...  % shortwave radiation (W/m^2)
'QL                  dat      0  136'; ...  % longwave radiation (W/m^2)
'QH                  dat      0  137'; ...  % latent heat flux (W/m^2)
'QB                  dat      0  138'; ...  % sensible heat flux (W/m^2)
'PR                  dat      0  139'; ...  % precipitation
'F11                 dat      0  150'; ...  % freon 11 (micro moles/kg)
'F12                 dat      0  151'; ...  % freon 12 (micro moles/kg)
'NO3                 dat      0  182'; ...  % nitrate (micro moles/l)
'NO2                 dat      0  184'; ...  % nitrite (micro moles/l)
'NN                  dat      0  185'; ...  % nitrate + nitrite (micro moles/l)
'PO4                 dat      0  186'; ...  % phosphate (micro moles/l)
'SI                  dat      0  188'; ...  % silicate (micro moles/l)
'CCL4                dat      0  190'; ...  % ccl4, carbon tetrachloride (pico moles/kg)
'MTHCF               dat      0  191'; ...  % methylchloroform (pico moles/kg)
'IC                  dat      0  215'; ...  % ice concentration (percent)
'CS                  dat      0  300'; ...  % current speed (cm/s)
'CD                  dat      0  310'; ...  % current direction (deg)
'U                   dat      0  320'; ...  % zonal current (cm/s)
'V                   dat      0  321'; ...  % meridional current (cm/s)
'W                   dat      0  329'; ...  % vertical velocity (cm/s)
'X                   dat      0  370'; ...  % x position (m)
'Y                   dat      0  371'; ...  % y position (m)
'L                   dat      0  375'; ...  % length, length scale (m)
'WS                  dat      0  401'; ...  % wind speed (m/s)
'WD                  dat      0  410'; ...  % wind direction (deg) (0-360)
'WU                  dat      0  422'; ...  % wind speed (m/s) (going E)
'WV                  dat      0  423'; ...  % wind speed (m/s) (going N)
'TX                  dat      0  446'; ...  % zonal wind stress (N/m^2)
'TY                  dat      0  447'; ...  % meridional wind stress (N/m^2)
'LAT                 dat      0  500'; ...  % latitude (degrees north)
'LON                 dat      0  501'; ...  % longitude (degrees east)
'YY                  dat      0  601'; ...  % year
'YYYY                dat      0  601'; ...  % year
'MM                  dat      0  602'; ...  % month
'DD                  dat      0  603'; ...  % day
'HH                  dat      0  604'; ...  % hour (decimal)
'MIN                 dat      0  605'; ...  % Minute
'SEC                 dat      0  606'; ...  % Second 
'TIM                 dat      0  626'; ...  % time since start (s)
'SD                  dat      0  700'; ...  % standard deviation
'SK                  dat      0  701'; ...  % skewness
'PG                  dat      0  721'; ...  % ADCP percent good pings
'EV                  dat      0  722'; ...  % ADCP error velocity (cm/s)
'AMP                 dat      0  723'; ...  % ADCP amplitude
'SH                  dat      0  913'; ...  % Specific humidity (hPa)
'AGC                 dat      0 1202'; ...  % ADCP Echo intensity (counts)
'HDG                 dat      0 1215'; ...  % ADCP heading (deg)
'PIT                 dat      0 1216'; ...  % ADCP pitch (deg)
'ROL                 dat      0 1217'; ...  % ADCP roll (deg)
'TLT                 dat      0 1218'; ...  % FSI ACM Tilt (deg)
'F113                dat      0 1703'; ...  % freon 113 (pico moles/kg)
'TCO2                dat      0 1751'; ...  % total carbon dioxide (micro moles)
'DIC                 dat      0 1753'; ...  % Dissolved Inorganic Carbon (micromoles/kg)
'ALK                 dat      0 1756'; ...  % alkalinity (micro moles/kg)
'USD                 dat      0 2001'; ...  % stand. dev. zonal current (cm/s)
'VSD                 dat      0 2002'; ...  % stand. dev. merid current (cm/s)
'WSD                 dat      0 2003'; ...  % stand. dev. vert velocity (cm/s)
'CSSD                dat      0 2004'; ...  % stand. dv. current speed (cm/s)
'USS                 dat      0 2011'; ...  % signal strength zon. cur. (counts)
'VSS                 dat      0 2012'; ...  % signal strength mer. cur (counts)
'WSS                 dat      0 2013'; ...  % signal strength ver. cur. (counts)
'UNOISE              dat      0 2014'; ...  % signal noise zon cur. (counts) - NB: Not strictly accurate. Should be BEAM1NOISE
'VNOISE              dat      0 2015'; ...  % signal noise mer cur. (counts) - NB: Not strictly accurate. Should be BEAM2NOISE
'WNOISE              dat      0 2016'; ...  % signal noise ver cur. (counts) - NB: Not strictly accurate. Should be BEAM3NOISE
'USNR                dat      0 2017'; ...  % signal to noise ratio zon cur - NB: Not strictly accurate. Should be BEAM1SNR
'VSNR                dat      0 2018'; ...  % signal to noise ratio mer cur - NB: Not strictly accurate. Should be BEAM2SNR
'WSNR                dat      0 2019'; ...  % signal to noise ratio ver cur - NB: Not strictly accurate. Should be BEAM3SNR
'PGP                 dat      0 2021'; ...  % percent good pings (percent)
'PITSD               dat      0 2022'; ...  % Pitch stand. dev. [deg] for Argonaut
'ROLSD               dat      0 2023'; ...  % Roll stand. dev. [deg] for Argonaut
'HDGSD               dat      0 2024'; ...  % Heading stand. dev [deg] for Argonaut
'PSD                 dat      0 2031'; ...  % Pressure stand. dev.[dbar] for Argonaut
'IPOW                dat      0 2041'; ...  % input power level  [V] for Argonaut
'CBEG                dat      0 2042'; ...  % start loc.of sampl. volume [m] for Argonaut 
'CEND                dat      0 2043'; ...  % end loc. of sampl. volume [m] for Argonaut
'IREL                dat      0 2044'; ...  % release battery current drain [Amperes] for PIES
'ISYS                dat      0 2045'; ...  % system battery current drain [Amperes] for PIES
'VREL                dat      0 2046'; ...  % release battery voltage [V] for PIES
'VSYS                dat      0 2047'; ...  % system battery current drain [V] for PIES
'BEAM1SS             dat      0 2050'; ...  % Signal Strength for Beam 1 (counts)
'BEAM2SS             dat      0 2051'; ...  % Signal Strength for Beam 2 (counts)
'BEAM3SS             dat      0 2052'; ...  % Signal Strength for Beam 3 (counts)
'BEAM1NOISE          dat      0 2053'; ...  % Signal noise for Beam 1 (counts)
'BEAM2NOISE          dat      0 2054'; ...  % Signal noise for Beam 2 (counts)
'BEAM3NOISE          dat      0 2055'; ...  % Signal noise for Beam 3 (counts)
'BEAM1SNR            dat      0 2056'; ...  % Signal to Noise ratio for Beam 1 (dB)
'BEAM2SNR            dat      0 2057'; ...  % Signal to Noise ratio for Beam 2 (dB)
'BEAM3SNR            dat      0 2058'; ...  % Signal to Noise ratio for Beam 3 (dB)
'BEAM4SS             dat      0 2059'; ...  % Signal Strength for Beam 4 (counts)
'BEAM4NOISE          dat      0 2060'; ...  % Signal noise for Beam 4 (counts)
'BEAM4SNR            dat      0 2061'; ...  % Signal to Noise ratio for Beam 4 (dB)
'BEAM1COR            dat      0 2062'; ...  % Correlation for Beam 1 (counts)
'BEAM2COR            dat      0 2063'; ...  % Correlation for Beam 2 (counts)
'BEAM3COR            dat      0 2064'; ...  % Correlation for Beam 3 (counts)
'BEAM4COR            dat      0 2065'; ...  % Correlation for Beam 4 (counts)
'BEAM1PGP            dat      0 2066'; ...  % percent good pings for Beam 1 (percent)
'BEAM2PGP            dat      0 2067'; ...  % percent good pings for Beam 2 (percent)
'BEAM3PGP            dat      0 2068'; ...  % percent good pings for Beam 3 (percent)
'BEAM4PGP            dat      0 2069'; ...  % percent good pings for Beam 4 (percent)
'T1                  dat      0 9111'; ...  % temp TC sensor 1 (deg C)
'T2                  dat      0 9112'; ...  % temp TC sensor 2 (deg C)
'T3                  dat      0 9113'; ...  % temp TC sensor 3 (deg C)
'T4                  dat      0 9114'; ...  % temp TC sensor 4 (deg C)
'T5                  dat      0 9115'; ...  % temp TC sensor 5 (deg C)
'T6                  dat      0 9116'; ...  % temp TC sensor 6 (deg C)
'T7                  dat      0 9117'; ...  % temp TC sensor 7 (deg C)
'T8                  dat      0 9118'; ...  % temp TC sensor 8 (deg C)
'T9                  dat      0 9119'; ...  % temp TC sensor 9 (deg C)
'T10                 dat      0 9120'; ...  % temp TC sensor 10 (deg C)
'T11                 dat      0 9121'; ...  % temp TC sensor 11 (deg C)
'T12                 dat      0 9122'; ...  % temp TC sensor 12 (deg C)
'PT1                 dat      0 9131'; ...  % pot temp TC sensor 1 (deg C)
'PT2                 dat      0 9132'; ...  % pot temp TC sensor 2 (deg C)
'PT3                 dat      0 9133'; ...  % pot temp TC sensor 3 (deg C)
'PT4                 dat      0 9134'; ...  % pot temp TC sensor 4 (deg C)
'PT5                 dat      0 9135'; ...  % pot temp TC sensor 5 (deg C)
'PT6                 dat      0 9136'; ...  % pot temp TC sensor 6 (deg C)
'PT7                 dat      0 9137'; ...  % pot temp TC sensor 7 (deg C)
'PT8                 dat      0 9138'; ...  % pot temp TC sensor 8 (deg C)
'PT9                 dat      0 9139'; ...  % pot temp TC sensor 9 (deg C)
'PT10                dat      0 9140'; ...  % pot temp TC sensor 10 (deg C)
'PT11                dat      0 9141'; ...  % pot temp TC sensor 11 (deg C)
'PT12                dat      0 9142'; ...  % pot temp TC sensor 12 (deg C)
'DPT                 dat      0 9143'; ...  % Depth calc. with Moordyn (m) 
'                                   ';];

allkeys1 = allkeys;
allkeys(:,1:20) = upper(allkeys(:,1:20));
nall = abs(allkeys(:,1:20));

end


% delimiter
d = ':,';


%
% check for multiple keywords in 'ident'
% or for header keyword list
%

s = size(allkeys);
ident = upper(ident);
nident = abs(ident);
sident = size(ident);
if sident(1) == 1
  multi = [];
  for i = 1:length(d)
    multi = [multi,find(nident == abs(d(i)))];
  end
else
  multi = [];
end
if isempty(multi)
  [dummy,xma] = min(abs(nident' - abs(' ')));
  [dummy2,xma2] = min(abs(nident' - abs('=')));
  dummy = dummy + dummy2;
  xma = min([xma;xma2]);
  if all(dummy)
    xmi = 1;
    xma = sident(2);
    yyy = 1;
  else
    xmi = ones(1,length(xma));
    xma = xma - 1;
    yyy = [1:sident(1)];
  end    
  nkeys = sident(1);
else
  dels = [];
  for i = 1:length(d)
    dels = [dels,find(ident == abs(d(i)))];
  end
  xmi = [1,dels+1];
  xma = [dels-1,length(ident)];
  yyy = ones(1,length(xmi));
  nkeys = length(dels) + 1;
end     

%
% search for keywords in keyword list
%
names = char(32*ones(nkeys,20));
formats = char(32*ones(nkeys,8));
vars = NaN*ones(nkeys,1);
lengths = NaN*ones(nkeys,1);
for n = 1:nkeys

  i = find(nall(:,1) == nident(yyy(n),xmi(n)));
  if length(i) > 1
    m = length(i);
    s1 = ones(m,1)* ...
          [nident(yyy(n),xmi(n):xma(n)),32*ones(1,20-(xma(n)-xmi(n)+1))];
    j = find(sum(abs(s1 - nall(i,:))') == 0);
    i = i(j);
  end

  if ~isempty(i)
    names(n,:) = allkeys1(i,1:20);
    vars(n) = sscanf(allkeys(i,31:35),'%d',1);
    formats(n,:) = allkeys(i,21:28);
    lengths(n) = abs(allkeys(i,30) - 48);
  else
    vars(n) = -9999;
    formats(n,:) = blanks(8);
    lengths(n) = -9999;
  end

end


