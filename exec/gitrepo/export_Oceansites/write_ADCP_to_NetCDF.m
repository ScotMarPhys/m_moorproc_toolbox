
% ------------------------------------------------------------------
% Copyright, IFREMER, 2012
% Modified by Thomas Richter at NIOZ, 2014
% Adapted for OSNAP ADCP data by Feili Li at Duke University, 2016
% 
% Last updated: August 12, 2016
% Edit by Loic Houpert, SAMS 18/11/2016
% ------------------------------------------------------------------


function write_ADCP_to_NetCDF(filename, moorname, latitude, longitude, adcpinfo, bin_depth, instrument_depth, pressure, time, vel_east, vel_north, vel_vert, vel_err)
%    dev2, corr, echo, perg, pitch_series, roll_series, heading_series, temp_series)

%
% ----------------------------------------------------------------- 
% NOTES ON INPUT
%
%   outfile-   output nc file name
%   moorname-  name of the mooring
%   latitude-   latitude of measurement  [degree north]
%   longitude-  longitude of measurement  [degree east]
%   adcpinfo-  structure with adcp metadata  
%   bin_depth-  depths of velocity bins [m]
%   instrument_depth-	nominal depth of measurement  [m]
%   time-       time of measurement
%   pressure-   sea water pressure directly measured  [db]
%   vel_east-   eastward sea water veloicty [m/s]
%   vel_north-  northward sea water veloicty [m/s]
%   vel_vert-	vertical sea water veloicty [m/s]
%   vel_err-	eorror in sea water velocity [m/s]
%
% The following parameters are needed when available:
%   dev2-           magnetic compass correction [degree]
%   corr-           correlation magnitude [counts]
%   echo-           echo amplitude [counts]
%   perg-           percentage of good resolution [percent]
%   pitch_series-   instrument pitch [degree]
%   roll_series-    instrument pitch [degree]
%   heading_series- pressure of each measurement [db]
%   temp_series-    sea water temperature [degree C]
%   
% -----------------------------------------------------------------




% --------------------+
%   Begin             |
% --------------------+

fprintf('Begin %s\n', datestr(now));

% Dimensions
TimeDim = length(time);
binnber = length(bin_depth);
% Open the file to write
outfile = [filename,'.nc'];
nc = netcdf.create(outfile,'CLOBBER');
netcdf.close(nc)
 

% ----------------------------+
%   Write global attributes   |
% ----------------------------+


% === DISCOVERY AND IDENTIFICATION ===
ncwriteatt(outfile,'/','site_code', 'OSNAP');
ncwriteatt(outfile,'/','platform_code', ['OSNAP-' moorname]);	% CHANGE
ncwriteatt(outfile,'/','data_mode', 'D');
ncwriteatt(outfile,'/','title', ['OSNAP ADCP data ',datestr(time(1),'mm/yyyy'),'-',datestr(time(end),'mm/yyyy')]);  % CHANGE
ncwriteatt(outfile,'/','summary', ['current records for depth range ',num2str(min(bin_depth)),'m-',num2str(max(bin_depth)),'m']);  %CHANGE
ncwriteatt(outfile,'/','naming_authority', 'OceanSITES'); 
ncwriteatt(outfile,'/','id', filename);   
ncwriteatt(outfile,'/','source', 'subsurface mooring');
ncwriteatt(outfile,'/','principal_investigator', adcpinfo.principal_investigator);	% CHANGE
ncwriteatt(outfile,'/','principal_investigator_email', adcpinfo.principal_investigator_email);	% CHANGE
ncwriteatt(outfile,'/','principal_investigator_url', adcpinfo.principal_investigator_url);	% CHANGE
ncwriteatt(outfile,'/','institution', adcpinfo.institution);	% CHANGE
ncwriteatt(outfile,'/','project', adcpinfo.project);
ncwriteatt(outfile,'/','array', adcpinfo.array);  
ncwriteatt(outfile,'/','network', adcpinfo.network);  
% ncwriteatt(outfile,'/','comment', 'N/A');	% CHANGE 



    
% === GEO-SPATIAL-TEMPORAL ===
ncwriteatt(outfile,'/','area', 'North Atlantic Ocean');
ncwriteatt(outfile,'/','geospatial_lat_min', sprintf('%.1f',latitude));
ncwriteatt(outfile,'/','geospatial_lat_max', sprintf('%.1f',latitude));
ncwriteatt(outfile,'/','geospatial_lat_units', 'degree_north');
ncwriteatt(outfile,'/','geospatial_lon_min', sprintf('%.1f',longitude));
ncwriteatt(outfile,'/','geospatial_lon_max', sprintf('%.1f',longitude));
ncwriteatt(outfile,'/','geospatial_lon_units', 'degree_east');
ncwriteatt(outfile,'/','geospatial_vertical_min', min(bin_depth));
ncwriteatt(outfile,'/','geospatial_vertical_max', max(bin_depth));
ncwriteatt(outfile,'/','geospatial_vertical_positive', 'down');
ncwriteatt(outfile,'/','geospatial_vertical_units', 'meter');
ncwriteatt(outfile,'/','time_coverage_start', datestr(min(time),'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','time_coverage_end', datestr(max(time),'yyyy-mm-ddTHH:MM:SSZ'));

dDD = floor(time(end) - time(1));
dHH = floor((time(end) - time(1) - dDD)*24);
ncwriteatt(outfile,'/','time_coverage_duration', ['P',num2str(dDD),'D',num2str(dHH),'H']);

dMM = round(mean(diff(time)*24*60));
ncwriteatt(outfile,'/','time_coverage_resolution', ['PT',num2str(dMM),'M']);    

ncwriteatt(outfile,'/','cdm_data_type', 'Station');
ncwriteatt(outfile,'/','featureType','timeSeries');
ncwriteatt(outfile,'/','data_type', 'OceanSITES time-series data');




% === NC CONVENTIONS ===
ncwriteatt(outfile,'/','format_version', '1.3');
ncwriteatt(outfile,'/','Conventions','CF-1.6,OceanSITES-1.3,ACDD-1.2')
ncwriteatt(outfile,'/','netcdf_version','4.3');




% === PUBLICATION INFORMATION ===
ncwriteatt(outfile,'/','publisher_name', 'Yao Fu');  %CHANGE     
ncwriteatt(outfile,'/','publisher_email', 'yaofu@gatech.edu');  %CHANGE   
ncwriteatt(outfile,'/','publisher_url', 'http://www.o-snap.org');  %CHANGE  
ncwriteatt(outfile,'/','references', 'http://www.o-snap.org, http://www.oceansites.org');
ncwriteatt(outfile,'/','data_assembly_center', 'Georgia Institute of Technology');  %CHANGE  
ncwriteatt(outfile,'/','update_interval', 'void');
ncwriteatt(outfile,'/','license', 'Follows CLIVAR (Climate Variability and Predictability) standards, cf. http://www.clivar.org/data/data_policy.php. Data available free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.');
ncwriteatt(outfile,'/','citation', 'These data were collected and made freely available by the OSNAP project and the national programs that contribute to it.');
ncwriteatt(outfile,'/','acknowledgement', 'Funding source: the UK Natural Environment Research Council (NERC), UK OSNAP project');	% CHANGE




% === PROVENANCE ===
ncwriteatt(outfile,'/','date_created', datestr(now+8/24,'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','date_modified', datestr(now+8/24,'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','history', 'Delayed time processed quality controlled');
ncwriteatt(outfile,'/','processing_level',adcpinfo.processing_level);
ncwriteatt(outfile,'/','QC_indicator', adcpinfo.QC_indicator);
ncwriteatt(outfile,'/','contributor_name', adcpinfo.contributor_name);	% CHANGE
ncwriteatt(outfile,'/','contributor_role', adcpinfo.contributor_role);	% CHANGE
ncwriteatt(outfile,'/','contributor_email', adcpinfo.contributor_email);	% CHANGE
  
  



% --------------------------------------+
%   Write coordinate variables          |
% --------------------------------------+
nccreate(outfile,'TIME', 'Dimensions',{'TIME',TimeDim}, 'Datatype','double');
ncwrite(outfile,'TIME', double(time) - datenum(1950,1,1,0,0,0));
ncwriteatt(outfile,'TIME', 'standard_name','time');
ncwriteatt(outfile,'TIME', 'units','days since 1950-01-01T00:00:00Z');
ncwriteatt(outfile,'TIME', 'axis','T');
ncwriteatt(outfile,'TIME', 'long_name','time of measurement');
ncwriteatt(outfile,'TIME', 'conventions', 'Relative julian days with decimal part (as parts of the day)');
ncwriteatt(outfile,'TIME', 'valid_min', single(0));
ncwriteatt(outfile,'TIME', 'valid_max', single(90000));
ncwriteatt(outfile,'TIME', 'QC_indicator','good data');
ncwriteatt(outfile,'TIME', 'processing_level','Data manually reviewed');
ncwriteatt(outfile,'TIME', 'comment','end points of measurement interval');	% CHANGE



nccreate(outfile,'LATITUDE', 'Dimensions',{'LATITUDE',1}, 'Datatype','single');
ncwrite(outfile,'LATITUDE', single(latitude));
ncwriteatt(outfile,'LATITUDE', 'standard_name','latitude');
ncwriteatt(outfile,'LATITUDE', 'units','degree_north');
ncwriteatt(outfile,'LATITUDE', 'axis','Y');
ncwriteatt(outfile,'LATITUDE', 'long_name','latitude of measurement');
ncwriteatt(outfile,'LATITUDE', 'reference','WGS84');
ncwriteatt(outfile,'LATITUDE', 'coordinate_reference_frame','urn:ogc:crs:EPSG::4326');
ncwriteatt(outfile,'LATITUDE', 'valid_min',-90);
ncwriteatt(outfile,'LATITUDE', 'valid_max',90);
ncwriteatt(outfile,'LATITUDE', 'QC_indicator','good data');
ncwriteatt(outfile,'LATITUDE', 'processing_level','Data manually reviewed');
ncwriteatt(outfile,'LATITUDE', 'uncertainty',.004);	% CHANGE
% ncwriteatt(outfile,'LATITUDE', 'comment','N/A');	% CHANGE



nccreate(outfile,'LONGITUDE', 'Dimensions',{'LONGITUDE',1}, 'Datatype','single');
ncwrite(outfile,'LONGITUDE', single(longitude));
ncwriteatt(outfile,'LONGITUDE', 'standard_name','longitude');
ncwriteatt(outfile,'LONGITUDE', 'units','degree_east');
ncwriteatt(outfile,'LONGITUDE', 'axis','X');
ncwriteatt(outfile,'LONGITUDE', 'long_name','longitude of measurement');
ncwriteatt(outfile,'LONGITUDE', 'reference','WGS84');
ncwriteatt(outfile,'LONGITUDE', 'coordinate_reference_frame','urn:ogc:crs:EPSG::4326');
ncwriteatt(outfile,'LONGITUDE', 'valid_min',-180);
ncwriteatt(outfile,'LONGITUDE', 'valid_max',180);
ncwriteatt(outfile,'LONGITUDE', 'QC_indicator','good data');
ncwriteatt(outfile,'LONGITUDE', 'processing_level','Data manually reviewed');
ncwriteatt(outfile,'LONGITUDE', 'uncertainty',.004);	% CHANGE
% ncwriteatt(outfile,'LONGITUDE', 'comment','N/A');	% CHANGE



nccreate(outfile,'INSTRDEPTH', 'Dimensions',{'INSTRDEPTH',1}, 'Datatype','single');
ncwrite(outfile,'INSTRDEPTH', single(instrument_depth));
ncwriteatt(outfile,'INSTRDEPTH', 'standard_name','depth');
ncwriteatt(outfile,'INSTRDEPTH', 'units','meter');
ncwriteatt(outfile,'INSTRDEPTH', 'positive','down');
ncwriteatt(outfile,'INSTRDEPTH', 'axis','Z');
ncwriteatt(outfile,'INSTRDEPTH', 'reference','sea_level');
ncwriteatt(outfile,'INSTRDEPTH', 'coordinate_reference_frame','urn:ogc:crs:EPSG::5831');
ncwriteatt(outfile,'INSTRDEPTH', 'long_name','nominal depth of the instrument');
ncwriteatt(outfile,'INSTRDEPTH', 'valid_min',0);
ncwriteatt(outfile,'INSTRDEPTH', 'valid_max',12000);
ncwriteatt(outfile,'INSTRDEPTH', 'QC_indicator','nominal value');
ncwriteatt(outfile,'INSTRDEPTH', 'processing_level','Data manually reviewed');
ncwriteatt(outfile,'INSTRDEPTH', 'uncertainty',1);	% CHANGE
ncwriteatt(outfile,'INSTRDEPTH', 'comment','Depth calculated from mooring diagram. Use PRES to derive time-varing depths of instruments.');	% CHANGE
 

nccreate(outfile,'BINDEPTH', 'Dimensions',{'BINDEPTH',binnber}, 'Datatype','single');
ncwrite(outfile,'BINDEPTH', single(bin_depth));
ncwriteatt(outfile,'BINDEPTH', 'standard_name','depth');
ncwriteatt(outfile,'BINDEPTH', 'units','meter');
ncwriteatt(outfile,'BINDEPTH', 'positive','down');
ncwriteatt(outfile,'BINDEPTH', 'axis','Z');
ncwriteatt(outfile,'BINDEPTH', 'reference','sea_level');
ncwriteatt(outfile,'BINDEPTH', 'coordinate_reference_frame','urn:ogc:crs:EPSG::5831');
ncwriteatt(outfile,'BINDEPTH', 'long_name','ADCP bin nominal depth');
ncwriteatt(outfile,'BINDEPTH', 'valid_min',0);
ncwriteatt(outfile,'BINDEPTH', 'valid_max',12000);
ncwriteatt(outfile,'BINDEPTH', 'QC_indicator','nominal value');
ncwriteatt(outfile,'BINDEPTH', 'processing_level','Data manually reviewed');
ncwriteatt(outfile,'BINDEPTH', 'uncertainty',1);	% CHANGE


% -----------------------------------+
%   Write data variables             |
% -----------------------------------+

nccreate(outfile,'PRES',  'Dimensions',{'BINDEPTH', binnber, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'PRES', single(pressure));
ncwriteatt(outfile,'PRES', 'standard_name','sea_water_pressure');
ncwriteatt(outfile,'PRES', 'units','decibar');
ncwriteatt(outfile,'PRES', 'reference','sea_level');
ncwriteatt(outfile,'PRES', 'long_name','pressure of measurement');
ncwriteatt(outfile,'PRES', 'coordinates','TIME BINDEPTH');
ncwriteatt(outfile,'PRES', '_FillValue',single(99999));
ncwriteatt(outfile,'PRES', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'PRES', 'valid_max',single(12000));  % CHANGE
ncwriteatt(outfile,'PRES', 'QC_indicator','good data');
ncwriteatt(outfile,'PRES', 'processing_level','Data manually reviewed');
ncwriteatt(outfile,'PRES', 'uncertainty',1);	% CHANGE



nccreate(outfile,'VCUR', 'Dimensions',{'BINDEPTH', binnber, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'VCUR', single(vel_north));
ncwriteatt(outfile,'VCUR', 'standard_name','northward_sea_water_velocity');
ncwriteatt(outfile,'VCUR','units','meter/second');
ncwriteatt(outfile,'VCUR', '_FillValue',single(99999));
ncwriteatt(outfile,'VCUR', 'coordinates','TIME BINDEPTH');
ncwriteatt(outfile,'VCUR', 'long_name','current north component');
ncwriteatt(outfile,'VCUR', 'QC_indicator','good data');
ncwriteatt(outfile,'VCUR', 'processing_level',adcpinfo.proclvl);
%ncwriteatt(outfile,'VCUR', 'valid_min',single(-5));  % CHANGE	
%ncwriteatt(outfile,'VCUR', 'valid_max',single(5));  % CHANGE	
ncwriteatt(outfile,'VCUR', 'accuracy',adcpinfo.velaccuracy);	% CHANGE
ncwriteatt(outfile,'VCUR', 'resolution',adcpinfo.velresolution);	% CHANGE
%ncwriteatt(outfile,'VCUR', 'cell_methods', ('TIME: mean (interval: 60 min) BINDEPTH: mean (interval: 10m) LATITUDE: point LONGITUDE: point'));   %CHANGE 
ncwriteatt(outfile,'VCUR', 'DM_indicator','D');
ncwriteatt(outfile,'VCUR', 'sensor_model',adcpinfo.sensor_model);	% CHANGE
ncwriteatt(outfile,'VCUR', 'sensor_manufacturer',adcpinfo.sensor_manufacturer);	% CHANGE
ncwriteatt(outfile,'VCUR', 'sensor_reference',adcpinfo.sensor_ref);	% CHANGE
ncwriteatt(outfile,'VCUR', 'sensor_serial_number',adcpinfo.serial_num);
ncwriteatt(outfile,'VCUR', 'sensor_mount',adcpinfo.sensor_mount);
ncwriteatt(outfile,'VCUR', 'sensor_orientation',adcpinfo.sens_orientation);  %CHANGE
ncwriteatt(outfile,'VCUR', 'coordinate_system',adcpinfo.coordsyst);
ncwriteatt(outfile,'VCUR', 'Distance_to_first_bin',adcpinfo.dist1stbin);   %CHANGE 
ncwriteatt(outfile,'VCUR', 'Bin_size',adcpinfo.binsize);	% CHANGE
ncwriteatt(outfile,'VCUR', 'Number_of_bins',adcpinfo.nbbins);	% CHANGE
ncwriteatt(outfile,'VCUR', 'pings_per_ensemble',adcpinfo.pingperens);	% CHANGE
ncwriteatt(outfile,'VCUR', 'Time_between_ping_groups',adcpinfo.timepinggroup);	% CHANGE
ncwriteatt(outfile,'VCUR', 'comment','excludes_bins_affected_by_seafloor_interference_and_side_lobes');	% CHANGE



nccreate(outfile,'UCUR', 'Dimensions',{'BINDEPTH', binnber, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'UCUR', single(vel_east));
ncwriteatt(outfile,'UCUR', 'standard_name','eastward_sea_water_velocity');
ncwriteatt(outfile,'UCUR','units','meter/second');
ncwriteatt(outfile,'UCUR', '_FillValue',single(99999));
ncwriteatt(outfile,'UCUR', 'coordinates','TIME BINDEPTH');
ncwriteatt(outfile,'UCUR', 'long_name','current east component');
ncwriteatt(outfile,'UCUR', 'QC_indicator','good data');
ncwriteatt(outfile,'UCUR', 'processing_level',adcpinfo.proclvl);
%ncwriteatt(outfile,'UCUR', 'valid_min',single(-5));  % CHANGE	
%ncwriteatt(outfile,'UCUR', 'valid_max',single(5));  % CHANGE	
ncwriteatt(outfile,'UCUR', 'accuracy',adcpinfo.velaccuracy);	% CHANGE
ncwriteatt(outfile,'UCUR', 'resolution',adcpinfo.velresolution);	% CHANGE
%ncwriteatt(outfile,'UCUR', 'cell_methods', ('TIME: mean (interval: 60 min) BINDEPTH: mean (interval: 10m) LATITUDE: point LONGITUDE: point'));   %CHANGE 
ncwriteatt(outfile,'UCUR', 'DM_indicator','D');
ncwriteatt(outfile,'UCUR', 'sensor_model',adcpinfo.sensor_model);	% CHANGE
ncwriteatt(outfile,'UCUR', 'sensor_manufacturer',adcpinfo.sensor_manufacturer);	% CHANGE
ncwriteatt(outfile,'UCUR', 'sensor_reference',adcpinfo.sensor_ref);	% CHANGE
ncwriteatt(outfile,'UCUR', 'sensor_serial_number',adcpinfo.serial_num);
ncwriteatt(outfile,'UCUR', 'sensor_mount',adcpinfo.sensor_mount);
ncwriteatt(outfile,'UCUR', 'sensor_orientation',adcpinfo.sens_orientation);  %CHANGE
ncwriteatt(outfile,'UCUR', 'coordinate_system',adcpinfo.coordsyst);
ncwriteatt(outfile,'UCUR', 'Distance_to_first_bin',adcpinfo.dist1stbin);   %CHANGE 
ncwriteatt(outfile,'UCUR', 'Bin_size',adcpinfo.binsize);	% CHANGE
ncwriteatt(outfile,'UCUR', 'Number_of_bins',adcpinfo.nbbins);	% CHANGE
ncwriteatt(outfile,'UCUR', 'pings_per_ensemble',adcpinfo.pingperens);	% CHANGE
ncwriteatt(outfile,'UCUR', 'Time_between_ping_groups',adcpinfo.timepinggroup);	% CHANGE
ncwriteatt(outfile,'UCUR', 'comment','excludes_bins_affected_by_seafloor_interference_and_side_lobes');	% CHANGE



nccreate(outfile,'WCUR', 'Dimensions',{'BINDEPTH', binnber, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'WCUR', single(vel_vert));
ncwriteatt(outfile,'WCUR', 'standard_name','vertical_sea_water_velocity');
ncwriteatt(outfile,'WCUR','units','meter/second');
ncwriteatt(outfile,'WCUR', '_FillValue',single(99999));
ncwriteatt(outfile,'WCUR', 'coordinates','TIME BINDEPTH');
ncwriteatt(outfile,'WCUR', 'long_name','current vertical component');
ncwriteatt(outfile,'WCUR', 'QC_indicator','good data');
ncwriteatt(outfile,'WCUR', 'processing_level',adcpinfo.proclvl);
%ncwriteatt(outfile,'WCUR', 'valid_min',single(-5));  % CHANGE	
%ncwriteatt(outfile,'WCUR', 'valid_max',single(5));  % CHANGE	
ncwriteatt(outfile,'WCUR', 'accuracy',adcpinfo.velaccuracy);	% CHANGE
ncwriteatt(outfile,'WCUR', 'resolution',adcpinfo.velresolution);	% CHANGE
%ncwriteatt(outfile,'WCUR', 'cell_methods', ('TIME: mean (interval: 60 min) BINDEPTH: mean (interval: 10m) LATITUDE: point LONGITUDE: point'));   %CHANGE 
ncwriteatt(outfile,'WCUR', 'DM_indicator','D');
ncwriteatt(outfile,'WCUR', 'sensor_model',adcpinfo.sensor_model);	% CHANGE
ncwriteatt(outfile,'WCUR', 'sensor_manufacturer',adcpinfo.sensor_manufacturer);	% CHANGE
ncwriteatt(outfile,'WCUR', 'sensor_reference',adcpinfo.sensor_ref);	% CHANGE
ncwriteatt(outfile,'WCUR', 'sensor_serial_number',adcpinfo.serial_num);
ncwriteatt(outfile,'WCUR', 'sensor_mount',adcpinfo.sensor_mount);
ncwriteatt(outfile,'WCUR', 'sensor_orientation',adcpinfo.sens_orientation);  %CHANGE
ncwriteatt(outfile,'WCUR', 'coordinate_system',adcpinfo.coordsyst);
ncwriteatt(outfile,'WCUR', 'Distance_to_first_bin',adcpinfo.dist1stbin);   %CHANGE 
ncwriteatt(outfile,'WCUR', 'Bin_size',adcpinfo.binsize);	% CHANGE
ncwriteatt(outfile,'WCUR', 'Number_of_bins',adcpinfo.nbbins);	% CHANGE
ncwriteatt(outfile,'WCUR', 'pings_per_ensemble',adcpinfo.pingperens);	% CHANGE
ncwriteatt(outfile,'WCUR', 'Time_between_ping_groups',adcpinfo.timepinggroup);	% CHANGE
ncwriteatt(outfile,'WCUR', 'comment','excludes_bins_affected_by_seafloor_interference_and_side_lobes');	% CHANGE



nccreate(outfile,'ECUR', 'Dimensions',{'BINDEPTH', binnber, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'ECUR', single(vel_err));
ncwriteatt(outfile,'ECUR','units','meter/second');
ncwriteatt(outfile,'ECUR', '_FillValue',single(99999));
ncwriteatt(outfile,'ECUR', 'coordinates','TIME BINDEPTH');
ncwriteatt(outfile,'ECUR', 'long_name','error_sea_water_velocity');
ncwriteatt(outfile,'ECUR', 'QC_indicator','good data');
ncwriteatt(outfile,'ECUR', 'processing_level',adcpinfo.proclvl);
%ncwriteatt(outfile,'ECUR', 'valid_min',single(-5));  % CHANGE	
%ncwriteatt(outfile,'ECUR', 'valid_max',single(5));  % CHANGE	
ncwriteatt(outfile,'ECUR', 'DM_indicator','D');
ncwriteatt(outfile,'ECUR', 'sensor_model',adcpinfo.sensor_model);	% CHANGE
ncwriteatt(outfile,'ECUR', 'sensor_manufacturer',adcpinfo.sensor_manufacturer);	% CHANGE
ncwriteatt(outfile,'ECUR', 'sensor_reference',adcpinfo.sensor_ref);	% CHANGE
ncwriteatt(outfile,'ECUR', 'sensor_serial_number',adcpinfo.serial_num);
ncwriteatt(outfile,'ECUR', 'sensor_mount',adcpinfo.sensor_mount);
ncwriteatt(outfile,'ECUR', 'sensor_orientation',adcpinfo.sens_orientation);  %CHANGE
ncwriteatt(outfile,'ECUR', 'coordinate_system',adcpinfo.coordsyst);
ncwriteatt(outfile,'ECUR', 'Distance_to_first_bin',adcpinfo.dist1stbin);   %CHANGE 
ncwriteatt(outfile,'ECUR', 'Bin_size',adcpinfo.binsize);	% CHANGE
ncwriteatt(outfile,'ECUR', 'Number_of_bins',adcpinfo.nbbins);	% CHANGE
ncwriteatt(outfile,'ECUR', 'pings_per_ensemble',adcpinfo.pingperens);	% CHANGE
ncwriteatt(outfile,'ECUR', 'Time_between_ping_groups',adcpinfo.timepinggroup);	% CHANGE
ncwriteatt(outfile,'ECUR', 'comment','excludes_bins_affected_by_seafloor_interference_and_side_lobes');	% CHANGE





% nccreate(outfile,'MAGVAR', 'Dimensions',{'INSTR_DEPTH',InstrDepthDim, 'TIME',TimeDim}, 'Datatype','single');
% ncwrite(outfile,'MAGVAR', single(dev2));
% ncwriteatt(outfile,'MAGVAR', 'long_name','magnetic_compass_correction_applied');
% ncwriteatt(outfile,'MAGVAR','units','degree');
% ncwriteatt(outfile,'MAGVAR', '_FillValue',single(99999));
% ncwriteatt(outfile,'MAGVAR', 'coordinates','TIME INSTR_DEPTH');
% ncwriteatt(outfile,'MAGVAR', 'QC_indicator','good data');
% ncwriteatt(outfile,'MAGVAR', 'model','IGRF11 via NGDC');	% CHANGE
% ncwriteatt(outfile,'MAGVAR', 'processing_level','Data manually reviewed');
% ncwriteatt(outfile,'MAGVAR', 'valid_min',single(-180));
% ncwriteatt(outfile,'MAGVAR', 'valid_max',single(180));
% ncwriteatt(outfile,'MAGVAR', 'comment','East is positive. Monthly values obtained with geomag70 software, linear interpolation in between.');	% CHANGE
% 
% 
% 
% nccreate(outfile,'EAMP', 'Dimensions',{'BINDEPTH', binnber,  'BEAM', 4, 'TIME',TimeDim}, 'Datatype','single');
% ncwrite(outfile,'EAMP', single(echo));
% ncwriteatt(outfile,'EAMP','units','count');
% ncwriteatt(outfile,'EAMP', '_FillValue',single(99999));
% ncwriteatt(outfile,'EAMP', 'coordinates','TIME BEAM BINDEPTH');
% ncwriteatt(outfile,'EAMP', 'long_name','echo amplitude');
% ncwriteatt(outfile,'EAMP', 'QC_indicator','good data');
% ncwriteatt(outfile,'EAMP', 'processing_level','Data manually reviewed');
% ncwriteatt(outfile,'EAMP', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'EAMP', 'valid_max',single(200));  % CHANGE
% ncwriteatt(outfile,'EAMP', 'accuracy',1);	% CHANGE
% ncwriteatt(outfile,'EAMP', 'resolution',1);	% CHANGE
% ncwriteatt(outfile,'EAMP', 'cell_methods', ('TIME: mean (interval: 60 min) BINDEPTH: mean (interval: 10m) LATITUDE: point LONGITUDE: point'));	% CHANGE
% ncwriteatt(outfile,'EAMP', 'DM_indicator','D');
% ncwriteatt(outfile,'EAMP', 'sensor_model','Workhorse Long Ranger 75 kHz');	% CHANGE
% ncwriteatt(outfile,'EAMP', 'sensor_manufacturer','Teledyne RDI');	% CHANGE
% ncwriteatt(outfile,'EAMP', 'sensor_reference','http://www.rdinstruments.com/datasheets/long_ranger_datasheet_lr.pdf');	% CHANGE
% ncwriteatt(outfile,'EAMP', 'sensor_serial_number',serial_num);
% ncwriteatt(outfile,'EAMP', 'sensor_mount','mounted_on_mooring_line');
% ncwriteatt(outfile,'EAMP', 'sensor_orientation','upward');	% CHANGE
% ncwriteatt(outfile,'EAMP', 'Distance_to_first_bin',18.77);	% CHANGE
% ncwriteatt(outfile,'EAMP', 'Bin_size',10);	% CHANGE
% ncwriteatt(outfile,'EAMP', 'Number_of_bins',42);	% CHANGE
% ncwriteatt(outfile,'EAMP', 'comment','excludes_bins_affected_by_seafloor_interference_and_side_lobes');	% CHANGE
% ncwriteatt(outfile,'EAMP', 'pings_per_ensemble',30);	% CHANGE
% ncwriteatt(outfile,'EAMP', 'Time_between_ping_groups','120s');	% CHANGE
% 
% 
% 
% nccreate(outfile,'CMAG', 'Dimensions',{'BINDEPTH', binnber,  'BEAM', 4, 'TIME',TimeDim}, 'Datatype','single');
% ncwrite(outfile,'CMAG', single(corr));
% ncwriteatt(outfile,'CMAG','units','count');
% ncwriteatt(outfile,'CMAG', '_FillValue',single(99999));
% ncwriteatt(outfile,'CMAG', 'coordinates','TIME BEAM BINDEPTH');
% ncwriteatt(outfile,'CMAG', 'long_name','correlation magnitude');
% ncwriteatt(outfile,'CMAG', 'QC_indicator','good data');
% ncwriteatt(outfile,'CMAG', 'processing_level','Data manually reviewed');
% ncwriteatt(outfile,'CMAG', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'CMAG', 'valid_max',single(160));  % CHANGE
% ncwriteatt(outfile,'CMAG', 'accuracy',1);	% CHANGE
% ncwriteatt(outfile,'CMAG', 'resolution',1);	% CHANGE
% ncwriteatt(outfile,'CMAG', 'cell_methods', ('TIME: mean (interval: 60 min) BINDEPTH: mean (interval: 10m) LATITUDE: point LONGITUDE: point'));	% CHANGE
% ncwriteatt(outfile,'CMAG', 'DM_indicator','D');
% ncwriteatt(outfile,'CMAG', 'sensor_model','Workhorse Long Ranger 75 kHz');	% CHANGE
% ncwriteatt(outfile,'CMAG', 'sensor_manufacturer','Teledyne RDI');	% CHANGE
% ncwriteatt(outfile,'CMAG', 'sensor_reference','http://www.rdinstruments.com/datasheets/long_ranger_datasheet_lr.pdf');	% CHANGE
% ncwriteatt(outfile,'CMAG', 'sensor_serial_number',serial_num);
% ncwriteatt(outfile,'CMAG', 'sensor_mount','mounted_on_mooring_line');
% ncwriteatt(outfile,'CMAG', 'sensor_orientation','upward');	% CHANGE
% ncwriteatt(outfile,'CMAG', 'Distance_to_first_bin',18.77);	% CHANGE
% ncwriteatt(outfile,'CMAG', 'Bin_size',10);	% CHANGE
% ncwriteatt(outfile,'CMAG', 'Number_of_bins',42);	% CHANGE
% ncwriteatt(outfile,'CMAG', 'comment','excludes_bins_affected_by_seafloor_interference_and_side_lobes');	% CHANGE
% ncwriteatt(outfile,'CMAG', 'pings_per_ensemble',30);	% CHANGE
% ncwriteatt(outfile,'CMAG', 'Time_between_ping_groups','120s');	% CHANGE
% 
% 
% 
% nccreate(outfile,'PERG', 'Dimensions',{'BINDEPTH', binnber,  'BEAM', 4, 'TIME',TimeDim}, 'Datatype','single');
% ncwrite(outfile,'PERG', single(perg));
% ncwriteatt(outfile,'PERG','units','percent');
% ncwriteatt(outfile,'PERG', '_FillValue',single(99999));
% ncwriteatt(outfile,'PERG', 'coordinates','TIME BEAM BINDEPTH');
% ncwriteatt(outfile,'PERG', 'long_name','percentage_of_good_solutions');
% ncwriteatt(outfile,'PERG', 'QC_indicator','good data');
% ncwriteatt(outfile,'PERG', 'processing_level','Data manually reviewed');
% ncwriteatt(outfile,'PERG', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'PERG', 'valid_max',single(100));  % CHANGE
% ncwriteatt(outfile,'PERG', 'accuracy',1);	% CHANGE
% ncwriteatt(outfile,'PERG', 'resolution',1);	% CHANGE
% ncwriteatt(outfile,'PERG', 'cell_methods', ('TIME: mean (interval: 60 min) BINDEPTH: mean (interval: 10m) LATITUDE: point LONGITUDE: point'));	% CHANGE
% ncwriteatt(outfile,'PERG', 'DM_indicator','D');
% ncwriteatt(outfile,'PERG', 'sensor_model','Workhorse Long Ranger 75 kHz');	% CHANGE
% ncwriteatt(outfile,'PERG', 'sensor_manufacturer','Teledyne RDI');	% CHANGE
% ncwriteatt(outfile,'PERG', 'sensor_reference','http://www.rdinstruments.com/datasheets/long_ranger_datasheet_lr.pdf');	% CHANGE
% ncwriteatt(outfile,'PERG', 'sensor_serial_number',serial_num);
% ncwriteatt(outfile,'PERG', 'sensor_mount','mounted_on_mooring_line');
% ncwriteatt(outfile,'PERG', 'sensor_orientation','upward');	% CHANGE
% ncwriteatt(outfile,'PERG', 'Distance_to_first_bin',18.77);	% CHANGE
% ncwriteatt(outfile,'PERG', 'Bin_size',10);	% CHANGE
% ncwriteatt(outfile,'PERG', 'Number_of_bins',42);	% CHANGE
% ncwriteatt(outfile,'PERG', 'comment','excludes_bins_affected_by_seafloor_interference_and_side_lobes');	% CHANGE
% ncwriteatt(outfile,'PERG', 'pings_per_ensemble',30);	% CHANGE
% ncwriteatt(outfile,'PERG', 'Time_between_ping_groups','120s');	% CHANGE
% 
% 
% 
% nccreate(outfile,'PITCH', 'Dimensions',{'INSTR_DEPTH',InstrDepthDim, 'TIME',TimeDim}, 'Datatype','single');
% ncwrite(outfile,'PITCH', single(pitch_series));
% ncwriteatt(outfile,'PITCH', 'units','degree');
% ncwriteatt(outfile,'PITCH', '_FillValue',single(99999));
% ncwriteatt(outfile,'PITCH', 'coordinates','TIME INSTR_DEPTH');
% ncwriteatt(outfile,'PITCH', 'long_name','instrument_pitch')
% ncwriteatt(outfile,'PITCH', 'QC_indicator','good data');
% ncwriteatt(outfile,'PITCH', 'processing_level','Data manually reviewed');
% ncwriteatt(outfile,'PITCH', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'PITCH', 'valid_max',single(90));  % CHANGE
% ncwriteatt(outfile,'PITCH', 'uncertainty',1);	% CHANGE
% ncwriteatt(outfile,'PITCH', 'accuracy',0.5);	% CHANGE
% ncwriteatt(outfile,'PITCH', 'precision',1);	% CHANGE
% ncwriteatt(outfile,'PITCH', 'resolution',0.01);	% CHANGE
% ncwriteatt(outfile,'PITCH', 'DM_indicator','D');
% 
% 
% 
% nccreate(outfile,'ROLL', 'Dimensions',{'INSTR_DEPTH',InstrDepthDim, 'TIME',TimeDim}, 'Datatype','single');
% ncwrite(outfile,'ROLL', single(roll_series));
% ncwriteatt(outfile,'ROLL', 'units','degree');
% ncwriteatt(outfile,'ROLL', '_FillValue',single(99999));
% ncwriteatt(outfile,'ROLL', 'coordinates','TIME INSTR_DEPTH');
% ncwriteatt(outfile,'ROLL', 'long_name','instrument_pitch')
% ncwriteatt(outfile,'ROLL', 'QC_indicator','good data');
% ncwriteatt(outfile,'ROLL', 'processing_level','Data manually reviewed');
% ncwriteatt(outfile,'ROLL', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'ROLL', 'valid_max',single(90));  % CHANGE
% ncwriteatt(outfile,'ROLL', 'uncertainty',1);	% CHANGE
% ncwriteatt(outfile,'ROLL', 'accuracy',0.5);	% CHANGE
% ncwriteatt(outfile,'ROLL', 'precision',1);	% CHANGE
% ncwriteatt(outfile,'ROLL', 'resolution',0.01);	% CHANGE
% ncwriteatt(outfile,'ROLL', 'DM_indicator','D');
% 
% 
% 
% nccreate(outfile,'HEADING', 'Dimensions',{'INSTR_DEPTH',InstrDepthDim, 'TIME',TimeDim}, 'Datatype','single');
% ncwrite(outfile,'HEADING', single(heading_series));
% ncwriteatt(outfile,'HEADING', 'long_name','pressure of each measurement');
% ncwriteatt(outfile,'HEADING', 'standard_name','sea_water_pressure');
% ncwriteatt(outfile,'HEADING', 'units','decibar');
% ncwriteatt(outfile,'HEADING', 'positive','down');
% ncwriteatt(outfile,'HEADING', '_FillValue',single(99999));
% ncwriteatt(outfile,'HEADING', 'coordinates','TIME INSTR_DEPTH');
% ncwriteatt(outfile,'HEADING', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'HEADING', 'valid_max',single(12000));  % CHANGE
% ncwriteatt(outfile,'HEADING', 'QC_indicator',1);	% CHANGE
% ncwriteatt(outfile,'HEADING', 'QC_procedure',5);	% CHANGE
% ncwriteatt(outfile,'HEADING', 'uncertainty',50);	% CHANGE
% ncwriteatt(outfile,'HEADING', 'axis','Z');
% ncwriteatt(outfile,'HEADING', 'reference','mean_sea_level');
% % ncwriteatt(outfile,'HEADING', 'comment','N/A');   % CHANGE
% 
% 
% 
% nccreate(outfile,'TEMP', 'Dimensions',{'INSTR_DEPTH',InstrDepthDim, 'TIME',TimeDim}, 'Datatype','single');
% ncwrite(outfile,'TEMP', single(temp_series));
% ncwriteatt(outfile,'TEMP', 'standard_name','sea_water_temperature');
% ncwriteatt(outfile,'TEMP','units','degree_Celsius');
% ncwriteatt(outfile,'TEMP', '_FillValue',single(99999));
% ncwriteatt(outfile,'TEMP', 'coordinates','TIME INSTR_DEPTH');
% ncwriteatt(outfile,'TEMP', 'long_name','sea water temperature');
% ncwriteatt(outfile,'TEMP', 'QC_indicator','good data');
% ncwriteatt(outfile,'TEMP', 'processing_level','Data manually reviewed');
% ncwriteatt(outfile,'TEMP', 'valid_min',single(-2.5));	% CHANGE;
% ncwriteatt(outfile,'TEMP', 'valid_max',single(40));	% CHANGE
% ncwriteatt(outfile,'TEMP', 'uncertainty',0.01);	% CHANGE
% ncwriteatt(outfile,'TEMP', 'accuracy',0.005);	% CHANGE
% ncwriteatt(outfile,'TEMP', 'precision',0.4);	% CHANGE
% ncwriteatt(outfile,'TEMP', 'resolution',0.01);	% CHANGE
% ncwriteatt(outfile,'TEMP', 'DM_indicator','D');
% ncwriteatt(outfile,'TEMP', 'reference_scale','ITS-90');
% ncwriteatt(outfile,'TEMP', 'sensor_model','Workhorse Long Ranger 75 kHz');	% CHANGE
% ncwriteatt(outfile,'TEMP', 'sensor_manufacturer','Teledyne RDI');	% CHANGE
% ncwriteatt(outfile,'TEMP', 'sensor_reference','http://www.rdinstruments.com/datasheets/long_ranger_datasheet_lr.pdf');	% CHANGE
% ncwriteatt(outfile,'TEMP', 'sensor_serial_number',serial_num);
% ncwriteatt(outfile,'TEMP', 'sensor_mount','mounted_on_mooring_line');
% ncwriteatt(outfile,'TEMP', 'sensor_orientation','horizontal');
% % ncwriteatt(outfile,'TEMP', 'comment','N/A');   % CHANGE








% -------------------+
%   End              |
% -------------------+

fprintf('Finish %s\n', datestr(now));

end
