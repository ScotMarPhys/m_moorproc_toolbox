% ------------------------------------------------------------------
% Copyright, IFREMER, 2012
% Modified by Thomas Richter at NIOZ, 2014
% Adapted for OSNAP CM data by Feili Li at Duke University, 2016
% 
% Last updated: August 12, 2016
% Edited by Loic Houpert, SAMS 18/11/2016
% ------------------------------------------------------------------


function write_CM_to_NetCDF(filename, moorname, latitude, longitude, nortekinfo, depth, time, pressure, vel_east, vel_north, vel_vert)
%
% ------------------------------------------------------------------
% NOTES ON INPUT
%
%   filename-   output nc file name (without file extention '.nc')
%   moorname-  name of the mooring
%   latitude-   latitude of measurement  [degree north]
%   longitude-  longitude of measurement  [degree east]
%   nortekinfo-  structure with nortek metadata  
%   depth-      nominal depth of measurement  [m]
%   time-       time of measurement
%   pressure-   sea water pressure  [db]
%   vel_east-   eastward sea water velocity [m/s]
%   vel_north-  northward sea water velocity [m/s]
%   vel_vert-	vertical sea water veloicty [m/s]
%
% The following parameter(s) is needed when available:
%
%   dev2-	magnetic compass correction [degree]
% ------------------------------------------------------------------



% --------------------+
%   Begin             |
% --------------------+

fprintf('Begin %s\n', datestr(now));

% Dimensions
TimeDim = length(time);
DepthDim = length(depth);

% Open the file to write
outfile = [filename,'.nc'];
nc = netcdf.create(outfile,'CLOBBER');
netcdf.close(nc)

  
  
% ------------------------------+
%   Write global attributes     |
% ------------------------------+


% === DISCOVERY AND IDENTIFICATION ===
ncwriteatt(outfile,'/','site_code', 'OSNAP');  
ncwriteatt(outfile,'/','platform_code', ['OSNAP-' moorname]);	% CHANGE
ncwriteatt(outfile,'/','data_mode', 'D');   
ncwriteatt(outfile,'/','title', ['OSNAP CM Data ',datestr(time(1),'mm/yyyy'),'-',datestr(time(end),'mm/yyyy')]); 
ncwriteatt(outfile,'/','summary', ['sea water velocity at nominal depth ',num2str(depth),'m']); 
ncwriteatt(outfile,'/','naming_authority', 'OceanSITES'); 
ncwriteatt(outfile,'/','id', filename);   
ncwriteatt(outfile,'/','source', 'subsurface mooring');
ncwriteatt(outfile,'/','principal_investigator', nortekinfo.principal_investigator);	% CHANGE
ncwriteatt(outfile,'/','principal_investigator_email', nortekinfo.principal_investigator_email);	% CHANGE
ncwriteatt(outfile,'/','principal_investigator_url', nortekinfo.principal_investigator_url);	% CHANGE
ncwriteatt(outfile,'/','institution', nortekinfo.institution);	% CHANGE
ncwriteatt(outfile,'/','project', nortekinfo.project);
ncwriteatt(outfile,'/','array', nortekinfo.array);  
ncwriteatt(outfile,'/','network', nortekinfo.network);  
% ncwriteatt(outfile,'/','comment', 'N/A');  %CHANGE 


    
% === GEO-SPATIAL-TEMPORAL INFO ===    
ncwriteatt(outfile,'/','area', 'North Atlantic Ocean');
ncwriteatt(outfile,'/','geospatial_lat_min', sprintf('%.1f',latitude));
ncwriteatt(outfile,'/','geospatial_lat_max', sprintf('%.1f',latitude));
ncwriteatt(outfile,'/','geospatial_lat_units', 'degree_north');
ncwriteatt(outfile,'/','geospatial_lon_min', sprintf('%.1f',longitude));
ncwriteatt(outfile,'/','geospatial_lon_max', sprintf('%.1f',longitude));
ncwriteatt(outfile,'/','geospatial_lon_units', 'degree_east');
ncwriteatt(outfile,'/','geospatial_vertical_min', depth);
ncwriteatt(outfile,'/','geospatial_vertical_max', depth);
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
ncwriteatt(outfile,'/','netcdf_version','4.3')



% === PUBLICATION INFO ===
ncwriteatt(outfile,'/','publisher_name', 'Feili Li');  %CHANGE     
ncwriteatt(outfile,'/','publisher_email', 'feili.li@duke.edu');  %CHANGE   
ncwriteatt(outfile,'/','publisher_url', 'http://www.o-snap.org');  %CHANGE  
ncwriteatt(outfile,'/','references', 'http://www.o-snap.org, http://www.oceansites.org');
ncwriteatt(outfile,'/','data_assembly_center', 'Duke University');  %CHANGE    
ncwriteatt(outfile,'/','update_interval', 'void');
ncwriteatt(outfile,'/','license', 'Follows CLIVAR (Climate Variability and Predictability) standards, cf. http://www.clivar.org/data/data_policy.php. Data available free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.');
ncwriteatt(outfile,'/','citation', 'These data were collected and made freely available by the OSNAP project and the national programs that contribute to it.');
ncwriteatt(outfile,'/','acknowledgement', 'Funding source: the UK Natural Environment Research Council (NERC), UK OSNAP project');  %CHANGE




% === PROVENANCE ===
ncwriteatt(outfile,'/','date_created', datestr(now+8/24,'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','date_modified', datestr(now+8/24,'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','history', 'Delayed time processed quality controlled');
ncwriteatt(outfile,'/','processing_level',nortekinfo.processing_level);
ncwriteatt(outfile,'/','QC_indicator', nortekinfo.QC_indicator);
ncwriteatt(outfile,'/','contributor_name', nortekinfo.contributor_name);	% CHANGE
ncwriteatt(outfile,'/','contributor_role', nortekinfo.contributor_role);	% CHANGE
ncwriteatt(outfile,'/','contributor_email', nortekinfo.contributor_email);	% CHANGE
   


% --------------------------------------+
%   Write coordinate variables          |
% --------------------------------------+

nccreate(outfile,'TIME', 'Dimensions',{'TIME',TimeDim}, 'Datatype','double');
ncwrite(outfile,'TIME', double(time) - datenum(1950,1,1,0,0,0));
ncwriteatt(outfile,'TIME', 'standard_name','time');
ncwriteatt(outfile,'TIME', 'long_name','time of measurement');
ncwriteatt(outfile,'TIME', 'units','days since 1950-01-01T00:00:00Z');
ncwriteatt(outfile,'TIME', 'axis','T');
ncwriteatt(outfile,'TIME', 'conventions', 'Relative julian days with decimal part (as parts of the day)');
ncwriteatt(outfile,'TIME', 'valid_min', single(0));
ncwriteatt(outfile,'TIME', 'valid_max', single(90000));
ncwriteatt(outfile,'TIME', 'QC_indicator','good data');
ncwriteatt(outfile,'TIME', 'processing_level','Data manually reviewed');



nccreate(outfile,'LATITUDE', 'Dimensions',{'LATITUDE',1}, 'Datatype','single');
ncwrite(outfile,'LATITUDE', single(latitude));
ncwriteatt(outfile,'LATITUDE', 'standard_name','latitude');
ncwriteatt(outfile,'LATITUDE', 'long_name','latitude of measurement');
ncwriteatt(outfile,'LATITUDE', 'units','degree_north');
ncwriteatt(outfile,'LATITUDE', 'axis','Y');
ncwriteatt(outfile,'LATITUDE', 'reference','WGS84');
ncwriteatt(outfile,'LATITUDE', 'coordinate_reference_frame','urn:ogc:crs:EPSG::4326');
ncwriteatt(outfile,'LATITUDE', 'valid_min',single(-90));
ncwriteatt(outfile,'LATITUDE', 'valid_max',single(90));
ncwriteatt(outfile,'LATITUDE', 'QC_indicator','good data');
ncwriteatt(outfile,'LATITUDE', 'processing_level','Data manually reviewed');
ncwriteatt(outfile,'LATITUDE', 'uncertainty',.004);  %CHANGE



nccreate(outfile,'LONGITUDE', 'Dimensions',{'LONGITUDE',1}, 'Datatype','single');
ncwrite(outfile,'LONGITUDE', single(longitude));
ncwriteatt(outfile,'LONGITUDE', 'standard_name','longitude');
ncwriteatt(outfile,'LONGITUDE', 'long_name','longitude of measurement');
ncwriteatt(outfile,'LONGITUDE', 'units','degree_east');
ncwriteatt(outfile,'LONGITUDE', 'axis','X');
ncwriteatt(outfile,'LONGITUDE', 'reference','WGS84');
ncwriteatt(outfile,'LONGITUDE', 'coordinate_reference_frame','urn:ogc:crs:EPSG::4326');
ncwriteatt(outfile,'LONGITUDE', 'valid_min',single(-180));
ncwriteatt(outfile,'LONGITUDE', 'valid_max',single(180));
ncwriteatt(outfile,'LONGITUDE', 'QC_indicator','good data');
ncwriteatt(outfile,'LONGITUDE', 'processing_level','Data manually reviewed');
ncwriteatt(outfile,'LONGITUDE', 'uncertainty',.004);  %CHANGE


nccreate(outfile,'DEPTH', 'Dimensions',{'DEPTH',DepthDim}, 'Datatype','single');
ncwrite(outfile,'DEPTH', single(depth));
ncwriteatt(outfile,'DEPTH', 'standard_name','depth');
ncwriteatt(outfile,'DEPTH', 'long_name','depth of measurement');
ncwriteatt(outfile,'DEPTH', 'units','meter');
ncwriteatt(outfile,'DEPTH', 'positive','down');
ncwriteatt(outfile,'DEPTH', 'axis','Z');
ncwriteatt(outfile,'DEPTH', 'reference','sea_level');
ncwriteatt(outfile,'DEPTH', 'coordinate_reference_frame','urn:ogc:crs:EPSG::5831');
ncwriteatt(outfile,'DEPTH', 'valid_min',single(0));
ncwriteatt(outfile,'DEPTH', 'valid_max',single(12000));
ncwriteatt(outfile,'DEPTH', 'QC_indicator','nominal value');
ncwriteatt(outfile,'DEPTH', 'processing_level','Data manually reviewed');
ncwriteatt(outfile,'DEPTH', 'uncertainty',1);  %CHANGE
ncwriteatt(outfile,'DEPTH', 'comment','These are nominal values. Use PRES to derive time-varing depths of instruments.');  %CHANGE


% -----------------------------------+
%   Write data variables             |
% -----------------------------------+

nccreate(outfile,'PRES', 'Dimensions',{'DEPTH',DepthDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'PRES', single(pressure));
ncwriteatt(outfile,'PRES', 'standard_name','sea_water_pressure');
ncwriteatt(outfile,'PRES', 'long_name','pressure of measurement');
ncwriteatt(outfile,'PRES', 'units','decibar');
ncwriteatt(outfile,'PRES', '_FillValue',single(99999));
ncwriteatt(outfile,'PRES', 'coordinates','TIME DEPTH');
ncwriteatt(outfile,'PRES', 'DM_indicator','D');
ncwriteatt(outfile,'PRES', 'axis','Z');
ncwriteatt(outfile,'PRES', 'positive','down');
ncwriteatt(outfile,'PRES', 'reference_scale','mean_sea_level');
ncwriteatt(outfile,'PRES', 'QC_indicator','good data');
ncwriteatt(outfile,'PRES', 'processing_level', nortekinfo.proclvl  );
ncwriteatt(outfile,'PRES', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'PRES', 'valid_max',single(12000));  % CHANGE
ncwriteatt(outfile,'PRES', 'sensor_model',nortekinfo.sensor_model );   %CHANGE 
ncwriteatt(outfile,'PRES', 'sensor_manufacturer',nortekinfo.sensor_manufacturer);   %CHANGE 
ncwriteatt(outfile,'PRES', 'sensor_reference',nortekinfo.sensor_ref);   %CHANGE
ncwriteatt(outfile,'PRES', 'accuracy',nortekinfo.presaccuracy);     %CHANGE 
ncwriteatt(outfile,'PRES', 'resolution',nortekinfo.presresolution);	%CHANGE 
ncwriteatt(outfile,'PRES', 'sensor_serial_number',nortekinfo.serial_num);
ncwriteatt(outfile,'PRES', 'sensor_mount',nortekinfo.sensor_mount);
ncwriteatt(outfile,'PRES', 'sensor_orientation',nortekinfo.sens_orientation);
% ncwriteatt(outfile,'PRES', 'comment','N/A');  %CHANGE



nccreate(outfile,'VCUR', 'Dimensions',{'DEPTH',DepthDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'VCUR', single(vel_north));
ncwriteatt(outfile,'VCUR', 'standard_name','northward_sea_water_velocity');
ncwriteatt(outfile,'VCUR', 'long_name','current north component');
ncwriteatt(outfile,'VCUR','units','meter/second');
ncwriteatt(outfile,'VCUR', '_FillValue',single(99999));
ncwriteatt(outfile,'VCUR', 'coordinates','TIME DEPTH');
ncwriteatt(outfile,'VCUR', 'QC_indicator','good data');
ncwriteatt(outfile,'VCUR', 'processing_level', nortekinfo.proclvl  );
ncwriteatt(outfile,'VCUR', 'valid_min',single(-5));  % CHANGE
ncwriteatt(outfile,'VCUR', 'valid_max',single(5));  % CHANGE  
ncwriteatt(outfile,'VCUR', 'cell_methods', ['TIME: mean DEPTH: mean', 'LATITUDE: point LONGITUDE: point']);  % CHANGE 
ncwriteatt(outfile,'VCUR', 'DM_indicator','D');
ncwriteatt(outfile,'VCUR', 'sensor_model',nortekinfo.sensor_model );   %CHANGE 
ncwriteatt(outfile,'VCUR', 'sensor_manufacturer',nortekinfo.sensor_manufacturer);   %CHANGE 
ncwriteatt(outfile,'VCUR', 'sensor_reference',nortekinfo.sensor_ref);   %CHANGE
ncwriteatt(outfile,'VCUR', 'accuracy',nortekinfo.velaccuracy);     %CHANGE 
ncwriteatt(outfile,'VCUR', 'resolution',nortekinfo.velresolution);	%CHANGE 
ncwriteatt(outfile,'VCUR', 'sensor_serial_number',nortekinfo.serial_num);
ncwriteatt(outfile,'VCUR', 'sensor_mount',nortekinfo.sensor_mount);
ncwriteatt(outfile,'VCUR', 'sensor_orientation',nortekinfo.sens_orientation);
% ncwriteatt(outfile,'VCUR', 'comment','N/A');  % CHANGE



nccreate(outfile,'UCUR', 'Dimensions',{'DEPTH',DepthDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'UCUR', single(vel_east));
ncwriteatt(outfile,'UCUR', 'standard_name','eastward_sea_water_velocity');
ncwriteatt(outfile,'UCUR', 'long_name','current east component');
ncwriteatt(outfile,'UCUR','units','meter/second');
ncwriteatt(outfile,'UCUR', '_FillValue',single(99999));
ncwriteatt(outfile,'UCUR', 'coordinates','TIME DEPTH');
ncwriteatt(outfile,'UCUR', 'QC_indicator','good data');
ncwriteatt(outfile,'UCUR', 'processing_level', nortekinfo.proclvl  );
ncwriteatt(outfile,'UCUR', 'valid_min',single(-5));  % CHANGE
ncwriteatt(outfile,'UCUR', 'valid_max',single(5));  % CHANGE  
ncwriteatt(outfile,'UCUR', 'cell_methods', ['TIME: mean DEPTH: mean', 'LATITUDE: point LONGITUDE: point']);  % CHANGE 
ncwriteatt(outfile,'UCUR', 'DM_indicator','D');
ncwriteatt(outfile,'UCUR', 'sensor_model',nortekinfo.sensor_model );   %CHANGE 
ncwriteatt(outfile,'UCUR', 'sensor_manufacturer',nortekinfo.sensor_manufacturer);   %CHANGE 
ncwriteatt(outfile,'UCUR', 'sensor_reference',nortekinfo.sensor_ref);   %CHANGE
ncwriteatt(outfile,'UCUR', 'accuracy',nortekinfo.velaccuracy);     %CHANGE 
ncwriteatt(outfile,'UCUR', 'resolution',nortekinfo.velresolution);	%CHANGE 
ncwriteatt(outfile,'UCUR', 'sensor_serial_number',nortekinfo.serial_num);
ncwriteatt(outfile,'UCUR', 'sensor_mount',nortekinfo.sensor_mount);
ncwriteatt(outfile,'UCUR', 'sensor_orientation',nortekinfo.sens_orientation);
% ncwriteatt(outfile,'UCUR', 'comment','N/A');  % CHANGE

nccreate(outfile,'WCUR', 'Dimensions',{'DEPTH',DepthDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'WCUR', single(vel_vert));
ncwriteatt(outfile,'WCUR', 'standard_name','vertical_sea_water_velocity');
ncwriteatt(outfile,'WCUR','units','meter/second');
ncwriteatt(outfile,'WCUR', '_FillValue',single(99999));
ncwriteatt(outfile,'WCUR', 'coordinates','TIME BINDEPTH');
ncwriteatt(outfile,'WCUR', 'long_name','current vertical component');
ncwriteatt(outfile,'WCUR', 'QC_indicator','good data');
ncwriteatt(outfile,'WCUR', 'processing_level', nortekinfo.proclvl  );
ncwriteatt(outfile,'WCUR', 'valid_min',single(-5));  % CHANGE
ncwriteatt(outfile,'WCUR', 'valid_max',single(5));  % CHANGE  
ncwriteatt(outfile,'WCUR', 'cell_methods', ['TIME: mean DEPTH: mean', 'LATITUDE: point LONGITUDE: point']);  % CHANGE 
ncwriteatt(outfile,'WCUR', 'DM_indicator','D');
ncwriteatt(outfile,'WCUR', 'sensor_model',nortekinfo.sensor_model );   %CHANGE 
ncwriteatt(outfile,'WCUR', 'sensor_manufacturer',nortekinfo.sensor_manufacturer);   %CHANGE 
ncwriteatt(outfile,'WCUR', 'sensor_reference',nortekinfo.sensor_ref);   %CHANGE
ncwriteatt(outfile,'WCUR', 'accuracy',nortekinfo.velaccuracy);     %CHANGE 
ncwriteatt(outfile,'WCUR', 'resolution',nortekinfo.velresolution);	%CHANGE 
ncwriteatt(outfile,'WCUR', 'sensor_serial_number',nortekinfo.serial_num);
ncwriteatt(outfile,'WCUR', 'sensor_mount',nortekinfo.sensor_mount);
ncwriteatt(outfile,'WCUR', 'sensor_orientation',nortekinfo.sens_orientation);


% nccreate(outfile,'MAGVAR', 'Dimensions',{'DEPTH',DepthDim, 'TIME',TimeDim}, 'Datatype','single');
% ncwrite(outfile,'MAGVAR', single(dev2));
% ncwriteatt(outfile,'MAGVAR', 'long_name','magnetic_compass_correction_applied');
% ncwriteatt(outfile,'MAGVAR','units','degree');
% ncwriteatt(outfile,'MAGVAR', '_FillValue',single(99999));
% ncwriteatt(outfile,'MAGVAR', 'coordinates','TIME DEPTH');
% ncwriteatt(outfile,'MAGVAR', 'QC_indicator','good data');
% ncwriteatt(outfile,'MAGVAR', 'model','IGRF11 via NGDC');  % CHANGE
% ncwriteatt(outfile,'MAGVAR', 'processing_level','Data manually reviewed');
% ncwriteatt(outfile,'MAGVAR', 'valid_min',single(-180));
% ncwriteatt(outfile,'MAGVAR', 'valid_max',single(180));
% ncwriteatt(outfile,'MAGVAR', 'comment','East is positive. Monthly values obtained with geomag70 software, linear interpolation in between');  % CHANGE
% 






% -------------------+
%   End              |
% -------------------+

fprintf('Finish %s\n', datestr(now));






end
