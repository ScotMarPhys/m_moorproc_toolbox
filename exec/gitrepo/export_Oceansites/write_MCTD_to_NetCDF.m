% -------------------------------------------------------------------------
% Copyright, IFREMER, 2012
% Modified by Thomas Richter at NIOZ, 2014
% Adapted for OSNAP MCTD data by Feili Li at Duke University, 2016
% 
% Last updated: August 12, 2016
% Edited by Loic Houpert, SAMS 18/11/2016
% -------------------------------------------------------------------------


function write_MCTD_to_NetCDF(filename, moorname, latitude, longitude, mcatinfo, depth, time, pressure, salinity, temperature, conductivity)
% 
% ------------------------------------------------------------------
% NOTES ON INPUT
%
%   filename-   output nc file name (without file extention '.nc')
%   moorname-  name of the mooring
%   latitude-   latitude of measurement  [degree north]
%   longitude-  longitude of measurement  [degree east]
%   mcatinfo-  structure with microcat metadata  
%   depth-      nominal depth of measurement  [m]
%   time-       time of measurement
%   pressure-	sea water pressure  [db]
%   salinty-	practial salinity derived from conductivity
%   tempearture-	sea water temperature  [degree C]
%   conductivity- sea water electrical conductivity [mS/cm]
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
ncwriteatt(outfile,'/','platform_code', ['OSNAP-' moorname]);  %CHANGE
ncwriteatt(outfile,'/','data_mode', 'D');   
ncwriteatt(outfile,'/','title', ['OSNAP MCTD Data ',datestr(time(1),'mm/yyyy'),'-',datestr(time(end),'mm/yyyy')]); 
ncwriteatt(outfile,'/','summary', ['Water temperature and salinity at nominal depth ',num2str(depth),'m']); 
ncwriteatt(outfile,'/','naming_authority', 'OceanSITES'); 
ncwriteatt(outfile,'/','id', filename);   
ncwriteatt(outfile,'/','source', 'subsurface mooring');
ncwriteatt(outfile,'/','principal_investigator', mcatinfo.principal_investigator);	% CHANGE
ncwriteatt(outfile,'/','principal_investigator_email', mcatinfo.principal_investigator_email);	% CHANGE
ncwriteatt(outfile,'/','principal_investigator_url', mcatinfo.principal_investigator_url);	% CHANGE
ncwriteatt(outfile,'/','institution', mcatinfo.institution);	% CHANGE
ncwriteatt(outfile,'/','project', mcatinfo.project);
ncwriteatt(outfile,'/','array', mcatinfo.array);  
ncwriteatt(outfile,'/','network', mcatinfo.network);  
% ncwriteatt(outfile,'/','comment', 'N/A');  %CHANGE 
    


% === GEO-SPATIAL-TEMPORAL ===
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



% === CONVENTIONS USED ===
ncwriteatt(outfile,'/','format_version', '1.3');
ncwriteatt(outfile,'/','Conventions','CF-1.6 OceanSITES-1.3 ACDD-1.2')
ncwriteatt(outfile,'/','netcdf_version','4.3');



% === PUBLICATION INFORMATION ===
ncwriteatt(outfile,'/','publisher_name', 'Yao Fu');  %CHANGE     
ncwriteatt(outfile,'/','publisher_email', 'yaofu@gatech.edu');  %CHANGE   
ncwriteatt(outfile,'/','publisher_url', 'http://www.o-snap.org');  %CHANGE  
ncwriteatt(outfile,'/','references', 'http://www.o-snap.org, http://www.oceansites.org');
ncwriteatt(outfile,'/','data_assembly_center', 'Georgia Institute of Technology');  %CHANGE    
ncwriteatt(outfile,'/','update_interval', 'void');
ncwriteatt(outfile,'/','license', 'Follows CLIVAR (Climate Variability and Predictability) standards, cf. http://www.clivar.org/data/data_policy.php. Data available free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.');
ncwriteatt(outfile,'/','citation', 'These data were collected and made freely available by the CLASS project and the national programs that contribute to it.');
ncwriteatt(outfile,'/','acknowledgement', 'Funding source: the UK Natural Environment Research Council (NERC), CLASS Project NE/R015953/1');  %CHANGE


% === PROVENANCE ===
ncwriteatt(outfile,'/','date_created', datestr(now+8/24,'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','date_modified', datestr(now+8/24,'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','history', 'Delayed time processed quality controlled');
ncwriteatt(outfile,'/','processing_level',mcatinfo.processing_level);
ncwriteatt(outfile,'/','QC_indicator', mcatinfo.QC_indicator);
ncwriteatt(outfile,'/','contributor_name', mcatinfo.contributor_name);	% CHANGE
ncwriteatt(outfile,'/','contributor_role', mcatinfo.contributor_role);	% CHANGE
ncwriteatt(outfile,'/','contributor_email', mcatinfo.contributor_email);	% CHANGE
   


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
ncwriteatt(outfile,'LATITUDE', 'coordinate_reference_frame','urn:ogc:def:crs:EPSG::4326');
ncwriteatt(outfile,'LATITUDE', 'valid_min', single(-90));
ncwriteatt(outfile,'LATITUDE', 'valid_max', single(90));
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
ncwriteatt(outfile,'LONGITUDE', 'coordinate_reference_frame','urn:ogc:def:crs:EPSG::4326');
ncwriteatt(outfile,'LONGITUDE', 'valid_min', single(-180));
ncwriteatt(outfile,'LONGITUDE', 'valid_max', single(180));
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
ncwriteatt(outfile,'DEPTH', 'valid_min', single(0));
ncwriteatt(outfile,'DEPTH', 'valid_max', single(12000));
ncwriteatt(outfile,'DEPTH', 'QC_indicator','nominal value');
ncwriteatt(outfile,'DEPTH', 'processing_level','Data manually reviewed');
ncwriteatt(outfile,'DEPTH', 'uncertainty',1);  %CHANGE
ncwriteatt(outfile,'DEPTH', 'comment','These are nominal values. Use PRES to derive time-varing depths of instrument.');  %CHANGE


  
% -----------------------------------+
%   Write data variables             |
% -----------------------------------+


nccreate(outfile,'PRES', 'Dimensions',{'DEPTH',DepthDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'PRES', single(pressure));
ncwriteatt(outfile,'PRES', 'standard_name','sea_water_pressure');
ncwriteatt(outfile,'PRES', 'units','decibar');
ncwriteatt(outfile,'PRES', '_FillValue',single(99999));
ncwriteatt(outfile,'PRES', 'coordinates','TIME DEPTH');
ncwriteatt(outfile,'PRES', 'long_name','pressure of measurement');
ncwriteatt(outfile,'PRES', 'DM_indicator','D');
ncwriteatt(outfile,'PRES', 'axis','Z');
ncwriteatt(outfile,'PRES', 'positive','down');
ncwriteatt(outfile,'PRES', 'reference_scale','mean_sea_level');
ncwriteatt(outfile,'PRES', 'QC_indicator','good data');
ncwriteatt(outfile,'PRES', 'processing_level', mcatinfo.proclvl  );
ncwriteatt(outfile,'PRES', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'PRES', 'valid_max',single(12000));  % CHANGE
ncwriteatt(outfile,'PRES', 'sensor_model',mcatinfo.sensor_model );   %CHANGE 
ncwriteatt(outfile,'PRES', 'sensor_manufacturer',mcatinfo.sensor_manufacturer);   %CHANGE 
ncwriteatt(outfile,'PRES', 'sensor_reference',mcatinfo.sensor_ref);   %CHANGE
ncwriteatt(outfile,'PRES', 'sensor_serial_number',mcatinfo.serial_num);
ncwriteatt(outfile,'PRES', 'sensor_mount',mcatinfo.sensor_mount);
ncwriteatt(outfile,'PRES', 'sensor_orientation',mcatinfo.sens_orientation);
ncwriteatt(outfile,'PRES', 'resolution',mcatinfo.presresolution);	%CHANGE 
ncwriteatt(outfile,'PRES', 'uncertainty',mcatinfo.presuncertainty);   %CHANGE 
ncwriteatt(outfile,'PRES', 'accuracy',mcatinfo.presaccuracy);   %CHANGE 


% ncwriteatt(outfile,'PRES', 'comment','N/A');  %CHANGE




nccreate(outfile,'TEMP', 'Dimensions',{'DEPTH',DepthDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'TEMP', single(temperature));
ncwriteatt(outfile,'TEMP', 'standard_name','sea_water_temperature');
ncwriteatt(outfile,'TEMP', 'units','degree_Celsius');
ncwriteatt(outfile,'TEMP', '_FillValue', single(99999));
ncwriteatt(outfile,'TEMP', 'coordinates','TIME DEPTH');
ncwriteatt(outfile,'TEMP', 'long_name','in situ sea water temperature');
ncwriteatt(outfile,'TEMP', 'DM_indicator','D');
ncwriteatt(outfile,'TEMP', 'reference_scale','ITS-90');
ncwriteatt(outfile,'TEMP', 'QC_indicator','good data');
ncwriteatt(outfile,'TEMP', 'processing_level', mcatinfo.proclvl);
ncwriteatt(outfile,'TEMP', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'TEMP', 'valid_max',single(100));  % CHANGE
ncwriteatt(outfile,'TEMP', 'sensor_model',mcatinfo.sensor_model );   %CHANGE 
ncwriteatt(outfile,'TEMP', 'sensor_manufacturer',mcatinfo.sensor_manufacturer);   %CHANGE 
ncwriteatt(outfile,'TEMP', 'sensor_reference',mcatinfo.sensor_ref);   %CHANGE
ncwriteatt(outfile,'TEMP', 'sensor_serial_number',mcatinfo.serial_num);
ncwriteatt(outfile,'TEMP', 'sensor_mount',mcatinfo.sensor_mount);
ncwriteatt(outfile,'TEMP', 'sensor_orientation',mcatinfo.sens_orientation);
ncwriteatt(outfile,'TEMP', 'resolution',mcatinfo.tempresolution);	%CHANGE 
ncwriteatt(outfile,'TEMP', 'uncertainty',mcatinfo.tempuncertainty);   %CHANGE 
ncwriteatt(outfile,'TEMP', 'accuracy',mcatinfo.tempaccuracy);  
ncwriteatt(outfile,'TEMP', 'precision',mcatinfo.tempprecision);  %CHANGE 




nccreate(outfile,'PSAL', 'Dimensions',{'DEPTH',DepthDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'PSAL', single(salinity));
ncwriteatt(outfile,'PSAL', 'standard_name','sea_water_practical_salinity');
ncwriteatt(outfile,'PSAL', 'units','psu');
ncwriteatt(outfile,'PSAL', '_FillValue', single(99999));
ncwriteatt(outfile,'PSAL', 'coordinates','TIME DEPTH');
ncwriteatt(outfile,'PSAL', 'long_name','sea water salinity');
ncwriteatt(outfile,'PSAL','DM_indicator','D');
ncwriteatt(outfile,'PSAL','reference_scale','PSS-78');
ncwriteatt(outfile,'PSAL', 'QC_indicator','good data');
ncwriteatt(outfile,'PSAL', 'processing_level', mcatinfo.proclvl);
ncwriteatt(outfile,'PSAL', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'PSAL', 'valid_max',single(40));  % CHANGE
ncwriteatt(outfile,'PSAL', 'sensor_model',mcatinfo.sensor_model );   %CHANGE 
ncwriteatt(outfile,'PSAL', 'sensor_manufacturer',mcatinfo.sensor_manufacturer);   %CHANGE 
ncwriteatt(outfile,'PSAL', 'sensor_reference',mcatinfo.sensor_ref);   %CHANGE
ncwriteatt(outfile,'PSAL', 'sensor_serial_number',mcatinfo.serial_num);
ncwriteatt(outfile,'PSAL', 'sensor_mount',mcatinfo.sensor_mount);
ncwriteatt(outfile,'PSAL', 'sensor_orientation',mcatinfo.sens_orientation);
ncwriteatt(outfile,'PSAL', 'resolution',mcatinfo.salresolution);	%CHANGE 
ncwriteatt(outfile,'PSAL', 'uncertainty',mcatinfo.saluncertainty);   %CHANGE 
ncwriteatt(outfile,'PSAL', 'accuracy',mcatinfo.salaccuracy);  
ncwriteatt(outfile,'PSAL', 'precision',mcatinfo.salprecision);  %CHANGE 



nccreate(outfile,'CNDC', 'Dimensions',{'DEPTH',DepthDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'CNDC', single(conductivity));
ncwriteatt(outfile,'CNDC', 'standard_name','sea_water_electrical_conductivity');
ncwriteatt(outfile,'CNDC','units','mS/cm');
ncwriteatt(outfile,'CNDC', '_FillValue', single(99999));
ncwriteatt(outfile,'CNDC', 'coordinates','TIME DEPTH');
ncwriteatt(outfile,'CNDC', 'long_name','sea water conductivity');
ncwriteatt(outfile,'CNDC', 'DM_indicator','D');
ncwriteatt(outfile,'CNDC', 'reference_scale',' ');
ncwriteatt(outfile,'CNDC', 'QC_indicator','good data');
ncwriteatt(outfile,'CNDC', 'processing_level', mcatinfo.proclvl);
ncwriteatt(outfile,'CNDC', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'CNDC', 'valid_max',single(7));  % CHANGE
ncwriteatt(outfile,'CNDC', 'sensor_model',mcatinfo.sensor_model );   %CHANGE 
ncwriteatt(outfile,'CNDC', 'sensor_manufacturer',mcatinfo.sensor_manufacturer);   %CHANGE 
ncwriteatt(outfile,'CNDC', 'sensor_reference',mcatinfo.sensor_ref);   %CHANGE
ncwriteatt(outfile,'CNDC', 'sensor_serial_number',mcatinfo.serial_num);
ncwriteatt(outfile,'CNDC', 'sensor_mount',mcatinfo.sensor_mount);
ncwriteatt(outfile,'CNDC', 'sensor_orientation',mcatinfo.sens_orientation);
ncwriteatt(outfile,'CNDC', 'resolution',mcatinfo.condresolution);	%CHANGE 
ncwriteatt(outfile,'CNDC', 'uncertainty',mcatinfo.conduncertainty);   %CHANGE 
ncwriteatt(outfile,'CNDC', 'accuracy',mcatinfo.condaccuracy);  
ncwriteatt(outfile,'CNDC', 'precision',mcatinfo.condprecision);  %CHANGE 




% ------------------+
%   End             |
% ------------------+
fprintf('Finish %s\n', datestr(now));


end
