%% GATHER AND CONVERT Stacked ROCKALL TROUGH MOORING DATA IN .MAT FORMAT AND CONVERT TO NETCDF
% Lewis Drysdale, SAMS, 2020
% close('all')
% --------
% Update: Kristin Burmeister, SAMS, Apr 2024 - rewrite for export of
% stacked UV and TS files


%% filename 
startdate = '201407'; % does not change!
enddate='202407';
version='v1';
filename =strcat('Rockall_Trough_mooring_stacked_TSUV_',startdate,'_',enddate,'_', version)
outdir                      = [basedir '/osnap/data/moor/THREDDS'];
if exist(outdir,'dir')==0;mkdir(outdir);end
outfile                    = fullfile(outdir, [filename '.nc']);

%% Input data
tsdir                       = [basedir '/osnap/data/moor/proc/hydro_grid_merged/'];
vldir                       = [basedir '/osnap/data/moor/proc/velocity_grid_merged/'];%[pathosnap '/data/moor/proc/velocity_grid_merged/'];

%% T S DATA 
% Load western array data
TSWB = load([tsdir 'Hydro_rtwb_stacked_201407_202407_v0.mat']);

% rename var and make NaN 99999
TS_WEST                 =TSWB.Tstacked;
SS_WEST                 =TSWB.Sstacked;
PS_WEST_TS                 =TSWB.Pstacked;
ZS_WEST_TS                 =TSWB.z_stacked;
time                    =TSWB.JG; 
    
%Load eastern array data
TSEB = load([tsdir 'Hydro_rteb_stacked_201407_202407_v0.mat']);
TS_EAST                 =TSEB.Tstacked;
SS_EAST                 =TSEB.Sstacked;
PS_EAST_TS                 =TSEB.Pstacked;
ZS_EAST_TS                =TSEB.z_stacked;

%  VELOCITY
% western boundary 1
ffile =['CM_rtwb1_stacked_',startdate,'_',enddate,'_', version '.mat'];
UVWB1 = load([vldir ffile]);
% rename vars
US_WEST_1             =UVWB1.Ustacked; 
VS_WEST_1             =UVWB1.Vstacked; 
PS_WEST_1_UV          =UVWB1.Pstacked;
ZS_WEST_1_UV          =UVWB1.z_stacked;

% western boundary 2
ffile =['CM_rtwb2_stacked_',startdate,'_',enddate,'_', version '.mat'];
UVWB2 = load([vldir ffile]);
% rename vars
US_WEST_2             =UVWB2.Ustacked; 
VS_WEST_2             =UVWB2.Vstacked; 
PS_WEST_2_UV          =UVWB2.Pstacked;
ZS_WEST_2_UV          =UVWB2.z_stacked;
            
% eastern boundary 
ffile =['CM_rteb_stacked_',startdate,'_',enddate,'_', version '.mat'];
UVEB1 = load([vldir ffile]);
% rename vars
US_EAST_1             =UVEB1.Ustacked; 
VS_EAST_1             =UVEB1.Vstacked; 
PS_EAST_1_UV          =UVEB1.Pstacked;
ZS_EAST_1_UV          =UVEB1.z_stacked;
            
% WRITE FILE 

fprintf('Begin %s\n', datestr(now));

% Dimensions
TimeDim                     = length(time);

% Open the file to write
nc                          = netcdf.create(outfile,'CLOBBER'); % CLOBBER overwrites 
netcdf.close(nc)

  
%  Write global attributes 
% These use the CF (Climate and Forecast) metadata conventions <http:// https://cfconventions.org/index.html  
% ><https://cfconventions.org/ https://cfconventions.org/>


ncwriteatt(outfile,'/','title', ['Stacked Rockall Trough mooring data ',datestr(time(1),'mm/yyyy'),'-',datestr(time(end),'mm/yyyy')]); 
ncwriteatt(outfile,'/','institution', ['Scottish Association for Marine Science, Scottish Marine Institute Oban, Argyll, PA37 1QA, UK']); 
ncwriteatt(outfile,'/','history', 'Stacked water temperature, salinity, and velocity data at 12-hour time intervals at depth as measured.'); 
ncwriteatt(outfile,'/','id', filename);   
ncwriteatt(outfile,'/','source', 'subsurface mooring');
ncwriteatt(outfile,'/','project','Climate Linked Atlantic Sector Science, AtlantiS, OSNAP');

% GEO-SPATIAL-TEMPORAL
ncwriteatt(outfile,'/','area', 'North Atlantic Ocean');
ncwriteatt(outfile,'/','geospatial_vertical_min', ZS_WEST_TS(1));
ncwriteatt(outfile,'/','geospatial_vertical_max', ZS_WEST_TS(end));
ncwriteatt(outfile,'/','geospatial_vertical_positive', 'down');
ncwriteatt(outfile,'/','geospatial_vertical_units', 'm');
ncwriteatt(outfile,'/','time_coverage_start', datestr(min(time),'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','time_coverage_end', datestr(max(time),'yyyy-mm-ddTHH:MM:SSZ'));

dDD = floor(time(end) - time(1));
dHH = floor((time(end) - time(1) - dDD)*24);
ncwriteatt(outfile,'/','time_coverage_duration', ['P',num2str(dDD),'D',num2str(dHH),'H']);
dMM = round(mean(diff(time)*24*60));
ncwriteatt(outfile,'/','time_coverage_resolution', ['PT',num2str(dMM),'M']);  

% CONVENTIONS USED 
ncwriteatt(outfile,'/','netcdf_version','4.3');

% ACKNOWLEDGMENT
ncwriteatt(outfile,'/','publisher_url', '');  %CHANGE  
ncwriteatt(outfile,'/','references', 'http://www.o-snap.org');
ncwriteatt(outfile,'/','citation', 'These data were collected and made freely available by the OSNAP project and the national programs that contribute to it.');
ncwriteatt(outfile,'/','acknowledgement', 'Funding source: the UK Natural Environment Research Council (NERC), UK OSNAP project');  %CHANGE


% PROVENANCE
ncwriteatt(outfile,'/','date_created', datestr(now+8/24,'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','date_modified', datestr(now+8/24,'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','history', 'Delayed time processed quality controlled');
ncwriteatt(outfile,'/','processing_level','');
ncwriteatt(outfile,'/','QC_indicator','');
   
% Write coordinate variables
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

nccreate(outfile,'ZS_WEST_TS', 'Dimensions',{'ZS_WEST_TS',length(ZS_WEST_TS)}, 'Datatype','single');
ncwrite(outfile,'ZS_WEST_TS', single(ZS_WEST_TS));
ncwriteatt(outfile,'ZS_WEST_TS', 'standard_name','sea_water_depth');
ncwriteatt(outfile,'ZS_WEST_TS', 'units','m');
ncwriteatt(outfile,'ZS_WEST_TS', '_FillValue','NaN');
ncwriteatt(outfile,'ZS_WEST_TS', 'coordinates','ZS_WEST_TS');
ncwriteatt(outfile,'ZS_WEST_TS', 'long_name','depth grid');
ncwriteatt(outfile,'ZS_WEST_TS', 'axis','Z');
ncwriteatt(outfile,'ZS_WEST_TS', 'positive','down');
ncwriteatt(outfile,'ZS_WEST_TS', 'QC_indicator','good data');

nccreate(outfile,'ZS_EAST_TS', 'Dimensions',{'ZS_EAST_TS',length(ZS_EAST_TS)}, 'Datatype','single');
ncwrite(outfile,'ZS_EAST_TS', single(ZS_EAST_TS));
ncwriteatt(outfile,'ZS_EAST_TS', 'standard_name','sea_water_depth');
ncwriteatt(outfile,'ZS_EAST_TS', 'units','m');
ncwriteatt(outfile,'ZS_EAST_TS', '_FillValue','NaN');
ncwriteatt(outfile,'ZS_EAST_TS', 'coordinates','ZS_EAST_TS');
ncwriteatt(outfile,'ZS_EAST_TS', 'long_name','depth grid');
ncwriteatt(outfile,'ZS_EAST_TS', 'axis','Z');
ncwriteatt(outfile,'ZS_EAST_TS', 'positive','down');
ncwriteatt(outfile,'ZS_EAST_TS', 'QC_indicator','good data');

nccreate(outfile,'ZS_WEST_1_UV', 'Dimensions',{'ZS_WEST_1_UV',length(ZS_WEST_1_UV)}, 'Datatype','single');
ncwrite(outfile,'ZS_WEST_1_UV', single(ZS_WEST_1_UV));
ncwriteatt(outfile,'ZS_WEST_1_UV', 'standard_name','sea_water_depth');
ncwriteatt(outfile,'ZS_WEST_1_UV', 'units','m');
ncwriteatt(outfile,'ZS_WEST_1_UV', '_FillValue','NaN');
ncwriteatt(outfile,'ZS_WEST_1_UV', 'coordinates','ZS_WEST_1_UV');
ncwriteatt(outfile,'ZS_WEST_1_UV', 'long_name','depth grid');
ncwriteatt(outfile,'ZS_WEST_1_UV', 'axis','Z');
ncwriteatt(outfile,'ZS_WEST_1_UV', 'positive','down');
ncwriteatt(outfile,'ZS_WEST_1_UV', 'QC_indicator','good data');

nccreate(outfile,'ZS_WEST_2_UV', 'Dimensions',{'ZS_WEST_2_UV',length(ZS_WEST_2_UV)}, 'Datatype','single');
ncwrite(outfile,'ZS_WEST_2_UV', single(ZS_WEST_2_UV));
ncwriteatt(outfile,'ZS_WEST_2_UV', 'standard_name','sea_water_depth');
ncwriteatt(outfile,'ZS_WEST_2_UV', 'units','m');
ncwriteatt(outfile,'ZS_WEST_2_UV', '_FillValue','NaN');
ncwriteatt(outfile,'ZS_WEST_2_UV', 'coordinates','ZS_WEST_2_UV');
ncwriteatt(outfile,'ZS_WEST_2_UV', 'long_name','depth grid');
ncwriteatt(outfile,'ZS_WEST_2_UV', 'axis','Z');
ncwriteatt(outfile,'ZS_WEST_2_UV', 'positive','down');
ncwriteatt(outfile,'ZS_WEST_2_UV', 'QC_indicator','good data');


nccreate(outfile,'ZS_EAST_1_UV', 'Dimensions',{'ZS_EAST_1_UV',length(ZS_EAST_1_UV)}, 'Datatype','single');
ncwrite(outfile,'ZS_EAST_1_UV', single(ZS_EAST_1_UV));
ncwriteatt(outfile,'ZS_EAST_1_UV', 'standard_name','sea_water_depth');
ncwriteatt(outfile,'ZS_EAST_1_UV', 'units','m');
ncwriteatt(outfile,'ZS_EAST_1_UV', '_FillValue','NaN');
ncwriteatt(outfile,'ZS_EAST_1_UV', 'coordinates','ZS_EAST_1_UV');
ncwriteatt(outfile,'ZS_EAST_1_UV', 'long_name','depth grid');
ncwriteatt(outfile,'ZS_EAST_1_UV', 'axis','Z');
ncwriteatt(outfile,'ZS_EAST_1_UV', 'positive','down');
ncwriteatt(outfile,'ZS_EAST_1_UV', 'QC_indicator','good data');

% Write data variables
nccreate(outfile,'TS_EAST', 'Dimensions',{'ZS_EAST_TS',length(ZS_EAST_TS), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'TS_EAST', single(TS_EAST));
ncwriteatt(outfile,'TS_EAST', 'standard_name','sea_water_conservative_temperature');
ncwriteatt(outfile,'TS_EAST', 'units','^{o} C');
ncwriteatt(outfile,'TS_EAST', '_FillValue', 'NaN');
ncwriteatt(outfile,'TS_EAST', 'coordinates','TIME ZS_EAST_TS');
ncwriteatt(outfile,'TS_EAST', 'long_name','water temperature at eastern boundary 57.1N/-9.6W');
ncwriteatt(outfile,'TS_EAST', 'reference_scale','ITS-90');
ncwriteatt(outfile,'TS_EAST', 'QC_indicator','good data');
ncwriteatt(outfile,'TS_EAST', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'TS_EAST', 'valid_max',single(100));  % CHANGE

nccreate(outfile,'TS_WEST', 'Dimensions',{'ZS_WEST_TS',length(ZS_WEST_TS), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'TS_WEST', single(TS_WEST));
ncwriteatt(outfile,'TS_WEST', 'standard_name','sea_water_conservative_temperature');
ncwriteatt(outfile,'TS_WEST', 'units','^{o} C');
ncwriteatt(outfile,'TS_WEST', '_FillValue', 'NaN');
ncwriteatt(outfile,'TS_WEST', 'coordinates','TIME ZS_WEST_TS');
ncwriteatt(outfile,'TS_WEST', 'long_name','water temperature at western boundary 57.5N/-12.3W');
ncwriteatt(outfile,'TS_WEST', 'reference_scale','ITS-90');
ncwriteatt(outfile,'TS_WEST', 'QC_indicator','good data');
ncwriteatt(outfile,'TS_WEST', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'TS_WEST', 'valid_max',single(100));  % CHANGE

nccreate(outfile,'SS_EAST', 'Dimensions',{'ZS_EAST_TS',length(ZS_EAST_TS), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'SS_EAST', single(SS_EAST));
ncwriteatt(outfile,'SS_EAST', 'standard_name','sea_water_absolute_salinity');
ncwriteatt(outfile,'SS_EAST', 'units','g/kg');
ncwriteatt(outfile,'SS_EAST', '_FillValue', 'NaN');
ncwriteatt(outfile,'SS_EAST', 'coordinates','TIME ZS_EAST_TS');
ncwriteatt(outfile,'SS_EAST', 'long_name','sea water salinity at eastern boundary 57.1N/-9.6W');
ncwriteatt(outfile,'SS_EAST','reference_scale','TEOS-10');
ncwriteatt(outfile,'SS_EAST', 'QC_indicator','good data');
ncwriteatt(outfile,'SS_EAST', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'SS_EAST', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'PS_EAST_TS', 'Dimensions',{'ZS_EAST_TS',length(ZS_EAST_TS), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'PS_EAST_TS', single(PS_EAST_TS));
ncwriteatt(outfile,'PS_EAST_TS', 'standard_name','sea_water_pressure');
ncwriteatt(outfile,'PS_EAST_TS', 'units','decibar');
ncwriteatt(outfile,'PS_EAST_TS', '_FillValue', 'NaN');
ncwriteatt(outfile,'PS_EAST_TS', 'coordinates','TIME ZS_EAST_TS');
ncwriteatt(outfile,'PS_EAST_TS', 'long_name','sea water pressure');
ncwriteatt(outfile,'PS_EAST_TS', 'QC_indicator','good data');

nccreate(outfile,'PS_WEST_TS', 'Dimensions',{'ZS_WEST_TS',length(ZS_WEST_TS), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'PS_WEST_TS', single(PS_WEST_TS));
ncwriteatt(outfile,'PS_WEST_TS', 'standard_name','sea_water_pressure');
ncwriteatt(outfile,'PS_WEST_TS', 'units','decibar');
ncwriteatt(outfile,'PS_WEST_TS', '_FillValue', 'NaN');
ncwriteatt(outfile,'PS_WEST_TS', 'coordinates','TIME ZS_WEST_TS');
ncwriteatt(outfile,'PS_WEST_TS', 'long_name','sea water pressure');
ncwriteatt(outfile,'PS_WEST_TS', 'QC_indicator','good data');

nccreate(outfile,'SS_WEST', 'Dimensions',{'ZS_WEST_TS',length(ZS_WEST_TS), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'SS_WEST', single(SS_WEST));
ncwriteatt(outfile,'SS_WEST', 'standard_name','sea_water_absolute_salinity');
ncwriteatt(outfile,'SS_WEST', 'units','g/kg');
ncwriteatt(outfile,'SS_WEST', '_FillValue', 'NaN');
ncwriteatt(outfile,'SS_WEST', 'coordinates','TIME ZS_WEST_TS');
ncwriteatt(outfile,'SS_WEST', 'long_name','sea water salinity at western boundary 57.5N/-12.3W');
ncwriteatt(outfile,'SS_WEST','reference_scale','TEOS-10');
ncwriteatt(outfile,'SS_WEST', 'QC_indicator','good data');
ncwriteatt(outfile,'SS_WEST', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'SS_WEST', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'US_WEST_1', 'Dimensions',{'ZS_WEST_1_UV',length(ZS_WEST_1_UV), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'US_WEST_1', single(US_WEST_1));
ncwriteatt(outfile,'US_WEST_1', 'standard_name','velocity');
ncwriteatt(outfile,'US_WEST_1', 'units','cm/s');
ncwriteatt(outfile,'US_WEST_1', '_FillValue', 'NaN');
ncwriteatt(outfile,'US_WEST_1', 'coordinates','TIME ZS_WEST_1_UV');
ncwriteatt(outfile,'US_WEST_1', 'long_name','current speed u-direction at western boundary 1 XXN/XXW');
ncwriteatt(outfile,'US_WEST_1', 'QC_indicator','good data');
% ncwriteatt(outfile,'US_WEST_1', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'US_WEST_1', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'VS_WEST_1', 'Dimensions',{'ZS_WEST_1_UV',length(ZS_WEST_1_UV), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'VS_WEST_1', single(VS_WEST_1));
ncwriteatt(outfile,'VS_WEST_1', 'standard_name','velocity');
ncwriteatt(outfile,'VS_WEST_1', 'units','cm/s');
ncwriteatt(outfile,'VS_WEST_1', '_FillValue', 'NaN');
ncwriteatt(outfile,'VS_WEST_1', 'coordinates','TIME ZS_WEST_1_UV');
ncwriteatt(outfile,'VS_WEST_1', 'long_name','current speed v-direction at western boundary 1 XXN/XXW');
ncwriteatt(outfile,'VS_WEST_1', 'QC_indicator','good data');
% ncwriteatt(outfile,'U_WEST_1', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'U_WEST_1', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'PS_WEST_1_UV', 'Dimensions',{'ZS_WEST_1_UV',length(ZS_WEST_1_UV), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'PS_WEST_1_UV', single(PS_WEST_1_UV));
ncwriteatt(outfile,'PS_WEST_1_UV', 'standard_name','sea_water_pressure');
ncwriteatt(outfile,'PS_WEST_1_UV', 'units','decibar');
ncwriteatt(outfile,'PS_WEST_1_UV', '_FillValue', 'NaN');
ncwriteatt(outfile,'PS_WEST_1_UV', 'coordinates','TIME ZS_WEST_1_UV');
ncwriteatt(outfile,'PS_WEST_1_UV', 'long_name','sea water pressure');
ncwriteatt(outfile,'PS_WEST_1_UV', 'QC_indicator','good data');

nccreate(outfile,'US_WEST_2', 'Dimensions',{'ZS_WEST_2_UV',length(ZS_WEST_2_UV), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'US_WEST_2', single(US_WEST_2));
ncwriteatt(outfile,'US_WEST_2', 'standard_name','velocity');
ncwriteatt(outfile,'US_WEST_2', 'units','cm/s');
ncwriteatt(outfile,'US_WEST_2', '_FillValue', 'NaN');
ncwriteatt(outfile,'US_WEST_2', 'coordinates','TIME ZS_WEST_2_UV');
ncwriteatt(outfile,'US_WEST_2', 'long_name','current speed u-direction at western boundary 2 XXN/XXW');
ncwriteatt(outfile,'US_WEST_2', 'QC_indicator','good data');
% ncwriteatt(outfile,'US_WEST_2', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'US_WEST_2', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'VS_WEST_2', 'Dimensions',{'ZS_WEST_2_UV',length(ZS_WEST_2_UV), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'VS_WEST_2', single(VS_WEST_2));
ncwriteatt(outfile,'VS_WEST_2', 'standard_name','velocity');
ncwriteatt(outfile,'VS_WEST_2', 'units','cm/s');
ncwriteatt(outfile,'VS_WEST_2', '_FillValue', 'NaN');
ncwriteatt(outfile,'VS_WEST_2', 'coordinates','TIME ZS_WEST_2_UV');
ncwriteatt(outfile,'VS_WEST_2', 'long_name','current speed v-direction at western boundary 2 XXN/XXW');
ncwriteatt(outfile,'VS_WEST_2', 'QC_indicator','good data');
% ncwriteatt(outfile,'US_WEST_2', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'US_WEST_2', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'PS_WEST_2_UV', 'Dimensions',{'ZS_WEST_2_UV',length(ZS_WEST_2_UV), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'PS_WEST_2_UV', single(PS_WEST_2_UV));
ncwriteatt(outfile,'PS_WEST_2_UV', 'standard_name','sea_water_pressure');
ncwriteatt(outfile,'PS_WEST_2_UV', 'units','decibar');
ncwriteatt(outfile,'PS_WEST_2_UV', '_FillValue', 'NaN');
ncwriteatt(outfile,'PS_WEST_2_UV', 'coordinates','TIME ZS_WEST_2_UV');
ncwriteatt(outfile,'PS_WEST_2_UV', 'long_name','sea water pressure');
ncwriteatt(outfile,'PS_WEST_2_UV', 'QC_indicator','good data');

nccreate(outfile,'US_EAST_1', 'Dimensions',{'ZS_EAST_1_UV',length(ZS_EAST_1_UV), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'US_EAST_1', single(US_EAST_1));
ncwriteatt(outfile,'US_EAST_1', 'standard_name','velocity');
ncwriteatt(outfile,'US_EAST_1', 'units','cm/s');
ncwriteatt(outfile,'US_EAST_1', '_FillValue', 'NaN');
ncwriteatt(outfile,'US_EAST_1', 'coordinates','TIME ZS_EAST_1_UV');
ncwriteatt(outfile,'US_EAST_1', 'long_name','current speed u-direction at eastern boundary 57.1N/-9.6W');
ncwriteatt(outfile,'US_EAST_1', 'QC_indicator','good data');
% ncwriteatt(outfile,'US_EAST_1', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'US_EAST_1', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'VS_EAST_1', 'Dimensions',{'ZS_EAST_1_UV',length(ZS_EAST_1_UV), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'VS_EAST_1', single(VS_EAST_1));
ncwriteatt(outfile,'VS_EAST_1', 'standard_name','velocity');
ncwriteatt(outfile,'VS_EAST_1', 'units','cm/s');
ncwriteatt(outfile,'VS_EAST_1', '_FillValue', 'NaN');
ncwriteatt(outfile,'VS_EAST_1', 'coordinates','TIME ZS_EAST_1_UV');
ncwriteatt(outfile,'VS_EAST_1', 'long_name','current speed v-direction at eastern boundary 57.1N/-9.6W');
ncwriteatt(outfile,'VS_EAST_1', 'QC_indicator','good data');
% ncwriteatt(outfile,'US_EAST_1', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'US_EAST_1', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'PS_EAST_1_UV', 'Dimensions',{'ZS_EAST_1_UV',length(ZS_EAST_1_UV), 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'PS_EAST_1_UV', single(PS_EAST_1_UV));
ncwriteatt(outfile,'PS_EAST_1_UV', 'standard_name','sea_water_pressure');
ncwriteatt(outfile,'PS_EAST_1_UV', 'units','decibar');
ncwriteatt(outfile,'PS_EAST_1_UV', '_FillValue', 'NaN');
ncwriteatt(outfile,'PS_EAST_1_UV', 'coordinates','TIME ZS_EAST_1_UV');
ncwriteatt(outfile,'PS_EAST_1_UV', 'long_name','sea water pressure');
ncwriteatt(outfile,'PS_EAST_1_UV', 'QC_indicator','good data');

fprintf('Finish %s\n', datestr(now));

clearvars -except outfile
finfo2=ncinfo(outfile);