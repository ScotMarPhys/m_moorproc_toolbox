% GATHER AND CONVERT ROCKALL TROUGH MOORING DATA IN .MAT FORMAT AND CONVERT
% TO NETCDF
%
% -------------------------------------------------------------------------
% Adapted by Lewis Drysdale SAMS, 2020 for Rockall Trough gridded mooring 
% data from a function write_MCAT_to_NetCDF.m with history:
%
% Copyright, IFREMER, 2012
% Modified by Thomas Richter at NIOZ, 2014
% Adapted for OSNAP MCTD data by Feili Li at Duke University, 2016
% Last updated: August 12, 2016
% Edited by Loic Houpert, SAMS 18/11/2016
%
% -------------------------------------------------------------------------
clearvars; close('all')

last_year = '_2020';

%%%%%%%%%%%%%%%%%%%%%%%%%% T S DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load western array data
indir = 'D:\Work_computer_sync\OSNAP_work\osnap\data\moor\proc\hydro_grid_merged\'; % SJ
%indir               =[getenv('OSNAP') '/data/moor/proc/hydro_grid_merged/'];
ffile               ='RTWB_merg_linear_interp_2018.mat';
load([indir ffile]);
% rename vars, convert degC to K, and make NaN 99999
TG_WEST             =RTWB_merg.TGfs2;
% convert to kelvin
TG_WEST             = convtemp(TG_WEST,'C','K');
TG_WEST(isnan(TG_WEST))=99999;

SG_WEST             =RTWB_merg.SGfs2; SG_WEST(isnan(SG_WEST))=99999;
pressure            =RTWB_merg.PGfs(:,1); 
time                =RTWB_merg.JG; 
clearvars -except TG_WEST SG_WEST pressure time indir
    
%Load eastern array data
ffile               ='RTEB_merg_linear_interp_2018.mat';
load([indir ffile]);
% rename vars
TG_EAST             =RTEB_merg.TGfs2; 
% convert to kelvin
TG_EAST             = convtemp(TG_EAST,'C','K');
TG_EAST(isnan(TG_EAST))=99999;
SG_EAST             =RTEB_merg.SGfs2; SG_EAST(isnan(SG_EAST))=99999;
clearvars -except TG_WEST SG_WEST TG_EAST SG_EAST pressure time indir

%%%%%%%%%%%%%%%%%%%%%%%%%%% VELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indir = 'D:\Work_computer_sync\OSNAP_work\osnap\data\moor\proc\velocity_grid_merged\'; % SJ
% indir               =[getenv('OSNAP') '/data/moor/proc/velocity_grid_merged/'];

% western boundary 1
ffile               ='RTWB1_merg_linear_interp_2020.mat';
load([indir ffile]);
% rename vars
U_WEST_1             =RTWB1_merg_CM.UGfs2; U_WEST_1(isnan(U_WEST_1))=99999;
V_WEST_1             =RTWB1_merg_CM.VGfs2; V_WEST_1(isnan(V_WEST_1))=99999;
W_WEST_1             =RTWB1_merg_CM.WGfs2; W_WEST_1(isnan(W_WEST_1))=99999;
clearvars -except U_WEST_1 V_WEST_1 W_WEST_1...
                TG_WEST SG_WEST...
                TG_EAST SG_EAST...
                pressure time indir

% western boundary 2
ffile               ='RTWB2_merg_linear_interp_2020.mat';
load([indir ffile]);
% rename vars
U_WEST_2             =RTWB2_merg_CM.UGfs2;U_WEST_2(isnan(U_WEST_2))=99999;
V_WEST_2             =RTWB2_merg_CM.VGfs2;V_WEST_2(isnan(V_WEST_2))=99999;
W_WEST_2             =RTWB2_merg_CM.WGfs2;W_WEST_2(isnan(W_WEST_2))=99999;
clearvars -except U_WEST_2 V_WEST_2 W_WEST_2...
                U_WEST_1 V_WEST_1 W_WEST_1...
                TG_WEST SG_WEST...
                TG_EAST SG_EAST...
                pressure time indir
            
% eastern boundary 
ffile               ='RTEB_merg_linear_interp_2020.mat';
load([indir ffile]);
% rename vars
U_EAST             =RTEB_merg_CM.UGfs2;U_EAST(isnan(U_EAST))=99999;
V_EAST             =RTEB_merg_CM.VGfs2;U_EAST(isnan(U_EAST))=99999;
W_EAST             =RTEB_merg_CM.WGfs2;U_EAST(isnan(U_EAST))=99999;
clearvars -except U_EAST V_EAST W_EAST...
                U_WEST_2 V_WEST_2 W_WEST_2...
                U_WEST_1 V_WEST_1 W_WEST_1...
                TG_WEST SG_WEST...
                TG_EAST SG_EAST...
                pressure time
            
            
%%%%%%%%%%%%%%%%%%%%% WRITE FILE TO DIR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdir                      = ['D:\Work_computer_sync\OSNAP_work\osnap\data\moor\THREDDS\']; % SJ
% outdir                      = [getenv('OSNAP') '/data/moor/THREDDS/'];

fprintf('Begin %s\n', datestr(now));

% Dimensions
TimeDim                     = length(time);
PressureDim                 = length(pressure);

% filename 
filename                    ='Rockall-Trough-Mooring-Time-Series-2020';

% Open the file to write
outfile                     = [outdir filename,'.nc'];
nc                          = netcdf.create(outfile,'CLOBBER'); % CLOBBER overwrites 
netcdf.close(nc)

  
%%%%%%%%%%%%%%%%%%%% Write global attributes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% these use the CF (Climate and Forecast) metadata conventions 
% https://cfconventions.org/index.html

ncwriteatt(outfile,'/','title', ['CLASS Rockall Trough mooring data ',datestr(time(1),'mm/yyyy'),'-',datestr(time(end),'mm/yyyy')]); 
ncwriteatt(outfile,'/','institution', ['Scottish Association for Marine Science, Scottish Marine Institute Oban, Argyll, PA37 1QA, UK']); 
ncwriteatt(outfile,'/','history', 'Gridded water temperature, salinity, and velocity data at 20db pressure and 12-hour time intervals.'); 
ncwriteatt(outfile,'/','id', filename);   
ncwriteatt(outfile,'/','source', 'subsurface mooring');
ncwriteatt(outfile,'/','project','Climate Linked Atlantic Sector Science');

%%%%%%%%%%%%%%%%%%%%%%% GEO-SPATIAL-TEMPORAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncwriteatt(outfile,'/','area', 'North Atlantic Ocean');
ncwriteatt(outfile,'/','geospatial_vertical_min', pressure(1));
ncwriteatt(outfile,'/','geospatial_vertical_max', pressure(end));
ncwriteatt(outfile,'/','geospatial_vertical_positive', 'down');
ncwriteatt(outfile,'/','geospatial_vertical_units', 'decibar');
ncwriteatt(outfile,'/','time_coverage_start', datestr(min(time),'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','time_coverage_end', datestr(max(time),'yyyy-mm-ddTHH:MM:SSZ'));

dDD = floor(time(end) - time(1));
dHH = floor((time(end) - time(1) - dDD)*24);
ncwriteatt(outfile,'/','time_coverage_duration', ['P',num2str(dDD),'D',num2str(dHH),'H']);
dMM = round(mean(diff(time)*24*60));
ncwriteatt(outfile,'/','time_coverage_resolution', ['PT',num2str(dMM),'M']);  

%%%%%%%%%%%%%%%%%%%%% CONVENTIONS USED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncwriteatt(outfile,'/','netcdf_version','4.3');

%%%%%%%%%%%%%%%%%%%%% ACKNOWLEDGMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncwriteatt(outfile,'/','publisher_url', '');  %CHANGE  
ncwriteatt(outfile,'/','references', 'http://www.o-snap.org');
ncwriteatt(outfile,'/','citation', 'These data were collected and made freely available by the OSNAP project and the national programs that contribute to it.');
ncwriteatt(outfile,'/','acknowledgement', 'Funding source: the UK Natural Environment Research Council (NERC), UK OSNAP project');  %CHANGE


%%%%%%%%%%%%%%%%%%%%% PROVENANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncwriteatt(outfile,'/','date_created', datestr(now+8/24,'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','date_modified', datestr(now+8/24,'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(outfile,'/','history', 'Delayed time processed quality controlled');
ncwriteatt(outfile,'/','processing_level','');
ncwriteatt(outfile,'/','QC_indicator','');
   
%%%%%%%%%%%%%%%%%%%%% Write coordinate variables %%%%%%%%%%%%%%%%%%%%%%%%%%

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


nccreate(outfile,'PRES', 'Dimensions',{'PRES',PressureDim}, 'Datatype','single');
ncwrite(outfile,'PRES', single(pressure));
ncwriteatt(outfile,'PRES', 'standard_name','sea_water_pressure');
ncwriteatt(outfile,'PRES', 'units','decibar');
ncwriteatt(outfile,'PRES', '_FillValue',single(99999));
ncwriteatt(outfile,'PRES', 'coordinates','PRES');
ncwriteatt(outfile,'PRES', 'long_name','pressure grid');
ncwriteatt(outfile,'PRES', 'axis','Z');
ncwriteatt(outfile,'PRES', 'positive','down');
ncwriteatt(outfile,'PRES', 'QC_indicator','good data');

%%%%%%%%%%%%%%%%%%%%% Write data variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

nccreate(outfile,'TG_EAST', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'TG_EAST', single(TG_EAST));
ncwriteatt(outfile,'TG_EAST', 'standard_name','sea_water_conservative_temperature');
ncwriteatt(outfile,'TG_EAST', 'units','K');
ncwriteatt(outfile,'TG_EAST', '_FillValue', single(99999));
ncwriteatt(outfile,'TG_EAST', 'coordinates','TIME PRES');
ncwriteatt(outfile,'TG_EAST', 'long_name','water temperature at eastern boundary 57.1N/-9.6W');
ncwriteatt(outfile,'TG_EAST', 'reference_scale','ITS-90');
ncwriteatt(outfile,'TG_EAST', 'QC_indicator','good data');
ncwriteatt(outfile,'TG_EAST', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'TG_EAST', 'valid_max',single(100));  % CHANGE

nccreate(outfile,'TG_WEST', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'TG_WEST', single(TG_WEST));
ncwriteatt(outfile,'TG_WEST', 'standard_name','sea_water_conservative_temperature');
ncwriteatt(outfile,'TG_WEST', 'units','K');
ncwriteatt(outfile,'TG_WEST', '_FillValue', single(99999));
ncwriteatt(outfile,'TG_WEST', 'coordinates','TIME PRES');
ncwriteatt(outfile,'TG_WEST', 'long_name','water temperature at western boundary 57.5N/-12.3W');
ncwriteatt(outfile,'TG_WEST', 'reference_scale','ITS-90');
ncwriteatt(outfile,'TG_WEST', 'QC_indicator','good data');
ncwriteatt(outfile,'TG_WEST', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'TG_WEST', 'valid_max',single(100));  % CHANGE

nccreate(outfile,'SG_EAST', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'SG_EAST', single(SG_EAST));
ncwriteatt(outfile,'SG_EAST', 'standard_name','sea_water_absolute_salinity');
ncwriteatt(outfile,'SG_EAST', 'units','g/kg');
ncwriteatt(outfile,'SG_EAST', '_FillValue', single(99999));
ncwriteatt(outfile,'SG_EAST', 'coordinates','TIME PRES');
ncwriteatt(outfile,'SG_EAST', 'long_name','sea water salinity at eastern boundary 57.1N/-9.6W');
ncwriteatt(outfile,'SG_EAST','reference_scale','TEOS-10');
ncwriteatt(outfile,'SG_EAST', 'QC_indicator','good data');
ncwriteatt(outfile,'SG_EAST', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'SG_EAST', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'SG_WEST', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'SG_WEST', single(SG_WEST));
ncwriteatt(outfile,'SG_WEST', 'standard_name','sea_water_absolute_salinity');
ncwriteatt(outfile,'SG_WEST', 'units','g/kg');
ncwriteatt(outfile,'SG_WEST', '_FillValue', single(99999));
ncwriteatt(outfile,'SG_WEST', 'coordinates','TIME PRES');
ncwriteatt(outfile,'SG_WEST', 'long_name','sea water salinity at western boundary 57.5N/-12.3W');
ncwriteatt(outfile,'SG_WEST','reference_scale','TEOS-10');
ncwriteatt(outfile,'SG_WEST', 'QC_indicator','good data');
ncwriteatt(outfile,'SG_WEST', 'valid_min',single(0));  % CHANGE
ncwriteatt(outfile,'SG_WEST', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'U_WEST_1', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'U_WEST_1', single(U_WEST_1));
ncwriteatt(outfile,'U_WEST_1', 'standard_name','velocity');
ncwriteatt(outfile,'U_WEST_1', 'units','cm/s');
ncwriteatt(outfile,'U_WEST_1', '_FillValue', single(99999));
ncwriteatt(outfile,'U_WEST_1', 'coordinates','TIME PRES');
ncwriteatt(outfile,'U_WEST_1', 'long_name','current speed u-direction at western boundary 1 XXN/XXW');
ncwriteatt(outfile,'U_WEST_1', 'QC_indicator','good data');
% ncwriteatt(outfile,'U_WEST_1', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'U_WEST_1', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'V_WEST_1', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'V_WEST_1', single(V_WEST_1));
ncwriteatt(outfile,'V_WEST_1', 'standard_name','velocity');
ncwriteatt(outfile,'V_WEST_1', 'units','cm/s');
ncwriteatt(outfile,'V_WEST_1', '_FillValue', single(99999));
ncwriteatt(outfile,'V_WEST_1', 'coordinates','TIME PRES');
ncwriteatt(outfile,'V_WEST_1', 'long_name','current speed v-direction at western boundary 1 XXN/XXW');
ncwriteatt(outfile,'V_WEST_1', 'QC_indicator','good data');
% ncwriteatt(outfile,'U_WEST_1', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'U_WEST_1', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'W_WEST_1', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'W_WEST_1', single(W_WEST_1));
ncwriteatt(outfile,'W_WEST_1', 'standard_name','velocity');
ncwriteatt(outfile,'W_WEST_1', 'units','cm/s');
ncwriteatt(outfile,'W_WEST_1', '_FillValue', single(99999));
ncwriteatt(outfile,'W_WEST_1', 'coordinates','TIME PRES');
ncwriteatt(outfile,'W_WEST_1', 'long_name','current speed w-direction at western boundary 1 XXN/XXW');
ncwriteatt(outfile,'W_WEST_1', 'QC_indicator','good data');
% ncwriteatt(outfile,'U_WEST_1', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'U_WEST_1', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'U_WEST_2', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'U_WEST_2', single(U_WEST_2));
ncwriteatt(outfile,'U_WEST_2', 'standard_name','velocity');
ncwriteatt(outfile,'U_WEST_2', 'units','cm/s');
ncwriteatt(outfile,'U_WEST_2', '_FillValue', single(99999));
ncwriteatt(outfile,'U_WEST_2', 'coordinates','TIME PRES');
ncwriteatt(outfile,'U_WEST_2', 'long_name','current speed u-direction at western boundary 2 XXN/XXW');
ncwriteatt(outfile,'U_WEST_2', 'QC_indicator','good data');
% ncwriteatt(outfile,'U_WEST_2', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'U_WEST_2', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'V_WEST_2', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'V_WEST_2', single(V_WEST_2));
ncwriteatt(outfile,'V_WEST_2', 'standard_name','velocity');
ncwriteatt(outfile,'V_WEST_2', 'units','cm/s');
ncwriteatt(outfile,'V_WEST_2', '_FillValue', single(99999));
ncwriteatt(outfile,'V_WEST_2', 'coordinates','TIME PRES');
ncwriteatt(outfile,'V_WEST_2', 'long_name','current speed v-direction at western boundary 2 XXN/XXW');
ncwriteatt(outfile,'V_WEST_2', 'QC_indicator','good data');
% ncwriteatt(outfile,'U_WEST_2', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'U_WEST_2', 'valid_max',single(40));  % CHANGE


nccreate(outfile,'W_WEST_2', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'W_WEST_2', single(W_WEST_2));
ncwriteatt(outfile,'W_WEST_2', 'standard_name','velocity');
ncwriteatt(outfile,'W_WEST_2', 'units','cm/s');
ncwriteatt(outfile,'W_WEST_2', '_FillValue', single(99999));
ncwriteatt(outfile,'W_WEST_2', 'coordinates','TIME PRES');
ncwriteatt(outfile,'W_WEST_2', 'long_name','current speed w-direction at western boundary 2 XXN/XXW');
ncwriteatt(outfile,'W_WEST_2', 'QC_indicator','good data');
% ncwriteatt(outfile,'U_WEST_2', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'U_WEST_2', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'U_EAST', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'U_EAST', single(U_EAST));
ncwriteatt(outfile,'U_EAST', 'standard_name','velocity');
ncwriteatt(outfile,'U_EAST', 'units','cm/s');
ncwriteatt(outfile,'U_EAST', '_FillValue', single(99999));
ncwriteatt(outfile,'U_EAST', 'coordinates','TIME PRES');
ncwriteatt(outfile,'U_EAST', 'long_name','current speed u-direction at eastern boundary 57.1N/-9.6W');
ncwriteatt(outfile,'U_EAST', 'QC_indicator','good data');
% ncwriteatt(outfile,'U_EAST', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'U_EAST', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'V_EAST', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'V_EAST', single(V_EAST));
ncwriteatt(outfile,'V_EAST', 'standard_name','velocity');
ncwriteatt(outfile,'V_EAST', 'units','cm/s');
ncwriteatt(outfile,'V_EAST', '_FillValue', single(99999));
ncwriteatt(outfile,'V_EAST', 'coordinates','TIME PRES');
ncwriteatt(outfile,'V_EAST', 'long_name','current speed v-direction at eastern boundary 57.1N/-9.6W');
ncwriteatt(outfile,'V_EAST', 'QC_indicator','good data');
% ncwriteatt(outfile,'U_EAST', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'U_EAST', 'valid_max',single(40));  % CHANGE

nccreate(outfile,'W_EAST', 'Dimensions',{'PRES',PressureDim, 'TIME',TimeDim}, 'Datatype','single');
ncwrite(outfile,'W_EAST', single(W_EAST));
ncwriteatt(outfile,'W_EAST', 'standard_name','velocity');
ncwriteatt(outfile,'W_EAST', 'units','cm/s');
ncwriteatt(outfile,'W_EAST', '_FillValue', single(99999));
ncwriteatt(outfile,'W_EAST', 'coordinates','TIME PRES');
ncwriteatt(outfile,'W_EAST', 'long_name','current speed w-direction at eastern boundary 57.1N/-9.6W');
ncwriteatt(outfile,'W_EAST', 'QC_indicator','good data');
% ncwriteatt(outfile,'U_EAST', 'valid_min',single(0));  % CHANGE
% ncwriteatt(outfile,'U_EAST', 'valid_max',single(40));  % CHANGE

fprintf('Finish %s\n', datestr(now));

clearvars -except outfile
finfo2=ncinfo(outfile);

