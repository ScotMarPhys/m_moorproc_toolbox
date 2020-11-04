%%  Experimental Development Code for the Merging of the Eastern Boundary
%%  MicroCAT Data
%
% function [MERG_REVISION,MERG_AUTHOR,MERG_DATE] =
%              eb_merging_v2_2012(TS_CLIMATOLOGY,TS_CLIMATOLOGY_TP,S_CLIMATOLOGY_NAME,EB_funct,EB_FILE,jg_end)
%  Paul Wright  March 2010
%  Updated Ben Moat March and Nov 2012
%
%  eb_merging_v1.m  achieves the following:  The selcted P, S and T data are obtained
%  from the mooring .microcat files and interpolated onto the full timeseries
%  time grid, sorted into pressure levels and then gridded.
%
%  uses meta data files to select the MicroCATs.
%
%  it would be possible to simplify this further but it is debatable
%  whether it may be better to leave a bit "open" for the ease of updating
%  on an annual basis.  Clearer too.
%
%  Notes specific to each deployment period are included within the loops.
%
% INPUT
%  TS_CLIMATOLOGY    -- climatology, eg. 'slope' or 'deep'
%  TS_CLIMATOLOGY_TP -- climatology, monthly or annual
%  TS_CLIMATOLOGY_NAME -- name of climatology, argo or hbase
%  EB_funct          -- name of this function. Store in structure for tracability
%  EB_FILE           -- name of the OUTPUT files from this function
%  jg_end            -- julian number of the end of the calulation
%
% OUTPUT
%  MERG_REVISION     -- subversion revision number for this script
%  MERG_AUTHOR       -- subversion author of this script
%  MERG_DATE	     -- date this script was banked with the subversion server
%
%  ADDITIONAL FILES LOADED INTO FUNCTION
%  e.g. :  eb_rapid_1.dat % these files contain the Mooring: SerialNumber: MicrocatNumber: Depth: Lon
%                         % used to make the eastern boundary profile.
%                         % they are stored in ./merging/mat_files/
%  e.g. :  ebh5_1_200402_grid.mat % gridded file from BODC for each mooring
%                         % they are stored in ./merging/gridded_profiles/
%                           C             5x4801
%                           P             5x4801
%                           Pfs           5x793  ***  pressure time series from 5 microcats
%                           S             5x4801
%                           SGfs         22x793
%                           Sfs           5x793  ***  salinity time series from 5 microcats
%                           T             5x4801
%                           TGfs         22x793
%                           Tfs           5x793  *** temperature time series from 5 microcats
%                           co            1x1
%                           jd            1x793  *** matlab time since deployment
%                           jd_grid       1x4801
%                           p_grid       22x1
% RAPID 1:
% /noc/users/pstar/rpdmoc/rapid/data/amoc/pauls_projects/merging/mat_files/eb_rapid_1.dat
% /noc/users/pstar/rpdmoc/rapid/data/amoc/pauls_projects/merging/gridded_profiles/ebh5_1_200402_grid.mat
% RAPID 2:
% /noc/users/pstar/rpdmoc/rapid/data/amoc/pauls_projects/merging/mat_files/eb_rapid_2.dat
% /noc/users/pstar/rpdmoc/rapid/data/amoc/grdat/eb_2005_grid_merge.mat'
% RAPID 3a:
% /noc/users/pstar/rpdmoc/rapid/data/amoc/pauls_projects/merging/mat_files/eb_rapid_3a.dat
%
%  Code History
% -------------
% 16 April 2010 - Code eb_merging_v1.m written and generally
% 	 functioning well.  Some issues at the joins between the deployment
% 	 periods and the occasional spkie to get rid of.  Issue with the gridding
% 	 not stopping the same way as hydro_grid_merge.m
% 16 April 2010 - changed code to read the gridded files rather than the
%	 microcat files, as they have been individually worked on and repaired,
% 	 corrected and de-spiked where necessary... now eb_merging_v2.m
% 19 April 2010 - removed final data point of each MicroCAT to attempt to
% 	remove the spiking at the joins between the deployments.
% 26 May 2010 - changed gridded function to con_tprof_v3.m which stops the
% 	gridding at the upper MicroCAT.  Interpolates over data gaps less than a
% 	year in order to match Kanzow's original code and present data for the
% 	MOC code.
% 9 June 2010 - PGW changed the climatology TSclim
% 28 June 2010 - changed climatology back to original.  Updated to include
% 	D344 data
% 11 August 2010 - Tidied up and improved comments slightly
% June 2011 - Updated merging Apr 2009 - Dec 2010
% Sept 2013 - added in name of cliamtology (e.g. hbase or argo)

function [MERG_REVISION,MERG_AUTHOR,MERG_DATE] = ...
    eb_merging_v2_2012(TS_CLIMATOLOGY,TS_CLIMATOLOGY_TP,TS_CLIMATOLOGY_NAME,EB_funct,EB_FILE,jg_end)

mfilename; function_name = mfilename('fullpath');

%Subversion code history
% do not alter the next three lines
REV             ='$Rev: 237 $';
REV_AUTHOR      ='$Author: bim $';
REV_DATE        ='$Date: 2016-03-10 14:02:10 +0000 (Thu, 10 Mar 2016) $';
MERG_REVISION    = REV(6:8);MERG_AUTHOR = REV_AUTHOR(9:end-1);MERG_DATE = REV_DATE(7:end-1);
%%%%%%%

%        clear all; close all; clc;
warning off
gridding        = true  ;  % turns off the hydro_gridding to save time!
bathy           = false ;  % turns on/off the bathy charts. off = flase
mc_check_plot   = false ;  % turns on/off the microcat check plots. off =flase

% convert the filenames to sttings
OUT_FILE1 = EB_FILE(1);OUT_FILE2 = EB_FILE(2);
OUT_FILE(1)=cellstr(OUT_FILE1);OUT_FILE(2)=cellstr(OUT_FILE2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  1.  DIRECTORIES AND PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get IP address to check if the machine is on the NOC
% network or on a ship
address = java.net.InetAddress.getLocalHost;
IPaddress = char(address.getHostAddress);

% set paths
if ismac % set paths for using a mac
    basedir = '/Volumes/rpdmoc/';
elseif isunix
    if IPaddress(1:7)==num2str(139.166) % AT NOC
        basedir = '/noc/mpoc/rpdmoc/';
    else
        basedir = '/noc/users/pstar/rpdmoc/'; % ON a SHIP
    end
end

% add in paths to other functions
addpath(genpath([basedir,'rapid/data/amoc/m/']));

% set paths
boundarydir = [basedir, 'rapid/data/amoc/boundary_files/']
grdatdir    = [basedir, 'rapid/data/amoc/grdat/'];
procdir     = [basedir, 'rapid/data/moor/proc/'];
hydrodir    = [basedir, 'rapid/data/amoc/gridded_profiles/'];
bathyfile   = [basedir, 'rapid/data/amoc/gridded_profiles/Discovery_bathy_data.mat'];

% name of the annual climatologies
% code to read in OLD annual files created by Torsen
% removed (but kept if required)  as not used
if strcmp(TS_CLIMATOLOGY_TP,'annual_vergrad')
    TSclim    = [basedir, 'rapid/data/moor/proc/hydro_grid/clim/' ...
        'hydro_vergrad_eb_' TS_CLIMATOLOGY '.mat'];
    EB_merg.TSclim=TSclim;
    display('RUNNING: OLD annual vergrad climatology')
else

          TSclim    = [basedir, 'rapid/data/moor/proc/hydro_grid/clim/' ...
                       'eb_' TS_CLIMATOLOGY '_' TS_CLIMATOLOGY_NAME '_' TS_CLIMATOLOGY_TP '_climatology.mat'];

% e.g.
%/noc/mpoc/rpdmoc/rapid/data/moor/proc/hydro_grid/clim/eb_slope_argo_seasonal_climatology.mat
%  dsdp_slope_argo               555x12             53280  double              
%  dtdp_slope_argo               555x12             53280  double              
%  hydrobase_slope_argo_lon        1x4                 32  double              
%  pressure_slope_argo           100x1                800  double              
%  salinity_slope_argo           100x12              9600  double              
%  temperature_slope_argo        100x12              9600  double              
%  tgrid_argo                      1x555             4440  double
%
% /noc/mpoc/rpdmoc/rapid/data/moor/proc/hydro_grid/clim/eb_slope_hbase_seasonal_climatology.mat
%  d2sdp2_slope_on_p_hbase      311x12             29856  double              
%  d2tdp2_slope_on_p_hbase      311x12             29856  double              
%  dsdp_slope_on_t_hbase        555x12             53280  double              
%  dtdp_slope_on_t_hbase        555x12             53280  double              
%  longitude_slope_hbase          4x1                 32  double              
%  pg_hbase                     311x1               2488  double              
%  pressure_slope_hbase         311x12             29856  double              
%  salinity_slope_hbase         311x12             29856  double              
%  temperature_slope_hbase      311x12             29856  double              
%  tgrid_hbase


end

display(['T/S climatology: ' TSclim]);
display(' ')

% T/S climatology variables
%       dsdp       1x237     dS/dP
%  	dtdp       1x237     dT/dP
%  	p          1x237     pressure
%  	s          1x237     salinity
%       tgrid      1x237     temperature grid ?


% temporal interpolation below idepth
        idepth = 220; % multiples of 20; 
                      
int_step        = 10;
dum             = -9999.000; % absent data value
jg_start        = julian(2004,3,1);
JG              = jg_start: 0.5: jg_end; % full time series using 2 samples per day
pg              = 0:20:4820; % depths in 20dbar bins
ratio_factor    = sw_c3515;
t90_68          = 1.00024 ;   % converts t68 to t90
col             = {'r','b','m','c','g','y','r--','b--','m--','c--','g--','y--','r:',...
    'b:','m:','c:','g:','y:','r.','b.','m.','c.','g.','y.','r'};
month 	        = ['JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'];

GTV             = gregorian(JG); % Gregorian time vector
% use to define the month for the seasonal climatology
start_date      = gregorian(jg_start);
start_month     = month(start_date(2),:);
end_date        = gregorian(jg_end);
end_month       = month(end_date(2),:);
disp('----------------------------');
disp('EASTERN BOUNDARY MERGING SEGMENT');
disp('  ');
disp(['START DATE  = ', num2str(start_date(3)),' ', start_month, ' ', num2str(start_date(1))]);
disp(['END DATE    = ', num2str(end_date(3)),' ', end_month, ' ', num2str(end_date(1))]);
disp('----------------------------');
disp('  ');

% ADDED, to use matlab time
mat_start = datenum(start_date(6),start_date(5),start_date(4),start_date(3),start_date(2),start_date(1));
mat_end = datenum(end_date(6),end_date(5),end_date(4),end_date(3),end_date(2),end_date(1));
mattime = mat_start: 0.5: mat_end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2a.  LOOP FOR RAPID 1 (D277 --> CD170)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Notes:
%  ------
%  The first deployment period has a number of issues.
%  1.  No shallow MicroCATs were recovered inshore so Tortsen Kanzow
%  extrapolated the upper measurement anomalies of EBH5 upwards to create two
%  "virtual" MicroCATs.  This was considered more realistic than just
%  allowing the ds/dp and dT/dp integration to extend on to the surface
%  2.  The temperature data from the BPR on EB1 has been included, although
%  it is considered to be less reliable than the MicroCAT temperature
%  sensor.  The salinity is estimated from the S-T relationship.


disp('---------  RAPID 1 (D277 --> CD170) ---------')
% metadat file holds the moorings to take data from
% to build the eastern boundary profile
% deepest moorings chosen the create a profile of all depths.
[mooring1, sn1, mc1, z1, lon1] = textread([boundarydir, 'eb_rapid_1.dat'],...
    '%s %d %d %d %6.2f','headerlines ',5);
% mooring1 - mooring name (e.g. ebh5_1_200402,
%                               ebh5 = array name
%                               1    = RAPID depolyment number
%                               2004 = year of deployment
%                               02   = UKCORS mooring number
% sn1      - microcat serial number
% mc1      - microcat number (lowest is the deepest)
% z1       - depth of microcat
% lon1     - position (longitude of the array)

T1 = zeros(length(sn1), length(JG)); % initialise temperature array
P1 = zeros(length(sn1), length(JG)); % initialise pressure array
S1 = zeros(length(sn1), length(JG)); % initialise salinity array

j = 1; length(sn1);
for ii = 1: length(sn1) % loop through each microcat at a specific depth
    
    % load in the gridded file from BODC
    % for each mooring
    infile = [hydrodir, mooring1{ii,:}, '_grid.mat'];
    load(infile);
    % jd(1) =2453066; gregorian(2453066) = [2004    3    1         0         0    0]
    % looks like the microcat data is every 2 hours
    % gridded data that comes back is erby 6 hours
    %  C             5x4801
    %  P             5x4801 *** Presure (2 hourly)
    %  Pfs           5x793  *** pressure time series from 5 microcats (6 hourly)
    %  S             5x4801 *** salinity
    %  SGfs         22x793  *** gridded salinity
    %  Sfs           5x793  ***  salinity time series from 5 microcats
    %  T             5x4801
    %  TGfs         22x793  *** gridded Temperature
    %  Tfs           5x793  *** temperature time series from 5 microcats
    %  co            1x1
    %  jd            1x793  *** matlab time since deployment
    %  jd_grid       1x4801
    %  p_grid       22x1  presure e.g 540 to 960 dbar
    
    sampling_rate = round(1./median(diff(jd))); %nominal sampling rate [per day]
    % interpolating the time series from each microcat onto the full time series array (JG)
    % do this for each intrument
    salinity    = interp1(jd(1:end-1), Sfs(mc1(ii),1:end-1)', JG)';
    temp        = interp1(jd(1:end-1), Tfs(mc1(ii),1:end-1)', JG)';
    pressure    = interp1(jd(1:end-1), Pfs(mc1(ii),1:end-1)', JG)';
    
    % store the time series as a row in an array whos depth is the rows.
    % i has the same index value as j, so use i instead
    %    Sfs1(j,:) = salinity;
    %    Pfs1(j,:) = pressure;
    %    Tfs1(j,:) = temp;
    Sfs1(ii,:) = salinity;
    Pfs1(ii,:) = pressure;
    Tfs1(ii,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(20041)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs1(ii, :), col{:,ii},'LineWidth',2);
        title('RAPID 1 - SALINITY')
        
        figure(20042)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs1(ii, :), col{:,ii},'LineWidth',2);
        title('RAPID 1 - TEMPERATURE')
        
        figure(20043)
        hold on; box on;
        plot(JG - julian(2004,1,1), Pfs1(ii, :), col{:,ii},'LineWidth',2);
        title('RAPID 1 - PRESSURE')
    end
end

disp('   ')
disp(' ')
disp('   ')
disp('*************************************');
disp(' For RAPID 1 (2004) KANZOW included the Seaguage temperature and pressure ')
disp(' data at the bottom of the EB1 mooring');
disp('   ');
disp(' loading from the gridded file : eb1_1_200409_grid.mat')
disp('*************************************');
load([procdir, 'hydro_grid/eb1_1_200409_grid.mat']);

% interpolating the time series from each microcat onto the full time series array (JG)

Sfs = interp1(jd, Sfs', JG)';
Pfs = interp1(jd, Pfs', JG)';
Tfs = interp1(jd, Tfs', JG)';


%add in a bottom of the EB1 mooring to the eastern boundary profile
% j is basically ii+1
Sfs1(j,:) = Sfs(7,:); % 7 is the deepest instrument (seaguage)
Pfs1(j,:) = Pfs(7,:); % 7 is the deepest instrument (seaguage)
Tfs1(j,:) = Tfs(7,:); % 7 is the deepest instrument (seaguage)

disp('  ')

disp('*************************************');
disp(' For RAPID 1 (2004) KANZOW extrapolated EBH5  ')
disp('  ');
disp(' Ref:  rpdmoc/users/tok/rapid/work/ebh5_extrapolation.m   ');
disp(' loading from the gridded file : ebh5_1_200402_grid_extrapolated.mat')
disp('*************************************');
% shallowest microcat is at a depth of 565 (see z1(1))
% need to extrapolate to close to the surface using EBH5

load([procdir, 'hydro_grid/ebh5_1_200402_grid_extrapolated.mat']);

% interpolating the time series from each microcat onto the full time series array (JG)
% only use the two shallowest instruments
Sfsx = interp1(jd(1:end-1), Sfs(1:2,[1:end-1])', JG)';
Pfsx = interp1(jd(1:end-1), Pfs(1:2,[1:end-1])', JG)';
Tfsx = interp1(jd(1:end-1), Tfs(1:2,[1:end-1])', JG)';

if mc_check_plot
    figure(20041)
    plot(JG - julian(2004,1,1), Sfsx, 'y--','LineWidth',2);
    figure(20042)
    plot(JG - julian(2004,1,1), Tfsx, 'y--','LineWidth',2);
    figure(20043)
    plot(JG - julian(2004,1,1), Pfsx, 'y--','LineWidth',2);
end

% why is the surface data put in the deep ??
% Sort on pressure later, so doesn't matter.
Tfs1 = [Tfs1;Tfsx];
Sfs1 = [Sfs1; Sfsx];
Pfs1 = [Pfs1; Pfsx];

%-------------------

load([basedir, 'rapid/data/amoc/grdat/eb_2004_extra_grid_merge.mat']);
jd1 = jd; SGfs1 = SGfs; TGfs1 = TGfs; PG1 = p_grid;
if mc_check_plot
    %    jd1 = jd; SGfs1 = SGfs; TGfs1 = TGfs; PG1 = p_grid;
    figure(20041)
    plot(jd - julian(2004,1,1), S, 'Color', [0.5, 0.5, 0.5])
    figure(20042)
    plot(jd - julian(2004,1,1), T, 'Color', [0.5, 0.5, 0.5])
    figure(20043)
    plot(jd - julian(2004,1,1), P, 'Color', [0.5, 0.5, 0.5])
    axis ij
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2b.  LOOP FOR RAPID 2 (CD170 --> CD177)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Notes:
%  --------
%  The second deployment period has been difficult to check as the MicroCAT
%  data was not included in Torsten's merged files.  The batteries failed
%  before planned recovery time as the sensors had been loaded with
%  alkaline cells rather than Lithium ones.  The established cut off point
%  to match the previous work by Torsten is 02 Feb 2006.


disp('---------  RAPID 2 (CD170 --> CD177) ---------')
disp('------ short record as batteries failed -----')


load([basedir, 'rapid/data/amoc/grdat/eb_2005_grid_merge.mat']);
jd2 = jd; SGfs2 = SGfs; TGfs2 = TGfs; PG2 = p_grid;
if mc_check_plot
    figure(20051)
    plot(jd - julian(2004,1,1), SGfs, 'Color', [0.5, 0.5, 0.5])
    figure(20052)
    plot(jd - julian(2004,1,1), TGfs, 'Color', [0.5, 0.5, 0.5])
    
end
[mooring2, sn2, mc2, z2, lon2] = textread([boundarydir, 'eb_rapid_2.dat'],'%s %d %d %d %6.2f','headerlines ',5);

T2 = zeros(length(sn2), length(JG));
P2 = zeros(length(sn2), length(JG));
S2 = zeros(length(sn2), length(JG));

i = 1; j = 1;
for i = 1: length(sn2)
    
    infile = [hydrodir, mooring2{i,:}, '_grid.mat'];
    load(infile);
    sampling_rate = round(1./median(diff(jd_grid)));  % nominal sampling rate [per day]
    salinity    = interp1(jd(1:end-1), Sfs(mc2(i),1:end-1)', JG)';
    temp        = interp1(jd(1:end-1), Tfs(mc2(i),1:end-1)', JG)';
    pressure    = interp1(jd(1:end-1), Pfs(mc2(i),1:end-1)', JG)';
    
    Sfs2(j,:) = salinity;
    Pfs2(j,:) = pressure;
    Tfs2(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(20051)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs2(i, :), col{:,i},'LineWidth',2);
        title('RAPID 2 - SALINITY')
        
        figure(20052)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs2(i, :), col{:,i},'LineWidth',2);
        title('RAPID 2 - TEMPERATURE')
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2c.  LOOP FOR RAPID 3a (D304 --> F343/4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Notes:
%  ------
%  RAPID 3 was complicated by having two cruises with overlapping
%  moorings...  Some batteries were discovered to have problems (RAPID2) but
%  this was not discovered until after a number of moorings had been
%  deployed...
%  There are two time periods - some moorings (EB1, EB2, EBH4 and EBH5)
%  were rebatteried so are a full record and extend to the end of RAPID3a
%  abd RAPID3b.  The middle range EBH2 --> EBHi were turned around on the
%  Poseidon cruise.   The inshore moorings started much later than the
%  deeper moorings.  Originally Torsten extrapolated them back to the start
%  point in a straight line - but I have not found the code he used for
%  this so I have left them as is and rely on the gridding dT/dp
%  integration, which looks more realistic anyway.


disp('---------  RAPID 3a (D304 --> F343/4)')
[mooring3a, sn3a, mc3a, z3a, lon3a] = textread([boundarydir, 'eb_rapid_3a.dat'],'%s %d %d %d %6.2f','headerlines ',5);

% mooring1 - mooring name (e.g. ebh5_1_200402,
%                               ebh5 = array name
%                               1    = RAPID depolyment number
%                               2004 = year of deployment
%                               02   = UKCORS mooring number
% sn1      - microcat serial number
% mc1      - microcat number (lowest is the deepest)
% z1       - depth of microcat
% lon1     - position (longitude of the array)
% eg. eb_rapid_3a.dat
%EASTERN BOUNDARY - RAPID 3a
%Start: Oct 2006
%End:  May 2006 and/or Oct 2007 (batteries failed)
%Columns =  Mooring: SerialNumber: MicrocatNumber: Depth: Lon
%----------------------------------------------
%ebh5_3_200610    3207    001    45    -13.36
%ebh5_3_200610    3208    002    90    -13.36
%ebh5_3_200610    3209    003    158   -13.36
%ebh5_3_200610    4473    004    234   -13.36
%ebh4_3_200611    3212    001    312   -13.53
% etc .........


T3a = zeros(length(sn3a), length(JG));
P3a = zeros(length(sn3a), length(JG));
S3a = zeros(length(sn3a), length(JG));


i = 1; j = 1;
for i = 1: length(sn3a)
    %    display( [hydrodir, mooring3a{i,:}, '_grid.mat'] )
    infile = [hydrodir, mooring3a{i,:}, '_grid.mat'];
    load(infile);
    sampling_rate = round(1./median(diff(jd_grid)));%nominal sampling rate [per day]
    salinity    = interp1(jd(1:end-1), Sfs(mc3a(i),1:end-1)', JG)';
    temp        = interp1(jd(1:end-1), Tfs(mc3a(i),1:end-1)', JG)';
    pressure    = interp1(jd(1:end-1), Pfs(mc3a(i),1:end-1)', JG)';
    
    Sfs3a(j,:)  = salinity;
    Pfs3a(j,:)  = pressure;
    Tfs3a(j,:)  = temp;
    M3a{j,:}    = [mooring3a{i,:}, '_',num2str(sn3a(i))];
    
    j = j + 1;
    
    if mc_check_plot
        figure(20061)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs3a(i, :), col{:,i},'LineWidth',2);
        title('RAPID 3 - SALINITY')
        
        figure(20062)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs3a(i, :), col{:,i},'LineWidth',2);
        title('RAPID 3 - TEMPERATURE')
        
        figure(20063)
        hold on; box on;
        plot(JG - julian(2004,1,1), Pfs3a(i, :), col{:,i},'LineWidth',2);
        title('RAPID 3 - PRESSURE')
    end
    
    
end

% what is this for?
% what file does it look for?
load([basedir, 'rapid/data/amoc/grdat/ebh_2006_grid_merge.mat']);
jd3a = jd; SGfs3a = SGfs; TGfs3a = TGfs; PG3a = p_grid;
if mc_check_plot
    figure(20061)
    plot(jd - julian(2004,1,1), S, 'k')
    figure(20062)
    plot(jd - julian(2004,1,1), T, 'k')
    figure(20063)
    plot(jd - julian(2004,1,1), P, 'k')
    axis ij
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2d.  LOOP FOR RAPID 3b (F343/4 --> D324)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Notes:
%  -------
%  RAPID3b again is complicated by the battery failure.  Only the
%  instruments that were started on the Poseidon cruises are included here.
%  The others continue on from RAPID 3a.

disp('---------  RAPID 3b (F343/4 --> D324)')
[mooring3b, sn3b, mc3b, z3b, lon3b] = textread([boundarydir, 'eb_rapid_3b.dat'],'%s %d %d %d %6.2f','headerlines ',5);  % note!!!

T3b = zeros(length(sn3b), length(JG));
P3b = zeros(length(sn3b), length(JG));
S3b = zeros(length(sn3b), length(JG));


i = 1; j = 1;
for i = 1: length(sn3b)
    
    infile = [hydrodir, mooring3b{i,:}, '_grid.mat'];
    load(infile);
    sampling_rate = round(1./median(diff(jd)));%nominal sampling rate [per day]
    salinity    = interp1(jd(1:end-1), Sfs(mc3b(i),1:end-1)', JG)';
    temp        = interp1(jd(1:end-1), Tfs(mc3b(i),1:end-1)', JG)';
    pressure    = interp1(jd(1:end-1), Pfs(mc3b(i),1:end-1)', JG)';
    
    Sfs3b(j,:) = salinity;
    Pfs3b(j,:) = pressure;
    Tfs3b(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(20061)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs3b(i, :), col{:,i},'LineWidth',2);
        title('RAPID 3 - SALINITY')
        
        figure(20062)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs3b(i, :), col{:,i},'LineWidth',2);
        title('RAPID 3 - TEMPERATURE')
        
        figure(20063)
        hold on; box on;
        plot(JG - julian(2004,1,1), Pfs3b(i, :), col{:,i},'LineWidth',2);
        title('RAPID 3 - PRESSURE')
    end
    
end


disp('************************************')
disp(' For RAPID 3b (2006) KANZOW extrapolated EBM  ')
disp('  ');
disp(' Ref: moor/proc/hydro_grid/ebm_2007_extrapolation.m   ');
disp(' loading from the gridded file : ebm_2007_grid_extrapolated.mat')
disp('*************************************');

load([procdir, 'hydro_grid/ebm_2007_grid_extrapolated.mat']);

Sfsx = interp1(jd(1:end-1), Sfs(1,1:end-1)', JG)';
Pfsx = interp1(jd(1:end-1), Pfs(1,1:end-1)', JG)';
Tfsx = interp1(jd(1:end-1), Tfs(1,1:end-1)', JG)';
if mc_check_plot
    figure(20061)
    plot(JG - julian(2004,1,1), Sfsx, 'y--','LineWidth',2);
    figure(20062)
    plot(JG - julian(2004,1,1), Tfsx, 'y--','LineWidth',2);
    figure(20063)
    plot(JG - julian(2004,1,1), Pfsx, 'y--','LineWidth',2);
end

Tfs3b = [Tfs3b; Tfsx'];
Sfs3b = [Sfs3b; Sfsx'];
Pfs3b = [Pfs3b; Pfsx'];


load([basedir, 'rapid/data/amoc/grdat/ebh_2007_grid_merge.mat']);
jd3b = jd; SGfs3b = SGfs; TGfs3b = TGfs; PG3b = p_grid;
if mc_check_plot
    figure(20061)
    plot(jd - julian(2004,1,1), S, 'Color', [0.5, 0.5, 0.5])
    figure(20062)
    plot(jd - julian(2004,1,1), T, 'Color', [0.5, 0.5, 0.5])
    figure(20063)
    plot(jd - julian(2004,1,1), P, 'Color', [0.5, 0.5, 0.5])
    axis ij
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2e.  LOOP FOR RAPID 4 (D324 --> D334)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Notes:
%  -------
%  This appears to be the first straightforward year!


disp('---------  RAPID 4 (D324 --> D334) ---------')
[mooring4, sn4, mc4, z4, lon4] = textread([boundarydir, 'eb_rapid_4.dat'],'%s %d %d %d %6.2f','headerlines ',5);

T4 = zeros(length(sn4), length(JG));
P4 = zeros(length(sn4), length(JG));
S4 = zeros(length(sn4), length(JG));

% spike in dynamic heights
% remove the last days of RAPID 4 and the start of RAPID 5
% was set at end-1, as 12 hour sampling increased this to -3
% to remove 1 days worth of data

i = 1; j = 1;
for i = 1: length(sn4)
    
    infile = [hydrodir, mooring4{i,:}, '_grid.mat'];
    load(infile);
    sampling_rate = round(1./median(diff(jd_grid)));%nominal sampling rate [per day]
    salinity    = interp1(jd(1:end-1), Sfs(mc4(i),1:end-1)', JG)'; %
    temp        = interp1(jd(1:end-1), Tfs(mc4(i),1:end-1)', JG)'; % 
    pressure    = interp1(jd(1:end-1), Pfs(mc4(i),1:end-1)', JG)'; % 
    
    Sfs4(j,:) = salinity;
    Pfs4(j,:) = pressure;
    Tfs4(j,:) = temp;
    
    j = j + 1;
    if mc_check_plot
        figure(20071)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs4(i, :), col{:,i},'LineWidth',2);
        title('RAPID 4 - SALINITY')
        
        figure(20072)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs4(i, :), col{:,i},'LineWidth',2);
        title('RAPID 4 - TEMPERATURE')
        
        figure(20073)
        hold on; box on;
        plot(JG - julian(2004,1,1), Pfs4(i, :), col{:,i},'LineWidth',2);
        title('RAPID 4 - PRESSURE')
        
    end
    
end

% two surface measurments from 1 mooring are causing spikes in the dynamic height calculation
% Remove the measurements
%Tfs4(1,3444)=nan;Tfs4(1,3443)=nan;
%Sfs4(1,3444)=nan;Sfs4(1,3443)=nan;
%Pfs4(1,3444)=nan;Pfs4(1,3443)=nan;
% DOESN'T IMPROVE THE DYNAMIC HIEGHT ANOMOLY
% NOT IMPLEMENTED

load([basedir, 'rapid/data/amoc/grdat/eb500_2007_grid_merge.mat']);
jd4 = jd; SGfs4 = SGfs; TGfs4 = TGfs; PG4 = p_grid;
if mc_check_plot
    figure(20071)
    plot(jd - julian(2004,1,1), S, 'k')
    figure(20072)
    plot(jd - julian(2004,1,1), T, 'k')
    figure(20073)
    plot(jd - julian(2004,1,1), P, 'k')
    axis ij
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2f. LOOP FOR RAPID 5 (D334 --> D344)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Notes:
%  ------
%  This appears to be the second straightforward year!


disp('---------  RAPID 5 (D334 --> D344) ---------')
[mooring5, sn5, mc5, z5, lon5] = textread([boundarydir, 'eb_rapid_5.dat'],'%s %d %d %d %6.2f','headerlines ',5);

T5 = zeros(length(sn5), length(JG));
P5 = zeros(length(sn5), length(JG));
S5 = zeros(length(sn5), length(JG));


i = 1; j = 1;
for i = 1: length(sn5)
    
    infile = [hydrodir, mooring5{i,:}, '_grid.mat'];
    load(infile);
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jd(1:end-1), Sfs(mc5(i),1:end-1)', JG)';% 
    temp        = interp1(jd(1:end-1), Tfs(mc5(i),1:end-1)', JG)';% 
    pressure    = interp1(jd(1:end-1), Pfs(mc5(i),1:end-1)', JG)';% 
    Sfs5(j,:) = salinity;
    Pfs5(j,:) = pressure;
    Tfs5(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(20081)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs5(i, :), col{:,i},'LineWidth',2);
        title('RAPID 5 - SALINITY')
        
        figure(20082)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs5(i, :), col{:,i},'LineWidth',2);
        title('RAPID 5 - TEMPERATURE')
        
        figure(20083)
        hold on; box on;
        plot(JG - julian(2004,1,1), Pfs5(i, :), col{:,i},'LineWidth',2);
        title('RAPID 5 - PRESSURE')
    end
    
end

% deep measurments from 1 mooring maybe causing spikes in the dynamic height calcualtion
% Remove the measurements from the start of the mooring
% deep measurements exist
% take it back a further * days so we go back to measurements from all moorings. no partial (no nans)
%Tfs5(:,3413:3456)=nan;
%Sfs5(:,3413:3456)=nan;
%Pfs5(:,3413:3456)=nan;
% DOESN'T IMPROVE THE DYNAMIC HIEGHT ANOMOLY
% NOT IMPLEMENTED

load([basedir, 'rapid/data/amoc/grdat/eb_2008_grid_merge.mat']);
jd5 = jd; SGfs5 = SGfs; TGfs5 = TGfs; PG5 = p_grid;
if mc_check_plot
    figure(20081)
    plot(jd - julian(2004,1,1), S, 'k')
    figure(20082)
    plot(jd - julian(2004,1,1), T, 'k')
    figure(20083)
    plot(jd - julian(2004,1,1), P, 'k')
    axis ij
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2g. LOOP FOR RAPID 6 (D344 --> D359)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Notes:
%  ------
%


disp('---------  RAPID 6 (D344 --> D359) ---------')
[mooring6, sn6, mc6, z6, lon6] = textread([boundarydir, 'eb_rapid_6.dat'],'%s %d %d %d %6.2f','headerlines ',5);

T6 = zeros(length(sn6), length(JG));
P6 = zeros(length(sn6), length(JG));
S6 = zeros(length(sn6), length(JG));


i = 1; j = 1;
for i = 1: length(sn6)
    
    infile = [hydrodir, mooring6{i,:}, '_grid.mat'];
    load(infile);
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jd(1:end-1), Sfs(mc6(i),1:end-1)', JG)';
    temp        = interp1(jd(1:end-1), Tfs(mc6(i),1:end-1)', JG)';
    pressure    = interp1(jd(1:end-1), Pfs(mc6(i),1:end-1)', JG)';
    Sfs6(j,:) = salinity;
    Pfs6(j,:) = pressure;
    Tfs6(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(20091)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs6(i, :), col{:,i},'LineWidth',2);
        title('RAPID 6 - SALINITY')
        
        figure(20092)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs6(i, :), col{:,i},'LineWidth',2);
        title('RAPID 6 - TEMPERATURE')
        
        figure(20093)
        hold on; box on;
        plot(JG - julian(2004,1,1), Pfs6(i, :), col{:,i},'LineWidth',2);
        title('RAPID 6 - PRESSURE')
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2h. LOOP FOR RAPID 7 (D359 --> JC064)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------
%
disp('---------  RAPID 7 (D359 --> JC064) ---------')
[mooring7, sn7, mc7, z7, lon7] = textread([boundarydir, 'eb_rapid_7.dat'],'%s %d %d %d %6.2f','headerlines ',5);

T7 = zeros(length(sn7), length(JG));
P7 = zeros(length(sn7), length(JG));
S7 = zeros(length(sn7), length(JG));

i = 1; j = 1;
for i = 1: length(sn7)
    
    infile = [hydrodir, mooring7{i,:}, '_grid.mat'];
    load(infile);
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jd(1:end-1), Sfs(mc7(i),1:end-1)', JG)';
    temp        = interp1(jd(1:end-1), Tfs(mc7(i),1:end-1)', JG)';
    pressure    = interp1(jd(1:end-1), Pfs(mc7(i),1:end-1)', JG)';
    Sfs7(j,:) = salinity;
    Pfs7(j,:) = pressure;
    Tfs7(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(20091)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs7(i, :), col{:,i},'LineWidth',2);
        title('RAPID 7 - SALINITY')
        
        figure(20092)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs7(i, :), col{:,i},'LineWidth',2);
        title('RAPID 7 - TEMPERATURE')
        
        figure(20093)
        hold on; box on;
        plot(JG - julian(2004,1,1), Pfs7(i, :), col{:,i},'LineWidth',2);
        title('RAPID 7 - PRESSURE')
    end
    
end


%%%%%WHAT FILE SHOULD BE LOADED IN HERE?
% load([basedir, 'rapid/data/amoc/grdat/eb_2008_grid_merge.mat']);
% jd6 = jd; SGfs6 = SGfs; TGfs6 = TGfs; PG6 = p_grid;
% if mc_check_plot
% figure(20091)
% plot(jd - julian(2004,1,1), S, 'k')
% figure(20092)
% plot(jd - julian(2004,1,1), T, 'k')
% figure(20093)
% plot(jd - julian(2004,1,1), P, 'k')
% axis ij
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2i. LOOP FOR RAPID 8 (JC064 --> Di382)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------
%
disp('---------  RAPID 8 (JC064 --> Di382) ---------')

[mooring8, sn8, mc8, z8, lon8] = textread([boundarydir, 'eb_rapid_8.dat'],'%s %d %d %d %6.2f','headerlines ',5);

T8 = zeros(length(sn8), length(JG));
P8 = zeros(length(sn8), length(JG));
S8 = zeros(length(sn8), length(JG));

i = 1; j = 1;
for i = 1: length(sn8)
    
    infile = [hydrodir, mooring8{i,:}, '_grid.mat'];
    load(infile);
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jd(1:end-1), Sfs(mc8(i),1:end-1)', JG)';
    temp        = interp1(jd(1:end-1), Tfs(mc8(i),1:end-1)', JG)';
    pressure    = interp1(jd(1:end-1), Pfs(mc8(i),1:end-1)', JG)';
    Sfs8(j,:) = salinity;
    Pfs8(j,:) = pressure;
    Tfs8(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(20091)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs7(i, :), col{:,i},'LineWidth',2);
        title('RAPID 8 - SALINITY')
        
        figure(20092)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs7(i, :), col{:,i},'LineWidth',2);
        title('RAPID 8 - TEMPERATURE')
        
        figure(20093)
        hold on; box on;
        plot(JG - julian(2004,1,1), Pfs7(i, :), col{:,i},'LineWidth',2);
        title('RAPID 8 - PRESSURE')
    end
    
end

s_date8         = gregorian(jd(1));
s_month8        = month(s_date8(2),:);
e_date8        = gregorian(jd(end));
e_month8        = month(e_date8(2),:);
disp('----------------------------');
disp('RAPID 8');
disp('  ');
disp(['START DATE  = ', num2str(s_date8(3)),' ', s_month8, ' ', num2str(s_date8(1))]);
disp(['END DATE    = ', num2str(e_date8(3)),' ', e_month8, ' ', num2str(e_date8(1))]);
disp('----------------------------');
disp('  ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2j. LOOP FOR RAPID 9 ( Di382 --> JC103)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------
%
disp('---------  RAPID 9 ( Di382 --> JC103 ) ---------')

[mooring9, sn9, mc9, z9, lon9] = textread([boundarydir, 'eb_rapid_9.dat'],'%s %d %d %d %6.2f','headerlines ',5);

T9 = zeros(length(sn9), length(JG));
P9 = zeros(length(sn9), length(JG));
S9 = zeros(length(sn9), length(JG));


i = 1; j = 1;
for i = 1: length(sn9)

    infile = [hydrodir, mooring9{i,:}, '_grid.mat'];
    load(infile);
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jd(1:end-1), Sfs(mc9(i),1:end-1)', JG)';
    temp        = interp1(jd(1:end-1), Tfs(mc9(i),1:end-1)', JG)';
    pressure    = interp1(jd(1:end-1), Pfs(mc9(i),1:end-1)', JG)';
    Sfs9(j,:) = salinity;
    Pfs9(j,:) = pressure;
    Tfs9(j,:) = temp;
%
    j = j + 1;

    if mc_check_plot
        figure(20091)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs9(i, :), col{:,i},'LineWidth',2);
        title('RAPID 9 - SALINITY')

        figure(20092)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs9(i, :), col{:,i},'LineWidth',2);
        title('RAPID 9 - TEMPERATURE')

        figure(20093)
        hold on; box on;
        plot(JG - julian(2004,1,1), Pfs9(i, :), col{:,i},'LineWidth',2);
        title('RAPID 9 - PRESSURE')
    end

end

%display(' NO calibrated RAPID 9 data yet. using hydro_grid_step1')
%[Sfs9, Tfs9, Pfs9, jd] = hydro_grid_step1(mooring9,sn9);
%
%    Sfs9    = interp1(jd(1:end-1), Sfs9(:,1:end-1)', JG)';
%    Tfs9    = interp1(jd(1:end-1), Tfs9(:,1:end-1)', JG)';
%    Pfs9    = interp1(jd(1:end-1), Pfs9(:,1:end-1)', JG)';

s_date9         = gregorian(jd(1));
s_month9        = month(s_date9(2),:);
e_date9        = gregorian(jd(end));
e_month9        = month(e_date9(2),:);
disp('----------------------------');
disp('RAPID 9');
disp('  ');
disp(['START DATE  = ', num2str(s_date9(3)),' ', s_month9, ' ', num2str(s_date9(1))]);
disp(['END DATE    = ', num2str(e_date9(3)),' ', e_month9, ' ', num2str(e_date9(1))]);
disp('----------------------------');
disp('  ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2k. LOOP FOR RAPID 10 ( JC103 --> DY039)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------
%
disp('---------  RAPID 10 ( JC103 --> DY039) ---------')

[mooring10, sn10, mc10, z10, lon10] = textread([boundarydir, 'eb_rapid_10.dat'],'%s %d %d %d %6.2f','headerlines ',5);

T10 = zeros(length(sn10), length(JG));
P10 = zeros(length(sn10), length(JG));
S10 = zeros(length(sn10), length(JG));

i = 1; j = 1;
for i = 1: length(sn10)

    infile = [hydrodir, mooring10{i,:}, '_grid.mat'];
    load(infile);
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jd(1:end-1), Sfs(mc10(i),1:end-1)', JG)';
    temp        = interp1(jd(1:end-1), Tfs(mc10(i),1:end-1)', JG)';
    pressure    = interp1(jd(1:end-1), Pfs(mc10(i),1:end-1)', JG)';
    Sfs10(j,:) = salinity;
    Pfs10(j,:) = pressure;
    Tfs10(j,:) = temp;
%
    j = j + 1;

    if mc_check_plot
        figure(20091)
        hold on; box on;
        plot(JG - julian(2004,1,1), Sfs10(i, :), col{:,i},'LineWidth',2);
        title('RAPID 10 - SALINITY')

        figure(20092)
        hold on; box on;
        plot(JG - julian(2004,1,1), Tfs10(i, :), col{:,i},'LineWidth',2);
        title('RAPID 10 - TEMPERATURE')

        figure(20093)
        hold on; box on;
        plot(JG - julian(2004,1,1), Pfs10(i, :), col{:,i},'LineWidth',2);
        title('RAPID 10 - PRESSURE')
    end

end

%display(' NO calibrated RAPID 9 data yet. using hydro_grid_step1')
%[Sfs10, Tfs10, Pfs10, jd] = hydro_grid_step1(mooring10,sn10);
%
%    Sfs10   = interp1(jd(1:end-1), Sfs10(:,1:end-1)', JG)';
%    Tfs10   = interp1(jd(1:end-1), Tfs10(:,1:end-1)', JG)';
%    Pfs10   = interp1(jd(1:end-1), Pfs10(:,1:end-1)', JG)';

s_date10         = gregorian(jd(1));
s_month10        = month(s_date10(2),:);
e_date10         = gregorian(jd(end));
e_month10        = month(e_date10(2),:);
disp('----------------------------');
disp('RAPID 10');
disp('  ');
disp(['START DATE  = ', num2str(s_date10(3)),' ', s_month10, ' ', num2str(s_date10(1))]);
disp(['END DATE    = ', num2str(e_date10(3)),' ', e_month10, ' ', num2str(e_date10(1))]);
disp('----------------------------');
disp('  ');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  3.  CONCATENATE AND ORDER THE MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  This step adds alll the matrices for the deployments together to form
%  large data sets.  These are then sorted into pressure order at every
%  time step before the gridding takes place.

%if mc_check_plot
    figure(1)   %  graph of the data to show that it is all there!
    clf;
    subplot(2,1,1);
    hold on; box on;
    plot(JG - julian(2004,1,1), Tfs1, 'k.')
    plot(JG - julian(2004,1,1), Tfs2, 'b.')
    plot(JG - julian(2004,1,1), Tfs3a, 'r.')
    plot(JG - julian(2004,1,1), Tfs3b, 'y.')
    plot(JG - julian(2004,1,1), Tfs4, 'b.')
    plot(JG - julian(2004,1,1), Tfs5, 'r.')
    plot(JG - julian(2004,1,1), Tfs6, 'g.')
    plot(JG - julian(2004,1,1), Tfs7, 'k.')
    plot(JG - julian(2004,1,1), Tfs8, 'b.')
    plot(JG - julian(2004,1,1), Tfs9, 'r.')
    plot(JG - julian(2004,1,1), Tfs10, 'g.')
    ylabel('TEMPERATURE')
    title('QUICK CHECK OF DATA')
    timeaxis([2004,1,1]);
    subplot(2,1,2);
    hold on; box on;
    plot(JG - julian(2004,1,1), Sfs1, 'k.')
    plot(JG - julian(2004,1,1), Sfs2, 'b.')
    plot(JG - julian(2004,1,1), Sfs3a, 'r.')
    plot(JG - julian(2004,1,1), Sfs3b, 'y.')
    plot(JG - julian(2004,1,1), Sfs4, 'b.')
    plot(JG - julian(2004,1,1), Sfs5, 'r.')
    plot(JG - julian(2004,1,1), Sfs6, 'g.')
    plot(JG - julian(2004,1,1), Sfs7, 'k.')
    plot(JG - julian(2004,1,1), Sfs8, 'b.')
    plot(JG - julian(2004,1,1), Sfs9, 'r.')
    plot(JG - julian(2004,1,1), Sfs10, 'g.')
    ylabel('SALINITY')
    title('QUICK CHECK OF DATA')
    timeaxis([2004,1,1]);
%end

% all the matrices for the deployments stacked together
Pfs     = [Pfs1;Pfs2;Pfs3a;Pfs3b;Pfs4;Pfs5;Pfs6;Pfs7;Pfs8;Pfs9;Pfs10];
Sfs     = [Sfs1;Sfs2;Sfs3a;Sfs3b;Sfs4;Sfs5;Sfs6;Sfs7;Sfs8;Sfs9;Sfs10];
Tfs     = [Tfs1;Tfs2;Tfs3a;Tfs3b;Tfs4;Tfs5;Tfs6;Tfs7;Tfs8;Tfs9;Tfs10];

% order the matrices at every time step to avoid too many NaNs creeping in
% 2004 removed....
P_sort = NaN .* ones(size(Pfs)); T_sort = NaN .* ones(size(Tfs)); S_sort = NaN .* ones(size(Sfs));
j = 1;
for ii = 1: length(JG)
    [P_variable, ix] = sort(Pfs(:, ii));
    P_sort(:,j) = Pfs(ix,ii);
    T_sort(:,j) = Tfs(ix,ii);
    S_sort(:,j) = Sfs(ix,ii);
    j = j + 1;
end

% removing unused rows of the sorted matrices

i = 1; j = 1;
for i = 1: length(P_sort(:,1))
    ix = find(isnan(P_sort(i,:)));
    if size(ix) < length(JG)
        Pfss(j,:) = P_sort(i, :);
        Tfss(j,:) = T_sort(i, :);
        Sfss(j,:) = S_sort(i, :);
        j = j + 1;
    end
end

clear Pfs1 Pfs2 Pfs3a Pfs3b Pfs4 Pfs5 Pfs6 Pfs7 Pfs8 Pfs9 Pfs10;
clear Sfs1 Sfs2 Sfs3a Sfs3b Sfs4 Sfs5 Sfs6 Sfs7 Sfs8 Sfs9 Sfs10;
clear Tfs1 Tfs2 Tfs3a Tfs3b Tfs4 Tfs5 Tfs6 Tfs7 Tfs8 Tfs9 Tfs10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  4.  GRIDDING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hydro-Gridding....

%  This is the most contentious part of this code.  I have edited the
%  gridding function con_tprof0.m to select the pressure grid at every time
%  step rather than just take the the common range throught out the whole
%  time series.  This then carries out the gridding at every 12 hour time
%  step based on the merged data matrix.  I also forced the gridding only
%  to run when there are at least 4 data points.  See notes on the griddig
%  in the document written on this project!
% Gridding is the process of vertically interpolating the T, S and P data onto a
% regular pressure grid.  It uses the ds/dp and dt/dp climatologies as a
% basis for the weighted integration.  Zoli is working on producing a
% somewhat better eastern boundary climatlogy as part of his PIES work
% (ref: D344 cruise report)

%  The climatology used is the existing eb_slope ones created by TOK.  I
%  have also experimented with using a mean climatology based on the 5
%  years of mooring data - i.e. using the rough 'climatology' from TOK as a
%  first iteration before re-gridding using the 5 year mean.

if gridding
    pmin = ceil(min(Pfss)/20) * 20;
    pmax = floor(max(Pfss)/20) * 20;
    
    if strcmp(TS_CLIMATOLOGY_TP,'annual')
        display('con_tprof0_annual')
        [TGfs, SGfs] = con_tprof0_annual(Tfss, Sfss, Pfss, pg', 0*ones(1, size(Tfss, 2)), int_step, TSclim, ...
                                          TS_CLIMATOLOGY,TS_CLIMATOLOGY_NAME);
    elseif strcmp(TS_CLIMATOLOGY_TP,'seasonal')
        display('con_tprof0_monthly')
        [TGfs, SGfs] = con_tprof0_monthly(Tfss, Sfss, Pfss, pg', GTV(:,2), int_step, TSclim,TS_CLIMATOLOGY, ...
                                          TS_CLIMATOLOGY_NAME );
    end
    
    % time is passed in as 0 in the above line?!?
    % pass the month of measurement into con_tprof0, instead of the empty time vector.
    %
    %	Tfss -- timeseries of discrete temperature measurements from the moorings (e.g. 25x5001)
    %       Tfss -- timeseries of discrete temperature measurements from the moorings.
    %       Tfss -- timeseries of discrete temperature measurements from the moorings.
    %       pg   -- pressure grid onto which the profiles will be integrated
    %       GTV  -- months when each measurement was made.
    %	int_step -- integration step size [dbar] between grid points
    %	TSclim   = path to the file containing the dT/dP and dS/dP climatology
    
    %%%%%%%%%%%%%%
    % remove the data gap from 2 Feb 2006 to 22 May 2006 (this will be
    % interpolated over prior to use by the MOC code.)
    % leaving the gap slightly smaller as implied by the existing data has
    % a large effect on the final result.  :(
    %  !!!   This requires looking into when the time series is next
    %  updated !!)
    ix = find(JG > julian(2006,2,15) & JG < julian(2006, 5, 16));
    SGfs(:,ix) = NaN;
    TGfs(:,ix) = NaN;
    
    display(['saving: EB_merged_' TS_CLIMATOLOGY_TP ' ' TS_CLIMATOLOGY '_2012.mat'])
    % Alloacte variables into a structure
    EB_merg.JG     = JG;
    EB_merg.pg     = pg;
    EB_merg.Tfs    = Tfs;
    EB_merg.Sfs    = Sfs;
    EB_merg.Pfs    = Pfs;
    EB_merg.TGfs   = TGfs;
    EB_merg.SGfs   = SGfs;
    EB_merg.P_sort = P_sort;
    EB_merg.T_sort = T_sort;
    EB_merg.S_sort = S_sort;
    EB_merg.TS_CLIMATOLOGY   = TS_CLIMATOLOGY;
    EB_merg.TS_CLIMATOLOGY_TP= TS_CLIMATOLOGY_TP;
    EB_merg.TS_CLIMATOLOGY_NAME= TS_CLIMATOLOGY_NAME;
    EB_merg.clim_file        = TSclim;
    EB_merg.OUT_FILE         = OUT_FILE;
    EB_merg.EB_creation_date = datestr(now);
    EB_merg.MERG_REVISION    = REV(6:8);            % store the revision number of this script
    EB_merg.MERG_AUTHOR      = REV_AUTHOR(9:end-1); % store the revision author of this script
    EB_merg.MERG_DATE        = REV_DATE(7:end-1);   % store the revision date of this script
    EB_merg.EB_path = grdatdir; % path to output file
    EB_merg.function_name = function_name; % the name and path to this function
    EB_merg.interpolation_depth=idepth-20;
    EB_merg.matlab_version = version;
    
% bim September 2014
% load in climatology to save in the data file
% read in and save as a structure
% load in climatology
	eval([' load ' TSclim])
% convert to structure

	eval(['matObj = matfile( '''  TSclim ''')' ])
        info=whos(matObj);
	for kk=1:length(info);
	  eval([ 'clim.' info(kk).name ' = ' info(kk).name])
	end
    EB_merg.climatology = clim; % structure of the climatology variables used to grid the data

    EB_merg
		
        
    
    % save data
    % original non-structure file format.
    % unused in later versions
    %     ['save ' grdatdir char(OUT_FILE(1)) ...
    %           ' JG pg Tfs Sfs Pfs TGfs SGfs P_sort T_sort S_sort']
    %     eval(['save ' grdatdir char(OUT_FILE(1)) ...
    %           ' JG pg Tfs Sfs Pfs TGfs SGfs P_sort T_sort S_sort']);
    
    ['save ' grdatdir char(OUT_FILE(1)) ' EB_merg']
    eval(['save ' grdatdir char(OUT_FILE(1)) ' EB_merg']);
    
    % interpolating over the shorter gaps in the data - but leaving the
    % larger ones close to the surface as they will be extrapolated once
    % the transports have been calculated by the MOC code.  The idea is
    % then that this product can be imported into the MOC code
    
    % A problem is that it wants to level out the gaps in the surface layer
    % that are best treated by extrapolating the derived transports in the
    % MOC code...  temp fix by ignoring the upper layer...
    
    %clear all; close all
    %clear
    %    load ../mat_files/EB_merged_data_no_interp_2010.mat
    %eval(['load ../mat_files/EB_merged_annual_' TS_CLIMATOLOGY '_2010.mat'])
    %JG -- julian day
    %pg -- pressure grid
    %Tfs -- original stacked Temperature data from the deployments (144x5001)
    %Sfs -- original stacked salinity data from the deployments (144x5001)
    %Pfs -- original stacked Pressure data from the deployments (144x5001)
    %TGfs -- temperature interpolated onto the pressure grid (PG)
    %SGfs -- salinity interpolated onto the pressure grid (PG)
    %P_sort -- Pfs sorted on pressure
    %T_sort -- Tfs sorted on pressure
    %S_sort -- Sfs sorted on pressure
    
    stddy_tol  = 3.5; % 4
    std_win    = 3.5; % 3.5 * std of the time series
    [nloop]    = 3; % 5
    graphics   = 'y';
    temp = []; salinity = [];
    
    %%%%%%
    %  despike time series
    %%%%%%
    for i = 1 : length(TGfs(:,1))   % loop through each depth level
        
        [temp(i,:),dx,ndx] = ddspike(TGfs(i,:),[-std_win*nanstd(TGfs(i,:)),...
            std_win*nanstd(TGfs(i,:))],stddy_tol,[nloop],'y',NaN);
        [salinity(i,:),dx,ndx] = ddspike(SGfs(i,:),[-std_win*nanstd(SGfs(i,:)),...
            std_win*nanstd(SGfs(i,:))],stddy_tol,[nloop],'y',NaN);
    end
    
    [m,n] = size(TGfs);
    TG_east = NaN * ones(m,n); SG_east = NaN * ones(m,n);
    
    % GDM, 9/4/2013
    % Need to interpolate horizontally over nans but...
    % Don't want to create fake values above knocked down moorings
    % So select a depth below which, interpoation happens
    
%    idepth = 220; % multiples of 20
    I = find(pg == idepth);
    
    % copy the top idepth into the new file
    % no temporal interpolation
    SG_east([1:I-1],:)=SGfs([1:I-1],:);
    TG_east([1:I-1],:)=TGfs([1:I-1],:);
    
    %    i = 1; j = 1;
    %    for i = 1: 6 % top 100m based on a grid of 20dbar
    %        SG_east(j,:) = SGfs(i, :);
    %        TG_east(j,:) = TGfs(i, :);
    %        j = j + 1;
    %    end
%
    
    % this piece of code has the ability to create spurious values of T & S ** tk
    % how large are the temporal gaps in the time series ?
    % linear should only be used for small gaps.
    
    %     i = 7; j = 7; % PG(7) = 120;
    %     i = 11; j = 11; % PG(11) = 120;

    i = I; j = I;
    for i = I: length(pg) % for each depth
        
        % locate all non nan values in the despiked temperature
        it = find(~isnan(temp(i,:)));
        % locate all non nan values in the despiked salinity
        is = find(~isnan(salinity(i,:)));
        % interpolate in time over the missing data
        SG_east(j,:) = interp1(JG(is), salinity(i, is), JG);
        % interpolate in time over the missing data
        TG_east(j,:) = interp1(JG(it), temp(i, it), JG);
        
        j = j + 1;
    end
    
    figure;
    % spurious values created in the interpolation?
    hold on; grid on; box on;
    plot(SG_east,TG_east,'k')
    plot(SGfs,TGfs,'r--')
    
    
    close all
    figure(1); clf
    subplot(3,1,1)
    contourf(JG - julian(2004,1,1), pg, TGfs); axis ij
    timeaxis([2004,1,1]); ylabel('TEMPERATURE')
    title('BEFORE DESPIKING AND INTEROLATION')
    subplot(3,1,2)
    contourf(JG - julian(2004,1,1), pg, temp); axis ij
    timeaxis([2004,1,1]); ylabel('TEMPERATURE');
    title('AFTER DESPIKING')
    subplot(3,1,3)
    contourf(JG - julian(2004,1,1), pg, TG_east); axis ij
    timeaxis([2004,1,1]); ylabel('TEMPERATURE')
    title('AFTER DESPIKING AND INTERPOLATION')
    
    figure(2);clf
    subplot(3,1,1)
    contourf(JG - julian(2004,1,1), pg, SGfs); axis ij
    timeaxis([2004,1,1]); ylabel('SALINITY')
    title('BEFORE DESPIKING AND INTEROLATION')
    subplot(3,1,2)
    contourf(JG - julian(2004,1,1), pg, salinity); axis ij
    timeaxis([2004,1,1]); ylabel('SALINITY')
    title('AFTER DESPIKING')
    subplot(3,1,3)
    contourf(JG - julian(2004,1,1), pg, SG_east); axis ij
    timeaxis([2004,1,1]); ylabel('SALINITY')
    title('AFTER DESPIKING AND INTERPOLATION')
    
    % Alloacte variables into a structure
    EB_merg.TG_east= TG_east;
    EB_merg.SG_east= SG_east;
    
    %    ['save ' grdatdir char(OUT_FILE(3)) ...
    %           ' JG pg Tfs Sfs Pfs TGfs SGfs TG_east SG_east P_sort T_sort S_sort']
    %    eval(['save ' grdatdir char(OUT_FILE(3)) ...
    %           ' JG pg Tfs Sfs Pfs TGfs SGfs TG_east SG_east P_sort T_sort S_sort'])
    
    ['save ' grdatdir char(OUT_FILE(2)) ' EB_merg']
    eval(['save ' grdatdir char(OUT_FILE(2)) ' EB_merg'])
    
    EB_merg
    
    %	JG -- julian day
    %	pg -- pressure grid
    %	Tfs -- original stacked Temperature data from the deployments (144x5001)
    %	Sfs -- original stacked salinity data from the deployments (144x5001)
    %	Pfs -- original stacked Pressure data from the deployments (144x5001)
    %	TG_east -- temperature interpolated onto the time grid (JG) after despiking
    %	SG_east -- salinity interpolated onto the time grid (JG) after despiking
    %	P_sort -- Pfs sorted on pressure. no despike. no temporal interp
    %	T_sort -- Tfs sorted on pressure. no despike. no temporal interp
    %	S_sort -- Sfs sorted on pressure. no despike. no temporal interp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wb_rapid_3b.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  5.  PLOTTING THE BATHYMETRIC CHARTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  charts to show the positions and serial numbers of the MicroCATS used
%  for the eastern boundary merging.  The bathymetry data comes from
%  Discovery cruise D279.


if bathy
    
    load(bathyfile);
    
    figure(31)
    clf
    hold on; box on;
    plot(lon, z, 'k')
    axis ij
    xlim([-25, -13]);
    xlabel('LONGITUDE','FontSize',14)
    ylabel('DEPTH [m]','FontSize',14)
    plot(lon1, z1, 'r*')
    text(lon1-0.2, z1, '\rightarrow')
    text(lon1-0.5, z1, num2str(sn1))
    title('RAPID 1 (D277 --> CD170)','FontSize',14)
    
    figure(32)
    clf
    hold on; box on;
    plot(lon, z, 'k')
    axis ij
    xlim([-25, -13]);
    xlabel('LONGITUDE','FontSize',14)
    ylabel('DEPTH [m]','FontSize',14)
    plot(lon2, z2, 'r*')
    text(lon2-0.2, z2, '\rightarrow')
    text(lon2-0.5, z2, num2str(sn2))
    title('RAPID 2 (CD170/CD177 --> D304)','FontSize',14)
    
    figure(331)
    clf
    hold on; box on;
    plot(lon, z, 'k')
    axis ij
    xlim([-16, -13]);
    xlabel('LONGITUDE','FontSize',14)
    ylabel('DEPTH [m]','FontSize',14)
    plot(lon3a, z3a, 'r*')
    text(lon3a-0.2, z3a, '\rightarrow')
    text(lon3a-0.5, z3a, num2str(sn3a))
    title('RAPID 3a (D304 --> F343) ZOOM','FontSize',14)
    
    figure(33)
    clf
    hold on; box on;
    plot(lon3a, z3a, 'r*')
    plot(lon, z, 'k')
    axis ij
    xlim([-25, -13]);
    xlabel('LONGITUDE','FontSize',14)
    ylabel('DEPTH [m]','FontSize',14)
    text(lon3a-0.2, z3a, '\rightarrow')
    text(lon3a-0.5, z3a, num2str(sn3a))
    title('RAPID 3a (D304 --> F343)','FontSize',14)
    
    figure(332)
    clf
    hold on; box on;
    plot(lon, z, 'k')
    axis ij
    xlim([-25, -13]);
    xlabel('LONGITUDE','FontSize',14)
    ylabel('DEPTH [m]','FontSize',14)
    plot(lon3b, z3b, 'r*')
    text(lon3b-0.2, z3b, '\rightarrow')
    text(lon3b-0.5, z3b, num2str(sn3b))
    
    plot(lon3a(11:21), z3a(11:21), 'b*')
    text(lon3a(11:21)-0.2, z3a(11:21), '\rightarrow')
    text(lon3a(11:21)-0.5, z3a(11:21), num2str(sn3a(11:21)))
    title('RAPID 3b (F343 --> D274) nb EBH0 to EBH3 carried from RAPID 3a ','FontSize',14)
    
    figure(34)
    clf
    hold on; box on;
    plot(lon, z, 'k')
    axis ij
    xlim([-25, -13]);
    xlabel('LONGITUDE','FontSize',14)
    ylabel('DEPTH [m]','FontSize',14)
    plot(lon4, z4, 'r*')
    text(lon4-0.2, z4, '\rightarrow')
    text(lon4-0.5, z4, num2str(sn4))
    title('RAPID 4 ( D324 --> D334)','FontSize',14)
    
    figure(35)
    clf
    hold on; box on;
    plot(lon, z, 'k')
    axis ij
    xlim([-25, -13]);
    xlabel('LONGITUDE','FontSize',14)
    ylabel('DEPTH [m]','FontSize',14)
    plot(lon5, z5, 'r*')
    text(lon5-0.2, z5, '\rightarrow')
    text(lon5-0.5, z5, num2str(sn5))
    title('RAPID 5 ( D334 --> D344)','FontSize',14)
    
    figure(36)
    clf
    hold on; box on;
    plot(lon, z, 'k')
    axis ij
    xlim([-25, -13]);
    xlabel('LONGITUDE','FontSize',14)
    ylabel('DEPTH [m]','FontSize',14)
    plot(lon6, z6, 'r*')
    text(lon6-0.2, z6, '\rightarrow')
    text(lon6-0.5, z6, num2str(sn6))
    title('RAPID 6 ( D344 --> D359)','FontSize',14)
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  6.  PLOTTING THE GRIDDED AND MERGED PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if gridding == 0
    load ../mat_files/EB_merged_data_no_interp_2010.mat
end

figure(1000)
clf
subplot(3,1,1)
size(jd1); size(PG1); size(SGfs1); size(TGfs1);
contourf(jd1 - julian(2004,1,1), PG1, SGfs1)
axis ij
hold on
contourf(jd2 - julian(2004,1,1), PG2, SGfs2)
contourf(jd3a - julian(2004,1,1), PG3a, SGfs3a)
contourf(jd3b - julian(2004,1,1), PG3b, SGfs3b)
contourf(jd4 - julian(2004,1,1), PG4, SGfs4)
contourf(jd5 - julian(2004,1,1), PG5, SGfs5)
%contourf(jd6 - julian(2004,1,1), PG6, SGfs6)
%contourf(jd7 - julian(2004,1,1), PG7, SGfs7)
%contourf(jd8 - julian(2004,1,1), PG8, SGfs8)
ylim([0,5000])
xlim([jd1(1) - julian(2004,1,1), JG(end) - julian(2004,1,1)])
title('CURRENT MERGING - SALINITY')
timeaxis([2004,1,1]);
subplot(3,1,2)
contourf(JG - julian(2004,1,1), pg, SGfs)
axis ij
ylim([0,5000])
xlim([jd1(1) - julian(2004,1,1), JG(end) - julian(2004,1,1)])
title('NEW MERGING BEFORE INTERPOLATION - SALINITY')
timeaxis([2004,1,1]);

figure(1001)
clf
subplot(3,1,1)
contourf(jd1 - julian(2004,1,1), PG1, TGfs1)
axis ij
hold on
contourf(jd2 - julian(2004,1,1), PG2, TGfs2)
contourf(jd3a - julian(2004,1,1), PG3a, TGfs3a)
contourf(jd3b - julian(2004,1,1), PG3b, TGfs3b)
contourf(jd4 - julian(2004,1,1), PG4, TGfs4)
contourf(jd5 - julian(2004,1,1), PG5, TGfs5)
%contourf(jd6 - julian(2004,1,1), PG6, TGfs6)
%contourf(jd7 - julian(2004,1,1), PG7, TGfs7)
%contourf(jd8 - julian(2004,1,1), PG8, TGfs8)
ylim([0,5000])
xlim([jd1(1) - julian(2004,1,1), JG(end) - julian(2004,1,1)])
title('CURRENT MERGING - TEMPERATURE')
timeaxis([2004,1,1]);
subplot(3,1,2)
contourf(JG - julian(2004,1,1), pg, TGfs)
axis ij
ylim([0,5000])
xlim([jd1(1) - julian(2004,1,1), JG(end) - julian(2004,1,1)])
title('NEW MERGING BEFORE INTERPOLATION- TEMPERATURE')
timeaxis([2004,1,1]);


if gridding == 0
    load ../mat_files/EB_merged_data_2010.mat
end

figure(1000)
subplot(3,1,3)
contourf(JG - julian(2004,1,1), pg, SG_east)
axis ij
ylim([0,5000])
xlim([jd1(1) - julian(2004,1,1), JG(end) - julian(2004,1,1)])
title('NEW MERGING (after interp - SALINITY')
timeaxis([2004,1,1]);

figure(1001)
subplot(3,1,3)
contourf(JG - julian(2004,1,1), pg, TG_east)
axis ij
ylim([0,5000])
xlim([jd1(1) - julian(2004,1,1), JG(end) - julian(2004,1,1)])
title('NEW MERGING (after interp - TEMPERATURE')
timeaxis([2004,1,1]);



