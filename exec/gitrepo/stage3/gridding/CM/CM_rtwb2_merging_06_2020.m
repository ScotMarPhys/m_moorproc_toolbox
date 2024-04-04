%%  Code for the Merging of the Western Boundary 2
%%  Current Meter Data
%
%  Code History
% -------------
% Jan 2017 - Adapted for OSNAP moorings (by L. Houpert)
% -------------
% -------------
% Mar 2024 - Started to clean code for current meter processing. Contained a
% lot of microcat processing which is not applicalbe so I deleted any 
% irrelvant lines of code and comments.
% Work in progress. (K. Burmeister)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PARAMETER PRAEMBLE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% Need to be updated for:
% - directories and file names
% - add additional deployment periods (add moorX = 'filename')
% - turn off/on checkplots
% - edit start/end of total time period
% - data version
% - set depth of shallowest instrument (idepth)
% - praembel for despiking
% - if you add a deployment period, you need to edit steps 1-3
%---------

close all
MOOR = 'RTWB2';
moor1 = 'CM_rtwb2_osnap_01_2014';
moor2 = 'CM_rtwb2_osnap_02_2015';
moor3 = 'CM_rtwb2_osnap_03_2016';
moor4 = 'CM_rtwb2_osnap_04_2017';
moor5 = 'CM_rtwb2_osnap_05_2018';
moor6 = 'CM_rtwb2_osnap_06_2020';

% in- and output directories
basedir  = pathosnap;
hydrodir = [basedir '/data/moor/proc/velocity_grid/'];
grdatdir = [pathgit '\data\processed\stage3_gridding_CM\'];%[basedir '/data/moor/proc/velocity_grid_merged/'];
boundarydir = [execdir 'gitrepo/stage3/gridding/CM/'];

% turns on/off check plots. off=false, on=true
cm_check_plot = false ;  

% add start and end of total time period
jg_start                = datenum(2014,07,01,00,00,00);
jg_end                  = datenum(2022,07,31,00,00,00);

jg_start_str = datestr(jg_start,'YYYYmm');
jg_end_str = datestr(jg_end,'YYYYmm');

% put out data version?
data_version = 'v0'; 

JG = jg_start: 0.5: jg_end; % full time series using 2 samples per day
pgg = 0:20:2000; % depths in 20dbar bins
idepth = 40; % depth of shallowest instrument, steps of 20 (i.e. if shallowest
             % instrument is at 50m, select 40m

% output filenames (should be automatic)
outputfile = [MOOR '_merg_linear_interp' jg_start_str '_' jg_end_str];
outputfile_stacked = ['CM_rtwb2_stacked_' jg_start_str '_' jg_end_str,...
                        '_' data_version];

% Praemble for 4.b despike time series
stddy_tol  = 10; % set max velocities renage (median pm stddy_tol)
std_win    = 3.5; % 3.5 * std of the time series (median pm std_wind)
[nloop]    = 5; % max number of despiking repetitions
graphics   = 'y'; %y if plots of despiking on, 'n' if off
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  LOAD AND VERTICALLY STACK CM FILES FOR EACH CRUISE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes: 
%  This need to be updated if you add a deployment period
%
%  KB 28 Mar 2023 - wrote a function to clean up script and ensure
%  changes are automatically propagated between RTEB, RTWB1, and RTWB2
%  scripts
%  ------

disp('---------  OSNAP 1 (KN221 --> PE399) ---------')
[Ufs1,Vfs1,Wfs1,Pfs1] = load_and_stacked_CM(boundarydir,hydrodir,moor1,JG,cm_check_plot);

disp('---------  OSNAP 2 (PE399 --> DY053) ---------')
[Ufs2,Vfs2,Wfs2,Pfs2] = load_and_stacked_CM(boundarydir,hydrodir,moor2,JG,cm_check_plot);

disp('---------  OSNAP 3 (DY053 --> DY078) ---------')
[Ufs3,Vfs3,Wfs3,Pfs3] = load_and_stacked_CM(boundarydir,hydrodir,moor3,JG,cm_check_plot);

disp('---------  OSNAP 4 (DY078 --> AR30) ---------')
[Ufs4,Vfs4,Wfs4,Pfs4] = load_and_stacked_CM(boundarydir,hydrodir,moor4,JG,cm_check_plot);

disp('---------  OSNAP 5 (AR30 --> DY120) ---------')
[Ufs5,Vfs5,Wfs5,Pfs5] = load_and_stacked_CM(boundarydir,hydrodir,moor5,JG,cm_check_plot);

disp('---------  OSNAP 6 (DY120 --> JC238) ---------')
[Ufs6,Vfs6,Wfs6,Pfs6] = load_and_stacked_CM(boundarydir,hydrodir,moor6,JG,cm_check_plot);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.  CONCATENATE MATRICES ALONG TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  This needs to be updated if you add a deployment period. Update also
%  internal function stack_vars() at the end of the script
%
%  This step concatenates all data along time axis. 
%
% KB Mar 24: Output stacked files with only 6 depths, need to output z in
% load_and_stacked_CM function to see how
% -----
aa=figure(1)   %  graph of the data to show that it is all there!
clf;
subplot(2,1,1);
hold on; box on;
plot(JG , Ufs1, 'k.')
plot(JG , Ufs2, 'b.')
plot(JG , Ufs3, 'g.')
plot(JG , Ufs4, 'r.')
plot(JG , Ufs5, 'm.')
plot(JG , Ufs6, 'k.')
ylabel('U')
datetick
title('QUICK CHECK OF DATA')

subplot(2,1,2);
hold on; box on;
plot(JG , Vfs1, 'k.')
plot(JG , Vfs2, 'b.')
plot(JG , Vfs3, 'g.')
plot(JG , Vfs4, 'r.')
plot(JG , Vfs5, 'm.')
plot(JG , Vfs6, 'k.')
ylabel('V')
datetick
title('QUICK CHECK OF DATA')

set(aa,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
if ~exist([grdatdir 'otherfigure' filesep],'dir');mkdir([grdatdir 'otherfigure' filesep]);end
print('-dpng',[grdatdir 'otherfigure' filesep MOOR '_merged_beforegrid_check'])

% all the matrices for the deployments stacked together
Ufs     = [Ufs1;Ufs2;Ufs3;Ufs4;Ufs5;Ufs6];
Vfs     = [Vfs1;Vfs2;Vfs3;Vfs4;Vfs5;Vfs6];
Wfs     = [Wfs1;Wfs2;Wfs3;Wfs4;Wfs5;Wfs6];
Pfs     = [Pfs1;Pfs2;Pfs3;Pfs4;Pfs5;Pfs6];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KB outputs stagged files
z_stacked = [1000;1350;1765];

%%% add NaN where instrument is missing
fill_values = NaN(1,length(Ufs5(1,:)));
u_stacked1 = [Ufs1(1,:);fill_values;Ufs1(2,:)];
v_stacked1 = [Vfs1(1,:);fill_values;Vfs1(2,:)];
w_stacked1 = [Wfs1(1,:);fill_values;Wfs1(2,:)];
p_stacked1 = [Pfs1(1,:);fill_values;Pfs1(2,:)];

Ustacked = stack_vars(u_stacked1,Ufs2,Ufs3,Ufs4,Ufs5,Ufs6);
Vstacked = stack_vars(v_stacked1,Vfs2,Vfs3,Vfs4,Vfs5,Vfs6);
Wstacked = stack_vars(w_stacked1,Wfs2,Wfs3,Wfs4,Wfs5,Wfs6);
Pstacked = stack_vars(p_stacked1,Pfs2,Pfs3,Pfs4,Pfs5,Pfs6);

% plot(JG,Pstacked')
% datetick
save([grdatdir outputfile_stacked '.mat'],...
    'Ustacked','Vstacked','Wstacked','Pstacked','JG','z_stacked')
%%%%%%%%%%%%%%%%%%

clear Pfs1 Pfs2 Pfs3 Pfs4 Pfs5 Pfs6
clear Ufs1 Ufs2 Ufs3 Ufs4 Ufs6 Ufs6
clear Vfs1 Vfs2 Vfs3 Vfs4 Vfs5 Vfs6
clear Wfs1 Wfs2 Wfs3 Wfs4 Wfs5 Wfs6
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  4.  GRIDDING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% This should all be automatic now
%
% KB, 28/03/2024
% Removed  legacies of hydro gridding as this is not needed for velocity
% -----


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.a Vertically interpolate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% KB, 28/03/2024
% I wrote vert interpolation as function so I only need to edit one script for all three
% moorings
%
% -----
RTWB2_merg_CM = stage3_4a_vertically_interp_data(Ufs,Vfs,Wfs,Pfs,pgg,JG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.b despike time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% KB, 28/03/2024
% Moved parameter settings to top to avoid the need of changes things
% here...
%
% -----

uuu = []; vvv = []; www =[];

for i = 1 : length(RTWB2_merg_CM.UGfs(:,1))   % loop through each depth level

    [uuu(i,:),dx,ndx] = ddspike(RTWB2_merg_CM.UGfs(i,:),...
                        [-std_win*nanstd(RTWB2_merg_CM.UGfs(i,:)),...
                        std_win*nanstd(RTWB2_merg_CM.UGfs(i,:))],...
                        stddy_tol,[nloop],'y',NaN);
    [vvv(i,:),dx,ndx] = ddspike(RTWB2_merg_CM.VGfs(i,:),...
                        [-std_win*nanstd(RTWB2_merg_CM.VGfs(i,:)),...
                        std_win*nanstd(RTWB2_merg_CM.VGfs(i,:))],...
                        stddy_tol,[nloop],'y',NaN);
    [www(i,:),dx,ndx] = ddspike(RTWB2_merg_CM.WGfs(i,:),...
                        [-std_win*nanstd(RTWB2_merg_CM.WGfs(i,:)),...
                        std_win*nanstd(RTWB2_merg_CM.WGfs(i,:))],...
                        stddy_tol,[nloop],'y',NaN);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.c Interpolate over time axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% KB, 28/03/2024
% I wrote time interpolation as function so I only need to edit one script for all three
% moorings
%
% -----
RTWB2_merg_CM = stage3_4c_time_interp_data(uuu,vvv,www,pgg,JG,idepth,RTWB2_merg_CM);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  5.  PLOTTING THE GRIDDED AND MERGED PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% This should all be automatic now
%
% KB, 28/03/2024
% Removed  legacies of hydro gridding as this is not needed for velocity
% -----

contourlimituv = [-40:10:40];
contourlimitw = [-5:1:5];    
close all
figure(1); clf
subplot(3,1,1)
contourf(JG , pgg, RTWB2_merg_CM.UGfs, contourlimituv); axis ij
hold on
plot(JG,Pstacked','k')
polarmap,colorbar
caxis([min(contourlimituv) max(contourlimituv)]);
datetick; ylabel('U')
title('BEFORE DESPIKING AND INTERPOLATION')
subplot(3,1,2)
contourf(JG , pgg, uuu, contourlimituv); axis ij
hold on
plot(JG,Pstacked','k')
polarmap,colorbar
caxis([min(contourlimituv) max(contourlimituv)]);    
datetick; ylabel('U');
title('AFTER DESPIKING')
subplot(3,1,3)
contourf(JG , pgg, RTWB2_merg_CM.UGfs2, contourlimituv); axis ij
hold on
plot(JG,Pstacked','k')
polarmap,colorbar
caxis([min(contourlimituv) max(contourlimituv)]);    
datetick; ylabel('U')
title('AFTER DESPIKING AND INTERPOLATION')

set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',[grdatdir outputfile '_uuu'])

figure(2);clf
subplot(3,1,1)
contourf(JG , pgg, RTWB2_merg_CM.VGfs, contourlimituv); axis ij
hold on
plot(JG,Pstacked','k')
polarmap,colorbar
datetick; ylabel('V')
caxis([min(contourlimituv) max(contourlimituv)]);    
title('BEFORE DESPIKING AND INTERPOLATION')
subplot(3,1,2)
contourf(JG , pgg, vvv, contourlimituv); axis ij
hold on
plot(JG,Pstacked','k')
polarmap,colorbar
caxis([min(contourlimituv) max(contourlimituv)]);    
datetick; ylabel('V')
title('AFTER DESPIKING')
subplot(3,1,3)
contourf(JG , pgg, RTWB2_merg_CM.VGfs2, contourlimituv); axis ij
hold on
plot(JG,Pstacked','k')
polarmap,colorbar
caxis([min(contourlimituv) max(contourlimituv)]);    
datetick; ylabel('V')
title('AFTER DESPIKING AND INTERPOLATION')

set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
print('-dpng',[grdatdir outputfile '_vvv'])

figure(3);clf
subplot(3,1,1)
contourf(JG , pgg, RTWB2_merg_CM.WGfs, contourlimitw); axis ij
hold on
plot(JG,Pstacked','k')
polarmap,colorbar
caxis([min(contourlimitw) max(contourlimitw)]);    
datetick; ylabel('W')
title('BEFORE DESPIKING AND INTERPOLATION')
subplot(3,1,2)
contourf(JG , pgg, www, contourlimitw); axis ij
hold on
plot(JG,Pstacked','k')
polarmap,colorbar
caxis([min(contourlimitw) max(contourlimitw)]);    
datetick; ylabel('W')
title('AFTER DESPIKING')
subplot(3,1,3)
contourf(JG , pgg, RTWB2_merg_CM.WGfs2, contourlimitw); axis ij
hold on
plot(JG,Pstacked','k')
polarmap,colorbar
caxis([min(contourlimitw) max(contourlimitw)]);    
datetick; ylabel('W')
title('AFTER DESPIKING AND INTERPOLATION')

set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
print('-dpng',[grdatdir outputfile '_www'])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  6.  SAVE DATA STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update data structure with interpolated variables and save
% Notes:
% This should all be automatic now
%
% KB, 28/03/2024
% Removed  legacies of hydro gridding as this is not needed for velocity
% -----

disp(['Save data in' grdatdir outputfile '.mat'])

RTWB2_merg_CM
   
save([grdatdir outputfile],'RTWB2_merg_CM');   
 
%% functions
function Ustacked = stack_vars(U1,U2,U3,U4,U5,U6)
    Ustacked = cat(3,U1,U2,U3,U4,U5,U6);
    Ustacked = sum(Ustacked,3,'omitnan');
    Ustacked(Ustacked==0)=NaN;
end

