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
%                           pgg       22x1
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
% Jan 2017 - Adapted for OSNAP moorings (by L. Houpert)

close all
moor1 = 'CM_rtwb1_osnap_01_2014';
moor2 = 'CM_rtwb1_osnap_02_2015';
moor3 = 'CM_rtwb1_osnap_03_2016';
moor4 = 'CM_rtwb1_osnap_04_2017';
moor5 = 'CM_rtwb1_osnap_05_2018';

basedir      = pathosnap
%basedir      = '/home/sa02lh/Data/Dropbox/Work/Postdoc_OSNAP/OSNAP_mooring/backup_mdrive';
hydrodir    = [basedir '/data/moor/proc/velocity_grid/'];
grdatdir    = [basedir '/data/moor/proc/velocity_grid_merged/'];
boundarydir = [execdir 'gitrepo/stage3/gridding/CM/'];

col             = {'r','b','m','c','g','y','r--','b--','m--','c--','g--','y--','r:',...
    'b:','m:','c:','g:','y:','r.','b.','m.','c.','g.','y.','r'};


%        clear all; close all; clc;
% warning off
gridding        = 1  ;  % 1: linear, 2: using climatological profiles
bathy           = false ;  % turns on/off the bathy charts. off = flase
cm_check_plot   = true; %false ;  % turns on/off the microcat check plots. off =flase

jg_start        = datenum(2014,6,01,00,00,00);
jg_end          = datenum(2020,10,20,00,00,00);
lastyeardata    = '_2020'; % for datafilename

JG              = jg_start: 0.5: jg_end; % full time series using 2 samples per day
pgg = 0:20:2000; % depths in 20dbar bins
depthminforhoriz_interp = 40; % in case no data are available at a specific time 
% (e.g.: mooring turn around, knock-down of the mooring head) don't interpolate on a time basis for level above 40 m 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  2a.  OSNAP 1 (PE399 --> DY053)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------

disp('---------  OSNAP 1 (KN221 --> PE399) ---------')
fileID1 = fopen([boundarydir, moor1, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
% % Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
% % Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
file1_data = textscan(fileID1, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
% % Close the text file.
fclose(fileID1);
mooring1 = file1_data{1};
sn1      = file1_data{2};
cm1      = file1_data{3};
z1       = file1_data{4};
lon1     = file1_data{5};

U1 = zeros(length(sn1), length(JG));
V1 = zeros(length(sn1), length(JG));
W1 = zeros(length(sn1), length(JG));
P1 = zeros(length(sn1), length(JG));

j=1;
for i = 1: length(sn1)
    
    infile = [hydrodir, mooring1{i,:}, '_velocity_grid.mat'];
    load(infile,'dnumi','ufi','vfi','wfi','pfi');
    jdnew = dnumi;    
    sampling_rate = round(1./median(diff(jdnew))); %nominal sampling rate [per day]
    uuu    = interp1(jdnew, ufi(cm1(i),:), JG);
    vvv    = interp1(jdnew, vfi(cm1(i),:), JG);
    www    = interp1(jdnew, wfi(cm1(i),:), JG);    
    ppp    = interp1(jdnew, pfi(cm1(i),:), JG);
    Ufs1(j,:) = uuu;
    Vfs1(j,:) = vvv;
    Wfs1(j,:) = www;    
    Pfs1(j,:) = ppp;
    
    j = j + 1;
    
    if cm_check_plot
        figure(20011)
        hold on; box on;
        plot(JG, Ufs1(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 1 - U')
        
        figure(20012)
        hold on; box on;
        plot(JG, Vfs1(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 1 - V')
        
        figure(20013)
        hold on; box on;
        plot(JG, Wfs1(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 1 - W')
        
        figure(20014)
        hold on; box on;
        plot(JG, Pfs1(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 1 - PRESSURE')
        
    end
    
end

% % If a merge product of RTWB1 is available for this time period: 
% jd1 = jdnew; SGfs1 = SGfs; TGfs1 = TGfs; PG1 = p_grid;    

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  2b.  OSNAP 2 (PE399 --> DY053)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------
%
disp('---------  OSNAP 2 (PE399 --> DY053) ---------')
fileID2 = fopen([boundarydir, moor2, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
% % Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
% % Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
file2_data = textscan(fileID2, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
% % Close the text file.
fclose(fileID2);
mooring2 = file2_data{1};
sn2      = file2_data{2};
cm2      = file2_data{3};
z2       = file2_data{4};
lon2     = file2_data{5};

U2 = zeros(length(sn2), length(JG));
V2 = zeros(length(sn2), length(JG));
W2 = zeros(length(sn2), length(JG));
P2 = zeros(length(sn2), length(JG));

i = 1; j = 1;
for i = 1: length(sn2)
    
    infile = [hydrodir, mooring2{i,:}, '_velocity_grid.mat'];
    load(infile,'dnumi','ufi','vfi','wfi','pfi');
    jdnew = dnumi;    
    sampling_rate = round(1./median(diff(jdnew))); %nominal sampling rate [per day]
    uuu    = interp1(jdnew(1:end-1), ufi(cm2(i),1:end-1)', JG)';
    vvv    = interp1(jdnew(1:end-1), vfi(cm2(i),1:end-1)', JG)';
    www    = interp1(jdnew(1:end-1), wfi(cm2(i),1:end-1)', JG)';    
    ppp    = interp1(jdnew(1:end-1), pfi(cm2(i),1:end-1)', JG)';
    Ufs2(j,:) = uuu;
    Vfs2(j,:) = vvv;
    Wfs2(j,:) = www;    
    Pfs2(j,:) = ppp;
    
    j = j + 1;
    
    if cm_check_plot
        figure(20021)
        hold on; box on;
        plot(JG, Ufs2(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 2 - U')
        
        figure(20022)
        hold on; box on;
        plot(JG, Vfs2(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 2 - V')

        figure(20023)
        hold on; box on;
        plot(JG, Wfs2(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 2 - W')
        
        figure(20024)
        hold on; box on;
        plot(JG, Pfs2(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 2 - P')        
        
    end
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  2c.  OSNAP 3 (DY053 --> DY078)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------
%
disp('---------  OSNAP 3 (DY053 --> DY078) ---------')
fileID3 = fopen([boundarydir, moor3, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
% % Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
% % Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
file3_data = textscan(fileID3, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
% % Close the text file.
fclose(fileID3);
mooring3 = file3_data{1};
sn3      = file3_data{2};
cm3      = file3_data{3};
z3      = file3_data{4};
lon3     = file3_data{5};

U3 = zeros(length(sn3), length(JG));
V3 = zeros(length(sn3), length(JG));
W3 = zeros(length(sn3), length(JG));
P3 = zeros(length(sn3), length(JG));

i = 1; j = 1;
for i = 1: length(sn3)
    
    infile = [hydrodir, mooring3{i,:}, '_velocity_grid.mat'];
    load(infile,'dnumi','ufi','vfi','wfi','pfi');
    jdnew = dnumi;    
    sampling_rate = round(1./median(diff(jdnew))); %nominal sampling rate [per day]
    uuu    = interp1(jdnew(1:end-1), ufi(cm3(i),1:end-1)', JG)';
    vvv    = interp1(jdnew(1:end-1), vfi(cm3(i),1:end-1)', JG)';
    www    = interp1(jdnew(1:end-1), wfi(cm3(i),1:end-1)', JG)';    
    ppp    = interp1(jdnew(1:end-1), pfi(cm3(i),1:end-1)', JG)';
    Ufs3(j,:) = uuu;
    Vfs3(j,:) = vvv;
    Wfs3(j,:) = www;    
    Pfs3(j,:) = ppp;
    
    j = j + 1;
    
    if cm_check_plot
        figure(30021)
        hold on; box on;
        plot(JG, Ufs3(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 3 - U')
        
        figure(30022)
        hold on; box on;
        plot(JG, Vfs3(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 3 - V')

        figure(30023)
        hold on; box on;
        plot(JG, Wfs3(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 3 - W')
        
        figure(30024)
        hold on; box on;
        plot(JG, Pfs3(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 3 - P')        
        
    end
    
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  2d.  OSNAP 4 (DY078 --> AR30)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------
%
disp('---------  OSNAP 4 (DY078 --> AR30) ---------')
fileID4 = fopen([boundarydir, moor4, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
% % Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
% % Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
file4_data = textscan(fileID4, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
% % Close the text file.
fclose(fileID4);
mooring4 = file4_data{1};
sn4      = file4_data{2};
cm4      = file4_data{3};
z4      = file4_data{4};
lon4     = file4_data{5};

U4 = zeros(length(sn4), length(JG));
V4 = zeros(length(sn4), length(JG));
W4 = zeros(length(sn4), length(JG));
P4 = zeros(length(sn4), length(JG));

i = 1; j = 1;
for i = 1: length(sn4)
    
    infile = [hydrodir, mooring4{i,:}, '_velocity_grid.mat'];
    load(infile,'dnumi','ufi','vfi','wfi','pfi');
    jdnew = dnumi;    
    sampling_rate = round(1./median(diff(jdnew))); %nominal sampling rate [per day]
    uuu    = interp1(jdnew(1:end-1), ufi(cm4(i),1:end-1)', JG)';
    vvv    = interp1(jdnew(1:end-1), vfi(cm4(i),1:end-1)', JG)';
    www    = interp1(jdnew(1:end-1), wfi(cm4(i),1:end-1)', JG)';    
    ppp    = interp1(jdnew(1:end-1), pfi(cm4(i),1:end-1)', JG)';
    Ufs4(j,:) = uuu;
    Vfs4(j,:) = vvv;
    Wfs4(j,:) = www;    
    Pfs4(j,:) = ppp;
    
    j = j + 1;
    
    if cm_check_plot
        figure(40021)
        hold on; box on;
        plot(JG, Ufs4(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 4 - U')
        
        figure(40022)
        hold on; box on;
        plot(JG, Vfs4(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 4 - V')

        figure(40023)
        hold on; box on;
        plot(JG, Wfs4(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 4 - W')
        
        figure(40024)
        hold on; box on;
        plot(JG, Pfs4(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 4 - P')        
        
    end
    
end





% % If a merge product of RTWB1 is available for this time period: 
% jd2 = jdnew; SGfs2 = SGfs; TGfs2 = TGfs; PG2 = p_grid; % only keep the grid for the last microcat


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  2d.  OSNAP 5 (AR30 --> DY120)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------
%
disp('---------  OSNAP 5 (AR30 --> DY120) ---------')
fileID5 = fopen([boundarydir, moor5, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
% % Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
% % Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
file5_data = textscan(fileID5, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
% % Close the text file.
fclose(fileID5);
mooring5 = file5_data{1};
sn5      = file5_data{2};
cm5      = file5_data{3};
z5      = file5_data{4};
lon5     = file5_data{5};

U5 = zeros(length(sn5), length(JG));
V5 = zeros(length(sn5), length(JG));
W5 = zeros(length(sn5), length(JG));
P5 = zeros(length(sn5), length(JG));

i = 1; j = 1;
for i = 1: length(sn5)
    
    infile = [hydrodir, mooring5{i,:}, '_velocity_grid.mat'];
    load(infile,'dnumi','ufi','vfi','wfi','pfi');
    jdnew = dnumi;    
    sampling_rate = round(1./median(diff(jdnew))); %nominal sampling rate [per day]
    uuu    = interp1(jdnew(1:end-1), ufi(cm5(i),1:end-1)', JG)';
    vvv    = interp1(jdnew(1:end-1), vfi(cm5(i),1:end-1)', JG)';
    www    = interp1(jdnew(1:end-1), wfi(cm5(i),1:end-1)', JG)';    
    ppp    = interp1(jdnew(1:end-1), pfi(cm5(i),1:end-1)', JG)';
    Ufs5(j,:) = uuu;
    Vfs5(j,:) = vvv;
    Wfs5(j,:) = www;    
    Pfs5(j,:) = ppp;
    
    j = j + 1;
    
    if cm_check_plot
        figure(50021)
        hold on; box on;
        plot(JG, Ufs5(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 5 - U')
        
        figure(50022)
        hold on; box on;
        plot(JG, Vfs5(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 5 - V')

        figure(50023)
        hold on; box on;
        plot(JG, Wfs5(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 5 - W')
        
        figure(50025)
        hold on; box on;
        plot(JG, Pfs5(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 5 - P')        
        
    end
    
end





% % If a merge product of RTWB1 is available for this time period: 
% jd2 = jdnew; SGfs2 = SGfs; TGfs2 = TGfs; PG2 = p_grid; % only keep the grid for the last microcat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  3.  CONCATENATE AND ORDER THE MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  This step adds alll the matrices for the deployments together to form
%  large data sets.  These are then sorted into pressure order at every
%  time step before the gridding takes place.

%if cm_check_plot
aa=figure(1)   %  graph of the data to show that it is all there!
clf;
subplot(2,1,1);
hold on; box on;
plot(JG , Ufs1, 'k.')
plot(JG , Ufs2, 'b.')
plot(JG , Ufs3, 'g.')
plot(JG , Ufs4, 'r.')
plot(JG , Ufs5, 'm.')
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
ylabel('V')
datetick
title('QUICK CHECK OF DATA')

set(aa,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
print('-dpng',[grdatdir 'otherfigure' filesep  'RTWB1merged_beforegrid_check'])


% all the matrices for the deployments stacked together
Ufs     = [Ufs1;Ufs2;Ufs3;Ufs4;Ufs5];
Vfs     = [Vfs1;Vfs2;Vfs3;Vfs4;Vfs5];
Wfs     = [Wfs1;Wfs2;Wfs3;Wfs4;Wfs5];
Pfs     = [Pfs1;Pfs2;Pfs3;Pfs4;Pfs5];

% order the matrices at every time step to avoid too many NaNs creeping in
% 2004 removed....
P_sort = NaN .* ones(size(Pfs)); U_sort = NaN .* ones(size(Ufs)); V_sort = NaN .* ones(size(Vfs)); W_sort = NaN .* ones(size(Wfs)); 
j = 1;
for ii = 1: length(JG)
    [P_variable, ix] = sort(Pfs(:, ii));
    P_sort(:,j) = Pfs(ix,ii);
    U_sort(:,j) = Ufs(ix,ii);
    V_sort(:,j) = Vfs(ix,ii);
    W_sort(:,j) = Wfs(ix,ii);    
    j = j + 1;
end

% removing unused rows of the sorted matrices
Pfss = nan(size(Pfs));
Ufss = nan(size(Pfs));
Vfss = nan(size(Pfs));
Wfss = nan(size(Pfs));
i = 1; j = 1;
for i = 1: length(P_sort(:,1))
    ix = find(isnan(P_sort(i,:)));
    if length(ix) < length(JG)
        Pfss(j,:) = P_sort(i, :);
        Ufss(j,:) = U_sort(i, :);
        Vfss(j,:) = V_sort(i, :);
        Wfss(j,:) = W_sort(i, :);        
        j = j + 1;
    end
end


clear Pfs1 Pfs2 Pfs3 Pfs4 Pfs5
clear Ufs1 Ufs2 Ufs3 Ufs4 Ufs5
clear Vfs1 Vfs2 Vfs3 Vfs4 Vfs5
clear Wfs1 Wfs2 Wfs3 Wfs4 Wfs5


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


    pmin = ceil(min(Pfss)/20) * 20;
    pmax = floor(max(Pfss)/20) * 20;
    
if gridding == 2 % using climatological profiles
    outputfile = ['RTWB1_merg_' TS_CLIMATOLOGY_TP ' ' TS_CLIMATOLOGY lastyeardata ];
    if strcmp(TS_CLIMATOLOGY_TP,'annual')
        display('con_tprof0_annual')
        [TGfs, SGfs] = con_tprof0_annual(Tfss, Sfss, Pfss, pgg', 0*ones(1, size(Tfss, 2)), int_step, TSclim, ...
                                          TS_CLIMATOLOGY,TS_CLIMATOLOGY_NAME);
    elseif strcmp(TS_CLIMATOLOGY_TP,'seasonal')
        display('con_tprof0_monthly')
        [TGfs, SGfs] = con_tprof0_monthly(Tfss, Sfss, Pfss, pgg', GTV(:,2), int_step, TSclim,TS_CLIMATOLOGY, ...
                                          TS_CLIMATOLOGY_NAME );
    end
    
    % time is passed in as 0 in the above line?!?
    % pass the month of measurement into con_tprof0, instead of the empty time vector.
    %
    %	Tfss -- timeseries of discrete temperature measurements from the moorings (e.g. 25x5001)
    %       Tfss -- timeseries of discrete temperature measurements from the moorings.
    %       Tfss -- timeseries of discrete temperature measurements from the moorings.
    %       pgg   -- pressure grid onto which the profiles will be integrated
    %       GTV  -- months when each measurement was made.
    %	int_step -- integration step size [dbar] between grid points
    %	TSclim   = path to the file containing the dT/dP and dS/dP climatology
        
   
elseif gridding == 1 % linear interpolation   
    
        outputfile = ['RTWB1_merg_linear_interp' lastyeardata ];
 
        maxdepth = max(nanmean(Pfss,2));
        inodata = find(pgg>maxdepth);
        
        UGfs = nan(length(pgg),length(JG));
        VGfs = nan(length(pgg),length(JG));
        WGfs = nan(length(pgg),length(JG)); 
        UGfs_akima = nan(length(pgg),length(JG));
        VGfs_akima = nan(length(pgg),length(JG));
        WGfs_akima = nan(length(pgg),length(JG));   
        UGfs_pchip = nan(length(pgg),length(JG));
        VGfs_pchip = nan(length(pgg),length(JG));
        WGfs_pchip = nan(length(pgg),length(JG));          
       for ijj=1:length(JG)
            iuok = find(~isnan(Ufss(:,ijj)));
            ivok = find(~isnan(Vfss(:,ijj)));
            iwok = find(~isnan(Wfss(:,ijj)));
            if length(iuok)>1
                UGfs(:,ijj) = interp1(Pfss(iuok,ijj),Ufss(iuok,ijj),pgg,'linear','extrap') ;
                UGfs(inodata,ijj) = nan;
                UGfs_pchip(:,ijj) = interp1(Pfss(iuok,ijj),Ufss(iuok,ijj),pgg,'pchip','extrap') ;
                UGfs_pchip(inodata,ijj) = nan;
                if length(iuok)>2
                    UGfs_akima(:,ijj) = akima(Pfss(iuok,ijj),Ufss(iuok,ijj),pgg) ;
                    UGfs_akima(inodata,ijj) = nan;
                end
            end
            if length(ivok)>1
                VGfs(:,ijj) = interp1(Pfss(ivok,ijj),Vfss(ivok,ijj),pgg,'linear','extrap') ;    
                VGfs(inodata,ijj) = nan;        
                VGfs_pchip(:,ijj) = interp1(Pfss(ivok,ijj),Vfss(ivok,ijj),pgg,'pchip','extrap') ;  
                VGfs_pchip(inodata,ijj) = nan;
                if length(ivok)>2
                    VGfs_akima(:,ijj) = akima(Pfss(iuok,ijj),Vfss(iuok,ijj),pgg) ;
                    VGfs_akima(inodata,ijj) = nan;
                end
            end
            if length(iwok)>1
                WGfs(:,ijj) = interp1(Pfss(iwok,ijj),Wfss(iwok,ijj),pgg,'linear','extrap') ;    
                WGfs(inodata,ijj) = nan;
                WGfs_pchip(:,ijj) = interp1(Pfss(iwok,ijj),Wfss(iwok,ijj),pgg,'pchip','extrap') ;   
                WGfs_pchip(inodata,ijj) = nan;
                if length(iwok)>2
                    WGfs_akima(:,ijj) = akima(Pfss(iuok,ijj),Wfss(iuok,ijj),pgg) ;
                    WGfs_akima(inodata,ijj) = nan;
                end
            end  
       end
end


%%%%%%%%%%%%%%

    
    display(['saving: ' outputfile '.mat'])
    % Alloacte variables into a structure
    RTWB1_merg_CM.JG     = JG;
    RTWB1_merg_CM.Ufs    = Ufs;
    RTWB1_merg_CM.Vfs    = Vfs;
    RTWB1_merg_CM.Wfs    = Wfs;
    RTWB1_merg_CM.Pfs    = Pfs;    
    RTWB1_merg_CM.PGfs   = pgg';    
    RTWB1_merg_CM.UGfs   = UGfs;
    RTWB1_merg_CM.VGfs   = VGfs;
    RTWB1_merg_CM.WGfs   = WGfs;
    RTWB1_merg_CM.UGfs_akima   = UGfs_akima;
    RTWB1_merg_CM.VGfs_akima   = VGfs_akima;
    RTWB1_merg_CM.WGfs_akima   = WGfs_akima;
    RTWB1_merg_CM.UGfs_pchip   = UGfs_pchip;
    RTWB1_merg_CM.VGfs_pchip   = VGfs_pchip;
    RTWB1_merg_CM.WGfs_pchip   = WGfs_pchip;    
if gridding == 2 % using climatological profiles    
    RTWB1_merg.TS_CLIMATOLOGY   = TS_CLIMATOLOGY;
    RTWB1_merg.TS_CLIMATOLOGY_TP= TS_CLIMATOLOGY_TP;
    RTWB1_merg.TS_CLIMATOLOGY_NAME= TS_CLIMATOLOGY_NAME;
    RTWB1_merg.clim_file        = TSclim;
end    
%     RTWB1_merg.OUT_FILE         = OUT_FILE;
%     RTWB1_merg.EB_creation_date = datestr(now);
%     RTWB1_merg.MERG_REVISION    = REV(6:8);            % store the revision number of this script
%     RTWB1_merg.MERG_AUTHOR      = REV_AUTHOR(9:end-1); % store the revision author of this script
%     RTWB1_merg.MERG_DATE        = REV_DATE(7:end-1);   % store the revision date of this script
%     RTWB1_merg.EB_path = grdatdir; % path to output file
%     RTWB1_merg.function_name = function_name; % the name and path to this function
%     RTWB1_merg.interpolation_depth=idepth-20;
%     RTWB1_merg.matlab_version = version;
%     
% % bim September 2014
% % load in climatology to save in the data file
% % read in and save as a structure
% % load in climatology
% 	eval([' load ' TSclim])
% % convert to structure
% 
% 	eval(['matObj = matfile( '''  TSclim ''')' ])
%         info=whos(matObj);
% 	for kk=1:length(info);
% 	  eval([ 'clim.' info(kk).name ' = ' info(kk).name])
% 	end
%     RTWB1_merg.climatology = clim; % structure of the climatology variables used to grid the data
% 
%     RTWB1_merg
		
    
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
    %Tfs -- original stacked Temperature data from the deployments (144x5001)
    %Sfs -- original stacked salinity data from the deployments (144x5001)
    %Pfs -- original stacked Pressure data from the deployments (144x5001)
    %PGfs -- pressure grid    
    %TGfs -- temperature interpolated onto the pressure grid (PG)
    %SGfs -- salinity interpolated onto the pressure grid (PG)
    %P_sort -- Pfs sorted on pressure
    %T_sort -- Tfs sorted on pressure
    %S_sort -- Sfs sorted on pressure
    
    stddy_tol  = 10; % 4
    std_win    = 3.5; % 3.5 * std of the time series
    [nloop]    = 5; % 5
    graphics   = 'y';
    uuu = []; vvv = []; www =[];
    
    %%%%%%
    %  despike time series
    %%%%%%
    for i = 1 : length(UGfs(:,1))   % loop through each depth level
        
        [uuu(i,:),dx,ndx] = ddspike(UGfs(i,:),[-std_win*nanstd(UGfs(i,:)),...
            std_win*nanstd(UGfs(i,:))],stddy_tol,[nloop],'y',NaN);
        [vvv(i,:),dx,ndx] = ddspike(VGfs(i,:),[-std_win*nanstd(VGfs(i,:)),...
            std_win*nanstd(VGfs(i,:))],stddy_tol,[nloop],'y',NaN);
        [www(i,:),dx,ndx] = ddspike(WGfs(i,:),[-std_win*nanstd(WGfs(i,:)),...
            std_win*nanstd(WGfs(i,:))],stddy_tol,[nloop],'y',NaN);
    end
    
    [m,n] = size(UGfs);
    UG_2 = NaN * ones(m,n); VG_2 = NaN * ones(m,n); WG_2 = NaN * ones(m,n);
    
    % GDM, 9/4/2013
    % Need to interpolate horizontally over nans but...
    % Don't want to create fake values above knocked down moorings
    % So select a depth below which, interpoation happens
    
    idepth = depthminforhoriz_interp; % multiples of 20
    I = find(pgg == idepth);
    
    % copy the top idepth into the new file
    % no temporal interpolation
    UG_2([1:I-1],:)=UGfs([1:I-1],:);
    VG_2([1:I-1],:)=VGfs([1:I-1],:);
    WG_2([1:I-1],:)=WGfs([1:I-1],:);
    
    %    i = 1; j = 1;
    %    for i = 1: 6 % top 100m based on a grid of 20dbar
    %        SG_2(j,:) = SGfs(i, :);
    %        TG_2(j,:) = TGfs(i, :);
    %        j = j + 1;
    %    end
%
    
    % this piece of code has the ability to create spurious values of T & S ** tk
    % how large are the temporal gaps in the time series ?
    % linear should only be used for small gaps.
    
    %     i = 7; j = 7; % PG(7) = 120;
    %     i = 11; j = 11; % PG(11) = 120;

    i = I; j = I;
    for i = I: length(pgg) % for each depth
        
        % locate all non nan values in the despiked zonal velocity
        iu = find(~isnan(uuu(i,:)));
        if length(iu) < 2
            continue
        end
        % locate all non nan values in the despiked meridional velocity
        iv = find(~isnan(vvv(i,:)));
        % locate all non nan values in the despiked vertical velocity
        iw = find(~isnan(www(i,:)));
          
        % interpolate in time over the missing data
        VG_2(j,:) = interp1(JG(iv), vvv(i, iv), JG);
        % interpolate in time over the missing data
        UG_2(j,:) = interp1(JG(iu), uuu(i, iu), JG);
        % interpolate in time over the missing data
        WG_2(j,:) = interp1(JG(iw), www(i, iw), JG);
        
        
        % set NAN for nans gap of more than 10 days             
        ibad = find([diff(iv)]>20); 
        if length(ibad)>=1
            for isss=1:length(ibad)  
                VG_2(i,iv(ibad(isss)):iv(ibad(isss)+1))=nan;
                UG_2(i,iv(ibad(isss)):iv(ibad(isss)+1))=nan;
                WG_2(i,iv(ibad(isss)):iv(ibad(isss)+1))=nan;
            end
        end
        
        j = j + 1;
    end
    
    
   contourlimituv = [-40:10:40];
    contourlimitw = [-5:1:5];    
    close all
    figure(1); clf
    subplot(3,1,1)
    contourf(JG , pgg, UGfs, contourlimituv); axis ij
    caxis([min(contourlimituv) max(contourlimituv)]);
    datetick; ylabel('U')
    title('BEFORE DESPIKING AND INTERPOLATION')
    subplot(3,1,2)
    contourf(JG , pgg, uuu, contourlimituv); axis ij
    caxis([min(contourlimituv) max(contourlimituv)]);    
    datetick; ylabel('U');
    title('AFTER DESPIKING')
    subplot(3,1,3)
    contourf(JG , pgg, UG_2, contourlimituv); axis ij
    caxis([min(contourlimituv) max(contourlimituv)]);    
    datetick; ylabel('U')
    title('AFTER DESPIKING AND INTERPOLATION')
    
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
    print('-dpng',[grdatdir outputfile '_uuu'])
    
    figure(2);clf
    subplot(3,1,1)
    contourf(JG , pgg, VGfs, contourlimituv); axis ij
    datetick; ylabel('V')
    caxis([min(contourlimituv) max(contourlimituv)]);    
    title('BEFORE DESPIKING AND INTERPOLATION')
    subplot(3,1,2)
    contourf(JG , pgg, vvv, contourlimituv); axis ij
    caxis([min(contourlimituv) max(contourlimituv)]);    
    datetick; ylabel('V')
    title('AFTER DESPIKING')
    subplot(3,1,3)
    contourf(JG , pgg, VG_2, contourlimituv); axis ij
    caxis([min(contourlimituv) max(contourlimituv)]);    
    datetick; ylabel('V')
    title('AFTER DESPIKING AND INTERPOLATION')
    
     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
    print('-dpng',[grdatdir outputfile '_vvv'])

    figure(2);clf
    subplot(3,1,1)
    contourf(JG , pgg, WGfs, contourlimitw); axis ij
    caxis([min(contourlimitw) max(contourlimitw)]);    
    datetick; ylabel('W')
    title('BEFORE DESPIKING AND INTERPOLATION')
    subplot(3,1,2)
    contourf(JG , pgg, www, contourlimitw); axis ij
    caxis([min(contourlimitw) max(contourlimitw)]);    
    datetick; ylabel('W')
    title('AFTER DESPIKING')
    subplot(3,1,3)
    contourf(JG , pgg, WG_2, contourlimitw); axis ij
    caxis([min(contourlimitw) max(contourlimitw)]);    
    datetick; ylabel('W')
    title('AFTER DESPIKING AND INTERPOLATION')
    
     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
    print('-dpng',[grdatdir outputfile '_www'])
    
    % Alloacte variables into a structure
    RTWB1_merg_CM.UGfs2= UG_2;
    RTWB1_merg_CM.VGfs2= VG_2;
    RTWB1_merg_CM.WGfs2= WG_2;
        
    %    ['save ' grdatdir char(OUT_FILE(3)) ...
    %           ' JG pgg Tfs Sfs Pfs TGfs SGfs TG_2 SG_2 P_sort T_sort S_sort']
    %    eval(['save ' grdatdir char(OUT_FILE(3)) ...
    %           ' JG pgg Tfs Sfs Pfs TGfs SGfs TG_2 SG_2 P_sort T_sort S_sort'])
    
    RTWB1_merg_CM.comment{1,1}= 'JG -- julian day';
    RTWB1_merg_CM.comment{2,1}= 'Ufs -- original stacked zonal velocity data from the deployments';
    RTWB1_merg_CM.comment{3,1}= 'Vfs -- original stacked meridional data from the deployments ';
    RTWB1_merg_CM.comment{4,1}= 'Wfs -- original stacked vertical data from the deployments ';    
    RTWB1_merg_CM.comment{5,1}= 'Pfs -- original stacked Pressure data from the deployments';    
    RTWB1_merg_CM.comment{6,1}= 'PGfs -- pressure grid ';
    RTWB1_merg_CM.comment{7,1}= 'UGfs -- zonal velocity interpolated onto the pressure grid (PGfs) with a linear extrapolation for surface values';
    RTWB1_merg_CM.comment{8,1}= 'VGfs -- meridional velocity interpolated onto the pressure grid (PGfs) with a linear extrapolation for surface values';
    RTWB1_merg_CM.comment{9,1}= 'WGfs -- vertical velocity interpolated onto the pressure grid (PGfs) with a linear extrapolation for surface values';    
    RTWB1_merg_CM.comment{10,1}= 'UGfs_akima -- zonal velocity interpolated onto the pressure grid (PGfs) with akima method';
    RTWB1_merg_CM.comment{11,1}= 'VGfs_akima -- meridional velocity interpolated onto the pressure grid (PGfs) with akima method';
    RTWB1_merg_CM.comment{12,1}= 'WGfs_akima -- vertical velocity interpolated onto the pressure grid (PGfs) with akima method';    
    RTWB1_merg_CM.comment{13,1}= 'UGfs_pchip -- zonal velocity interpolated onto the pressure grid (PGfs) with a pchip extrapolation for surface values';
    RTWB1_merg_CM.comment{14,1}= 'VGfs_pchip -- meridional velocity interpolated onto the pressure grid (PGfs) with a pchip extrapolation for surface values';
    RTWB1_merg_CM.comment{15,1}= 'WGfs_pchip -- vertical velocity interpolated onto the pressure grid (PGfs) with a pchip extrapolation for surface values';    
    RTWB1_merg_CM.comment{16,1}= 'UGfs2 -- zonal velocity interpolated onto the time grid (JG) after despiking with a linear extrapolation for surface values';  
    RTWB1_merg_CM.comment{17,1}= 'VGfs2 -- meridional velocity interpolated onto the time grid (JG) after despiking with a linear extrapolation for surface values';  
    RTWB1_merg_CM.comment{18,1}= 'WGfs2 -- vertical velocity interpolated onto the time grid (JG) after despiking with a linear extrapolation for surface values';  
    
   
    RTWB1_merg_CM

    

%     ['save ' grdatdir outputfile ' RTWB1_merg_CM']
%     eval(['save ' grdatdir outputfile ' RTWB1_merg_CM']);   
      save([grdatdir outputfile],'RTWB1_merg_CM');   

    
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%  6.  PLOTTING THE GRIDDED AND MERGED PROFILES
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% [c,h]=contourf(JG,PGfs,VGfs2,'LineColor','none');
% cmocean('balance','pivot',0);
% axis ij
% datetick('x',12, 'keeplimits')
% ylabel('Pressure (db)');
% C=colorbar;
% ylabel(C,' w velocity (cm s^{-1})')


% 
% figure(1000)
% clf
% subplot(3,1,1)
% size(JG); size(PG1); size(SGfs1); size(TGfs1); 
% contourf(jd1 , PG1, SGfs1) % doesnt work as there is no  merge RTWB mooring available
% axis ij
% hold on
% contourf(jd2 , PG2, SGfs2)
% % contourf(jd3a , PG3a, SGfs3a)
% % contourf(jd3b , PG3b, SGfs3b)
% % contourf(jd4 , PG4, SGfs4)
% % contourf(jd5 , PG5, SGfs5)
% %contourf(jd6 , PG6, SGfs6)
% %contourf(jd7 , PG7, SGfs7)
% %contourf(jd8 , PG8, SGfs8)
% ylim([0,2000])
% xlim([jd1(1) , JG(end) ])
% title('CURRENT MERGING - SALINITY')
% 
% subplot(3,1,2)
% contourf(JG , pgg, SGfs)
% axis ij
% ylim([0,2000])
% xlim([jd1(1) , JG(end) ])
% title('NEW MERGING BEFORE INTERPOLATION - SALINITY')
% 
% 
% figure(1001)
% clf
% subplot(3,1,1)
% contourf(jd1 , PG1, TGfs1)
% axis ij
% hold on
% contourf(jd2 , PG2, TGfs2)
% % contourf(jd3a , PG3a, TGfs3a)
% % contourf(jd3b , PG3b, TGfs3b)
% % contourf(jd4 , PG4, TGfs4)
% % contourf(jd5 , PG5, TGfs5)
% %contourf(jd6 , PG6, TGfs6)
% %contourf(jd7 , PG7, TGfs7)
% %contourf(jd8 , PG8, TGfs8)
% ylim([0,2000])
% xlim([jd1(1) , JG(end) ])
% title('CURRENT MERGING - TEMPERATURE')
% 
% subplot(3,1,2)
% contourf(JG , pgg, TGfs)
% axis ij
% ylim([0,2000])
% xlim([jd1(1) , JG(end) ])
% title('NEW MERGING BEFORE INTERPOLATION- TEMPERATURE')
% 
% 
% 
% if gridding == 0
%     load ../mat_files/EB_merged_data_2010.mat
% end
% 
% figure(1000)
% subplot(3,1,3)
% contourf(JG , pgg, SG_2)
% axis ij
% ylim([0,2000])
% xlim([jd1(1) , JG(end) ])
% title('NEW MERGING (after interp - SALINITY')
% 
% 
% figure(1001)
% subplot(3,1,3)
% contourf(JG , pgg, TG_2)
% axis ij
% ylim([0,2000])
% xlim([jd1(1) , JG(end) ])
% title('NEW MERGING (after interp - TEMPERATURE')
% 



