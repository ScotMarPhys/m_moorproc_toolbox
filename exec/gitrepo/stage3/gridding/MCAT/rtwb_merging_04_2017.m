%%  Experimental Development Code for the Merging of the Western Boundary
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
% Jan 2017 - Adapted for OSNAP mooring (by L. Houpert)


close all
moor1 = 'rtwb_osnap_01_2014';
moor2 = 'rtwb_osnap_02_2015';
moor3 = 'rtwb_osnap_03_2016';
moor4 = 'rtwb_osnap_04_2017';

basedir      = [pathosnap filesep];
hydrodir    = [basedir 'data/moor/proc/hydro_grid/'];
grdatdir    = [basedir 'data/moor/proc/hydro_grid_merged/'];
boundarydir = [basedir 'exec/ar30/stage3/gridding/MCAT/']; %[basedir 'users/loh/moor_merging_gridding/'];
%boundarydir = ['.' filesep]; %[basedir 'users/loh/moor_merging_gridding/'];

col             = {'r','b','m','c','g','y','r--','b--','m--','c--','g--','y--','r:',...
    'b:','m:','c:','g:','y:','r.','b.','m.','c.','g.','y.','r'};


%        clear all; close all; clc;
% warning off
gridding        = 1  ;  % 1: linear, 2: using climatological profiles
bathy           = false ;  % turns on/off the bathy charts. off = flase
mc_check_plot   = true; %false ;  % turns on/off the microcat check plots. off =false

jg_start        = datenum(2014,6,01,00,00,00);
jg_end          = datenum(2018,7,10,00,00,00);
lastyeardata    = '_2017'; % for datafilename

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
mc1      = file1_data{3};
z1       = file1_data{4};
lon1     = file1_data{5};

T1 = zeros(length(sn1), length(JG));
P1 = zeros(length(sn1), length(JG));
S1 = zeros(length(sn1), length(JG));

i = 1; j = 1;
for i = 1: length(sn1)
    
    infile = [hydrodir, mooring1{i,:}, '_grid.mat'];
    load(infile);
    jdnew = datenum(gregorian(jd));    
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jdnew(1:end-1), Sfs(mc1(i),1:end-1)', JG)';
    temp        = interp1(jdnew(1:end-1), Tfs(mc1(i),1:end-1)', JG)';
    pressure    = interp1(jdnew(1:end-1), Pfs(mc1(i),1:end-1)', JG)';
    Sfs1(j,:) = salinity;
    Pfs1(j,:) = pressure;
    Tfs1(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(20011)
        hold on; box on;
        plot(JG, Sfs1(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 2 - SALINITY')
        
        figure(20012)
        hold on; box on;
        plot(JG, Tfs1(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 2 - TEMPERATURE')
        
        figure(200120)
        hold on; box on;
        plot(JG, sw_pden(Sfs1(i, :),Tfs1(i, :),Pfs1(i, :),0)-1000, col{:,i},'LineWidth',2);
        title('OSNAP 1 - POTENTIAL DENSITY')      
        
        figure(20013)
        hold on; box on;
        plot(JG, Pfs1(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 1 - PRESSURE')
        
    end
    
end

% % If a merge product of RTEB is available for this time period: 
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
mc2      = file2_data{3};
z2       = file2_data{4};
lon2     = file2_data{5};

T2 = zeros(length(sn2), length(JG));
P2 = zeros(length(sn2), length(JG));
S2 = zeros(length(sn2), length(JG));

i = 1; j = 1;
for i = 1: length(sn2)
    
    infile = [hydrodir, mooring2{i,:}, '_grid.mat'];
    load(infile);
    jdnew = datenum(gregorian(jd));    
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jdnew(1:end-1), Sfs(mc2(i),1:end-1)', JG)';
    temp        = interp1(jdnew(1:end-1), Tfs(mc2(i),1:end-1)', JG)';
    pressure    = interp1(jdnew(1:end-1), Pfs(mc2(i),1:end-1)', JG)';
    Sfs2(j,:) = salinity;
    Pfs2(j,:) = pressure;
    Tfs2(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(20021)
        hold on; box on;
        plot(JG, Sfs2(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 2 - SALINITY')
        
        figure(20022)
        hold on; box on;
        plot(JG, Tfs2(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 2 - TEMPERATURE')
        
        figure(200220)
        hold on; box on;
        plot(JG, sw_pden(Sfs2(i, :),Tfs2(i, :),Pfs2(i, :),0)-1000, col{:,i},'LineWidth',2);
        title('OSNAP 2 - POTENTIAL DENSITY')      
    
        figure(20023)
        hold on; box on;
        plot(JG, Pfs2(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 2 - PRESSURE')
    end
    
end

% % If a merge product of RTEB is available for this time period: 
% jd2 = jdnew; SGfs2 = SGfs; TGfs2 = TGfs; PG2 = p_grid; % only keep the grid for the last microcat

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  2c.  OSNAP 3 (DY053 --> DY078)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------
%
disp('---------  OSNAP 3 (DY053--> DY078) ---------')
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
mc3      = file3_data{3};
z3       = file3_data{4};
lon3     = file3_data{5};

T3 = zeros(length(sn3), length(JG));
P3 = zeros(length(sn3), length(JG));
S3 = zeros(length(sn3), length(JG));

i = 1; j = 1;
for i = 1: length(sn3)
    
    infile = [hydrodir, mooring3{i,:}, '_grid.mat'];
    load(infile);
    jdnew = datenum(gregorian(jd));    
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jdnew(1:end-1), Sfs(mc3(i),1:end-1)', JG)';
    temp        = interp1(jdnew(1:end-1), Tfs(mc3(i),1:end-1)', JG)';
    pressure    = interp1(jdnew(1:end-1), Pfs(mc3(i),1:end-1)', JG)';
    Sfs3(j,:) = salinity;
    Pfs3(j,:) = pressure;
    Tfs3(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(30021)
        hold on; box on;
        plot(JG, Sfs3(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 3 - SALINITY')
        
        figure(30022)
        hold on; box on;
        plot(JG, Tfs3(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 3 - TEMPERATURE')
        
        figure(300220)
        hold on; box on;
        plot(JG, sw_pden(Sfs3(i, :),Tfs3(i, :),Pfs3(i, :),0)-1000, col{:,i},'LineWidth',2);
        title('OSNAP 3 - POTENTIAL DENSITY')      
    
        figure(30023)
        hold on; box on;
        plot(JG, Pfs3(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 3 - PRESSURE')
    end
    
end



% % If a merge product of RTEB is available for this time period: 
% jd2 = jdnew; SGfs2 = SGfs; TGfs2 = TGfs; PG2 = p_grid; % only keep the grid for the last microcat

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  2c.  OSNAP 4 (DY078 --> AR30)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Notes:
%  ------
%
disp('---------  OSNAP 4 (DY078--> AR30) ---------')
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
mc4      = file4_data{3};
z4       = file4_data{4};
lon4     = file4_data{5};

T4 = zeros(length(sn4), length(JG));
P4 = zeros(length(sn4), length(JG));
S4 = zeros(length(sn4), length(JG));

i = 1; j = 1;
for i = 1: length(sn4)
    
    infile = [hydrodir, mooring4{i,:}, '_grid.mat'];
    load(infile);
    jdnew = datenum(gregorian(jd));    
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jdnew(1:end-1), Sfs(mc4(i),1:end-1)', JG)';
    temp        = interp1(jdnew(1:end-1), Tfs(mc4(i),1:end-1)', JG)';
    pressure    = interp1(jdnew(1:end-1), Pfs(mc4(i),1:end-1)', JG)';
    Sfs4(j,:) = salinity;
    Pfs4(j,:) = pressure;
    Tfs4(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(40021)
        hold on; box on;
        plot(JG, Sfs4(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 4 - SALINITY')
        
        figure(40022)
        hold on; box on;
        plot(JG, Tfs4(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 4 - TEMPERATURE')
        
        figure(400220)
        hold on; box on;
        plot(JG, sw_pden(Sfs4(i, :),Tfs4(i, :),Pfs4(i, :),0)-1000, col{:,i},'LineWidth',2);
        title('OSNAP 4 - POTENTIAL DENSITY')      
    
        figure(40023)
        hold on; box on;
        plot(JG, Pfs4(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 4 - PRESSURE')
    end
    
end





close all

% % If a merge product of RTEB is available for this time period: 
% jd3 = jdnew; SGfs3 = SGfs; TGfs3 = TGfs; PG3 = p_grid; % only keep the grid for the last microcat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  3.  CONCATENATE AND ORDER THE MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  This step adds alll the matrices for the deployments together to form
%  large data sets.  These are then sorted into pressure order at every
%  time step before the gridding takes place.

%if mc_check_plot
    PDENfs1 = sw_pden(Sfs1,Tfs1,Pfs1,0) - 1000;
    PDENfs2 = sw_pden(Sfs2,Tfs2,Pfs2,0) - 1000;  
    PDENfs3 = sw_pden(Sfs3,Tfs3,Pfs3,0) - 1000;
    PDENfs4 = sw_pden(Sfs4,Tfs4,Pfs4,0) - 1000;

    
    figure(11)   %  graph of the data to show that it is all there!
    clf;
    subplot(2,1,1);
    hold on; box on;
    plot(JG , Tfs1, 'k.')
    plot(JG , Tfs2, 'b.')
    plot(JG , Tfs3, 'g.')
    plot(JG , Tfs4, 'r.')
    ylabel('C')
    datetick
    title('TEMPERATURE')

    subplot(2,1,2);
    hold on; box on;
    plot(JG , Sfs1, 'k.')
    plot(JG , Sfs2, 'b.')
    plot(JG , Sfs3, 'g.')  
    plot(JG , Sfs4, 'r.')  
%     plot(JG , Sfs3a, 'r.')
%     plot(JG , Sfs3b, 'y.')
%     plot(JG , Sfs4, 'b.')
%     plot(JG , Sfs5, 'r.')
%     plot(JG , Sfs6, 'g.')
%     plot(JG , Sfs7, 'k.')
%     plot(JG , Sfs8, 'b.')
%     plot(JG , Sfs9, 'r.')
%     plot(JG , Sfs10, 'g.')
    ylabel('SAL.')
    datetick
    title('SALINITY')
    
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
print('-dpng',[grdatdir 'otherfigure' filesep  'RTWBmerged_beforegrid_check1'])

   figure(10)   %  graph of the data to show that it is all there!
    clf;
    subplot(2,1,1);
    hold on; box on;
    plot(JG , Pfs1, 'k.')
    plot(JG , Pfs2, 'b.')
    plot(JG , Pfs3, 'g.') 
    plot(JG , Pfs4, 'r.')
    ylabel('dbar')
    datetick
    title('PRES')   
    
    subplot(2,1,2);
    hold on; box on;
    plot(JG , PDENfs1, 'k.')
    plot(JG , PDENfs2, 'b.')
    plot(JG , PDENfs3, 'g.')  
    plot(JG , PDENfs4, 'r.') 
    ylabel('kg/m3.')
    datetick
    title('POT. DENS')   


set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
print('-dpng',[grdatdir 'otherfigure' filesep  'RTWBmerged_beforegrid_check2'])

%end

% all the matrices for the deployments stacked together
Pfs     = [Pfs1;Pfs2;Pfs3;Pfs4];
Sfs     = [Sfs1;Sfs2;Sfs3;Sfs4];
Tfs     = [Tfs1;Tfs2;Tfs3;Tfs4];

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
Pfss = nan(size(Pfs));
Tfss = nan(size(Pfs));
Sfss = nan(size(Pfs));
i = 1; j = 1;
for i = 1: length(P_sort(:,1))
    ix = find(isnan(P_sort(i,:)));
    if length(ix) < length(JG)
        Pfss(j,:) = P_sort(i, :);
        Tfss(j,:) = T_sort(i, :);
        Sfss(j,:) = S_sort(i, :);
        j = j + 1;
    end
end

clear Pfs1 Pfs2 Pfs3 Pfs4
clear Sfs1 Sfs2 Sfs3 Sfs4
clear Tfs1 Tfs2 Tfs3 Tfs4
 

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
    outputfile = ['RTWB_merg_' TS_CLIMATOLOGY_TP ' ' TS_CLIMATOLOGY lastyeardata ];
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
    
        outputfile = ['RTWB_merg_linear_interp' lastyeardata ];
        TGfs = nan(length(pgg),length(JG));
        SGfs = nan(length(pgg),length(JG)); 
        for ijj=1:length(JG)
            itok = find(~isnan(Tfss(:,ijj)));
            isok = find(~isnan(Sfss(:,ijj)));
            if length(itok)>1
                TGfs(:,ijj) = interp1(Pfss(itok,ijj),Tfss(itok,ijj),pgg) ; 
            end
            if length(isok)>1
                SGfs(:,ijj) = interp1(Pfss(isok,ijj),Sfss(isok,ijj),pgg) ;     
            end
        end
       
    
end 




    %%%%%%%%%%%%%%

    
    display(['saving: ' outputfile '.mat'])
    % Alloacte variables into a structure
    RTWB_merg.JG     = JG;
    RTWB_merg.Tfs    = Tfs;
    RTWB_merg.Sfs    = Sfs;
    RTWB_merg.Pfs    = Pfs;
    RTWB_merg.PGfs    = pgg';    
    RTWB_merg.TGfs   = TGfs;
    RTWB_merg.SGfs   = SGfs;
%     RTEB_merg.P_sort = P_sort;
%     RTEB_merg.T_sort = T_sort;
%     RTEB_merg.S_sort = S_sort;
if gridding == 2 % using climatological profiles    
    RTWB_merg.TS_CLIMATOLOGY   = TS_CLIMATOLOGY;
    RTWB_merg.TS_CLIMATOLOGY_TP= TS_CLIMATOLOGY_TP;
    RTWB_merg.TS_CLIMATOLOGY_NAME= TS_CLIMATOLOGY_NAME;
    RTWB_merg.clim_file        = TSclim;
end    
%     RTEB_merg.OUT_FILE         = OUT_FILE;
%     RTEB_merg.EB_creation_date = datestr(now);
%     RTEB_merg.MERG_REVISION    = REV(6:8);            % store the revision number of this script
%     RTEB_merg.MERG_AUTHOR      = REV_AUTHOR(9:end-1); % store the revision author of this script
%     RTEB_merg.MERG_DATE        = REV_DATE(7:end-1);   % store the revision date of this script
%     RTEB_merg.EB_path = grdatdir; % path to output file
%     RTEB_merg.function_name = function_name; % the name and path to this function
%     RTEB_merg.interpolation_depth=idepth-20;
%     RTEB_merg.matlab_version = version;
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
%     RTEB_merg.climatology = clim; % structure of the climatology variables used to grid the data
% 
%     RTEB_merg
		
    
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
    
    idepth = depthminforhoriz_interp; % multiples of 20
    I = find(pgg == idepth);
    
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
    for i = I: length(pgg) % for each depth
        
        % locate all non nan values in the despiked temperature
        it = find(~isnan(temp(i,:)));
        if length(it) < 2
            continue
        end
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
    contourf(JG , pgg, TGfs,8); axis ij
    datetick; ylabel('TEMPERATURE')
    title('RTWB BEFORE DESPIKING AND INTEROLATION')
    subplot(3,1,2)
    contourf(JG , pgg, temp,8); axis ij
    datetick; ylabel('TEMPERATURE');
    title('AFTER DESPIKING')
    subplot(3,1,3)
    contourf(JG , pgg, TG_east,8); axis ij
    datetick; ylabel('TEMPERATURE')
    title('AFTER DESPIKING AND INTERPOLATION')
    
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
    print('-dpng',[grdatdir outputfile '_temperature'])
    
    figure(2);clf
    subplot(3,1,1)
    contourf(JG , pgg, SGfs,8); axis ij
    caxis([min(SGfs(:)) max(SGfs(:))])
    datetick; ylabel('SALINITY')
    title('RTWB BEFORE DESPIKING AND INTEROLATION')
    subplot(3,1,2)
    contourf(JG , pgg, salinity,8); axis ij
      caxis([min(SGfs(:)) max(SGfs(:))])
    datetick; ylabel('SALINITY')
    title('AFTER DESPIKING')
    subplot(3,1,3)
    contourf(JG , pgg, SG_east,8); axis ij
      caxis([min(SGfs(:)) max(SGfs(:))])
    datetick; ylabel('SALINITY')
    title('AFTER DESPIKING AND INTERPOLATION')
    
     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
    print('-dpng',[grdatdir outputfile '_salinity'])

    
    PDENGfs = sw_pden(SGfs,TGfs,pgg',0);
    pden = sw_pden(salinity,temp,pgg',0);   
    PDENG_east = sw_pden(SG_east,TG_east,pgg',0);   
    figure(3);clf
    subplot(3,1,1)
    contourf(JG , pgg, PDENGfs,10); axis ij
    caxis([min(PDENGfs(:)) max(PDENGfs(:))])
    datetick; ylabel('POT. DENS.')
    title('RTWB BEFORE DESPIKING AND INTEROLATION')
    subplot(3,1,2)
    contourf(JG , pgg, pden,10); axis ij
      caxis([min(PDENGfs(:)) max(PDENGfs(:))])
    datetick; ylabel('POT. DENS.')
    title('AFTER DESPIKING')
    subplot(3,1,3)
    contourf(JG , pgg,PDENG_east,10); axis ij
      caxis([min(PDENGfs(:)) max(PDENGfs(:))])
    datetick; ylabel('POT. DENS.')
    title('AFTER DESPIKING AND INTERPOLATION')
    
     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
    print('-dpng',[grdatdir outputfile '_potdens'])
    
    
    
    % Alloacte variables into a structure
    RTWB_merg.TGfs2= TG_east;
    RTWB_merg.SGfs2= SG_east;
    
    %    ['save ' grdatdir char(OUT_FILE(3)) ...
    %           ' JG pgg Tfs Sfs Pfs TGfs SGfs TG_east SG_east P_sort T_sort S_sort']
    %    eval(['save ' grdatdir char(OUT_FILE(3)) ...
    %           ' JG pgg Tfs Sfs Pfs TGfs SGfs TG_east SG_east P_sort T_sort S_sort'])
    
    
    RTWB_merg.comment{1,1}= 'JG -- julian day';
    RTWB_merg.comment{2,1}= 'Tfs -- original stacked Temperature data from the deployments';
    RTWB_merg.comment{3,1}= 'Sfs -- original stacked salinity data from the deployments ';
    RTWB_merg.comment{4,1}= 'Pfs -- original stacked Pressure data from the deployments';    
    RTWB_merg.comment{5,1}= 'PGfs -- pressure grid ';
    RTWB_merg.comment{6,1}= 'TGfs -- temperature interpolated onto the pressure grid (PGfs)';
    RTWB_merg.comment{7,1}= 'SGfs -- salinity interpolated onto the pressure grid (PGfs)';
    RTWB_merg.comment{8,1}= 'TGfs2 -- temperature interpolated onto the time grid (JG) after despiking';  
    RTWB_merg.comment{9,1}= 'SGfs2 -- salinity interpolated onto the time grid (JG) after despiking';  
   
    ['save ' grdatdir outputfile ' RTWB_merg'];
    eval(['save ' grdatdir outputfile ' RTWB_merg']);   
  
    
    RTWB_merg ; 
    
    
    %   JG -- julian day
    %   Tfs -- original stacked Temperature data from the deployments (144x5001)
    %   Sfs -- original stacked salinity data from the deployments (144x5001)
    %   Pfs -- original stacked Pressure data from the deployments (144x5001)
    %   PGfs -- pressure grid    
    %   TGfs -- temperature interpolated onto the pressure grid (PG)
    %   SGfs -- salinity interpolated onto the pressure grid (PG)
    %	TGfs2 -- temperature interpolated onto the time grid (JG) after despiking
    %	SGfs2 -- salinity interpolated onto the time grid (JG) after despiking
    %	P_sort -- Pfs sorted on pressure. no despike. no temporal interp
    %	T_sort -- Tfs sorted on pressure. no despike. no temporal interp
    %	S_sort -- Sfs sorted on pressure. no despike. no temporal interp




% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%  6.  PLOTTING THE GRIDDED AND MERGED PROFILES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
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
% contourf(JG , pgg, SG_east)
% axis ij
% ylim([0,2000])
% xlim([jd1(1) , JG(end) ])
% title('NEW MERGING (after interp - SALINITY')
% 
% 
% figure(1001)
% subplot(3,1,3)
% contourf(JG , pgg, TG_east)
% axis ij
% ylim([0,2000])
% xlim([jd1(1) , JG(end) ])
% title('NEW MERGING (after interp - TEMPERATURE')
% 



