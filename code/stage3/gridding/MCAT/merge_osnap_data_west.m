%% Code for the Merging of the Western Boundary
%% MicroCAT Data
% _*Original authors: Paul Wright, March 2010 (Updated Ben Moat, March and Nov 
% 2012)*_
% 
% The selcted P, S and T data are obtained from the mooring .microcat files 
% and interpolated onto the full timeseries time grid, sorted into pressure levels 
% and then gridded. uses meta data files to select the MicroCATs. it would be 
% possible to simplify this further but it is debatable whether it may be better 
% to leave a bit "open" for the ease of updating on an annual basis.  Clearer 
% too.
% 
% Notes specific to each deployment period are included within the loops.
% Code History
%% 
% * 19 April 2010 - removed final data point of each MicroCAT to attempt to 
% remove the spiking at the joins between the deployments.
% * 26 May 2010 - changed gridded function to con_tprof_v3.m which stops the 
% gridding at the upper MicroCAT.  Interpolates over data gaps less than a year 
% in order to match Kanzow's original code and present data for the MOC code.
% * 9 June 2010 - PGW changed the climatology TSclim
% * 28 June 2010 - changed climatology back to original.  Updated to include 
% D344 data
% * 11 August 2010 - Tidied up and improved comments slightly
% * June 2011 - Updated merging Apr 2009 - Dec 2010
% * Sept 2013 - added in name of cliamtology (e.g. hbase or argo)
% * Jan 2017 - Adapted for OSNAP mooring (by L. Houpert)
% * July 2020, further updates and tidy (Lewis Drysdale)
% * December 2020, update for DY120 (Lewis Drysdale)

close('all')
% MOORING FILENAMES
% if GRID_OSNAP_MCAT_DATA_INDEV has been used to get moor filenames, use them 
% else use list

if exist('moor','var')
else
    moor{1} = 'rtwb_osnap_01_2014';
    moor{2} = 'rtwb_osnap_02_2015';
    moor{3} = 'rtwb_osnap_03_2016';
    moor{4} = 'rtwb_osnap_04_2017';
    moor{5} = 'rtwb_osnap_05_2018';
    moor{6} = 'rtwb_osnap_06_2020';
end
% SET PLOT COLORS

col             = {'r','b','m','c','g','y','r--','b--','m--','c--','g--','y--','r:',...
    'b:','m:','c:','g:','y:','r.','b.','m.','c.','g.','y.','r'};
% CHECK NO DUPLICATES

[~,idx]=unique(cell2mat(moor),'rows');
moor =  moor(idx,:);

% DEFINE OUTPUTFILES
outputfile_stacked = ['Hydro_rtwb_stacked_' jg_start_str '_' jg_end_str,...
                        '_' data_version];

%% 1. READ MOORING DATA
% OSNAP 1 (KN221 --> PE399)

disp('---------  OSNAP 1 (KN221 --> PE399) ---------')
fileID1 = fopen([boundarydir, moor{1}, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
file1_data = textscan(fileID1, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID1);
mooring1 = file1_data{1};
sn1      = file1_data{2};
mc1      = file1_data{3};
z1       = file1_data{4};
lon1     = file1_data{5};

j = 1;
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
        title('OSNAP 1 - SALINITY')
        
        figure(20012)
        hold on; box on;
        plot(JG, Tfs1(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 1 - TEMPERATURE')
        
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

% OSNAP 2 (PE399 --> DY053)


disp('---------  OSNAP 2 (PE399 --> DY053) ---------')
fileID2 = fopen([boundarydir, moor{2}, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
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

% OSNAP 3 (DY053 --> DY078)


disp('---------  OSNAP 3 (DY053--> DY078) ---------')
fileID3 = fopen([boundarydir, moor{3}, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
file3_data = textscan(fileID3, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
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
% OSNAP 4 (DY078 --> AR30)


disp('---------  OSNAP 4 (DY078--> AR30) ---------')
fileID4 = fopen([boundarydir, moor{4}, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
file4_data = textscan(fileID4, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
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

% OSNAP 5 (AR30 --> DY120)

disp('---------  OSNAP 5 (AR30--> DY120) ---------')
fileID5 = fopen([boundarydir, moor{5}, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
file5_data = textscan(fileID5, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID5);

mooring5 = file5_data{1};
sn5      = file5_data{2};
mc5      = file5_data{3};
z5       = file5_data{4};
lon5     = file5_data{5};

T5 = zeros(length(sn5), length(JG));
P5 = zeros(length(sn5), length(JG));
S5 = zeros(length(sn5), length(JG));

i = 1; j = 1;
for i = 1: length(sn5)
    
    infile = [hydrodir, mooring5{i,:}, '_grid.mat'];
    load(infile);
    jdnew = datenum(gregorian(jd));    
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jdnew(1:end-1), Sfs(mc5(i),1:end-1)', JG)';
    temp        = interp1(jdnew(1:end-1), Tfs(mc5(i),1:end-1)', JG)';
    pressure    = interp1(jdnew(1:end-1), Pfs(mc5(i),1:end-1)', JG)';
    Sfs5(j,:) = salinity;
    Pfs5(j,:) = pressure;
    Tfs5(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(50021)
        hold on; box on;
        plot(JG, Sfs5(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 5 - SALINITY')
        
        figure(50022)
        hold on; box on;
        plot(JG, Tfs5(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 5 - TEMPERATURE')
        
        figure(500220)
        hold on; box on;
        plot(JG, sw_pden(Sfs5(i, :),Tfs5(i, :),Pfs5(i, :),0)-1000, col{:,i},'LineWidth',2);
        title('OSNAP 5 - POTENTIAL DENSITY')      
    
        figure(50023)
        hold on; box on;
        plot(JG, Pfs5(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 5 - PRESSURE')
    end
    
end


% OSNAP 6 (DY120 -->JC238)

disp('---------  OSNAP 6 (DY120 -->JC238) ---------')
fileID6 = fopen([boundarydir, moor{6}, '.dat']);
delimiter = {'\t',' '};
startRow = 6;
formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
file6_data = textscan(fileID6, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID6);

mooring6 = file6_data{1};
sn6      = file6_data{2};
mc6      = file6_data{3};
z6       = file6_data{4};
lon6     = file6_data{5};

T6 = zeros(length(sn6), length(JG));
P6 = zeros(length(sn6), length(JG));
S6 = zeros(length(sn6), length(JG));

i = 1; j = 1;
for i = 1: length(sn6)
    
    infile = [hydrodir, mooring6{i,:}, '_grid.mat'];
    load(infile);
    jdnew = datenum(gregorian(jd));    
    sampling_rate = round(1./median(diff(jd_grid))); %nominal sampling rate [per day]
    salinity    = interp1(jdnew(1:end-1), Sfs(mc6(i),1:end-1)', JG)';
    temp        = interp1(jdnew(1:end-1), Tfs(mc6(i),1:end-1)', JG)';
    pressure    = interp1(jdnew(1:end-1), Pfs(mc6(i),1:end-1)', JG)';
    Sfs6(j,:) = salinity;
    Pfs6(j,:) = pressure;
    Tfs6(j,:) = temp;
    
    j = j + 1;
    
    if mc_check_plot
        figure(60021)
        hold on; box on;
        plot(JG, Sfs6(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 6 - SALINITY')
        
        figure(60022)
        hold on; box on;
        plot(JG, Tfs6(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 6 - TEMPERATURE')
        
        figure(600220)
        hold on; box on;
        plot(JG, sw_pden(Sfs6(i, :),Tfs6(i, :),Pfs6(i, :),0)-1000, col{:,i},'LineWidth',2);
        title('OSNAP 6 - POTENTIAL DENSITY')      
    
        figure(60023)
        hold on; box on;
        plot(JG, Pfs6(i, :), col{:,i},'LineWidth',2);
        title('OSNAP 6 - PRESSURE')
    end
    
end

close all
%% 2. CONCATENATE AND ORDER THE MATRICES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This step adds alll the matrices for the deployments together to form
%  large data sets.  These are then sorted into pressure order at every
%  time step before the gridding takes place.

%if mc_check_plot
PDENfs1 = sw_pden(Sfs1,Tfs1,Pfs1,0) - 1000;
PDENfs2 = sw_pden(Sfs2,Tfs2,Pfs2,0) - 1000;  
PDENfs3 = sw_pden(Sfs3,Tfs3,Pfs3,0) - 1000;
PDENfs4 = sw_pden(Sfs4,Tfs4,Pfs4,0) - 1000;
PDENfs5 = sw_pden(Sfs5,Tfs5,Pfs5,0) - 1000;
PDENfs6 = sw_pden(Sfs6,Tfs6,Pfs6,0) - 1000;

    
figure(11)   %  graph of the data to show that it is all there!
clf;
subplot(2,1,1);
hold on; box on;
plot(JG , Tfs1, 'color',[0 0.4470 0.7410])
plot(JG , Tfs2, 'color',[0.8500 0.3250 0.0980])
plot(JG , Tfs3, 'color',[0.9290 0.6940 0.1250])
plot(JG , Tfs4, 'color',[0.4940 0.1840 0.5560])
plot(JG , Tfs5, 'color',[0.4660 0.6740 0.1880])
plot(JG , Tfs6, 'color',[0.3010 0.7450 0.9330])
ylabel('C')
datetick
title('TEMPERATURE')

subplot(2,1,2);
hold on; box on;
plot(JG , Sfs1, 'color',[0 0.4470 0.7410])
plot(JG , Sfs2, 'color',[0.8500 0.3250 0.0980])
plot(JG , Sfs3, 'color',[0.9290 0.6940 0.1250])
plot(JG , Sfs4, 'color',[0.4940 0.1840 0.5560])
plot(JG , Sfs5, 'color',[0.4660 0.6740 0.1880])
plot(JG , Sfs6, 'color',[0.3010 0.7450 0.9330])  
ylabel('SAL.')
datetick
title('SALINITY')
    
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
if ~exist([grdatdir 'otherfigure' filesep],'dir');mkdir([grdatdir 'otherfigure' filesep]);end
print('-dpng',[grdatdir 'otherfigure' filesep  'RTWBmerged_beforegrid_check1'])

figure(10)   %  graph of the data to show that it is all there!
clf;
subplot(2,1,1);
hold on; box on;
plot(JG , Pfs1, 'k.')
plot(JG , Pfs2, 'b.')
plot(JG , Pfs3, 'g.') 
plot(JG , Pfs4, 'r.')
plot(JG , Pfs5, 'm.')
plot(JG , Pfs6, 'm.')
ylabel('dbar')
datetick
title('PRES')   

subplot(2,1,2);
hold on; box on;
plot(JG , PDENfs1, 'k.')
plot(JG , PDENfs2, 'b.')
plot(JG , PDENfs3, 'g.')  
plot(JG , PDENfs4, 'r.') 
plot(JG , PDENfs5, 'm.') 
plot(JG , PDENfs6, 'm.') 
ylabel('kg/m3.')
datetick
title('POT. DENS')   

set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
print('-dpng',[grdatdir 'otherfigure' filesep  'RTWBmerged_beforegrid_check2'])

% all the matrices for the deployments stacked together
Pfs     = [Pfs1;Pfs2;Pfs3;Pfs4;Pfs5;Pfs6];
Sfs     = [Sfs1;Sfs2;Sfs3;Sfs4;Sfs5;Sfs6];
Tfs     = [Tfs1;Tfs2;Tfs3;Tfs4;Tfs5;Tfs6];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KB outputs stagged files
z_stacked = [50;75;100;250;500;750;1000;1250;1500;1650;1775];

%%% add NaN where instrument is missing
fvs = NaN(1,length(JG));
t_st1 = [fvs;Tfs1(1,:);fvs;Tfs1(2:4,:);fvs;Tfs1(5:6,:);fvs;Tfs1(7,:)];
s_st1 = [fvs;Sfs1(1,:);fvs;Sfs1(2:4,:);fvs;Sfs1(5:6,:);fvs;Sfs1(7,:)];
p_st1 = [fvs;Pfs1(1,:);fvs;Pfs1(2:4,:);fvs;Pfs1(5:6,:);fvs;Pfs1(7,:)];

t_st2 = [fvs;Tfs2(1,:);fvs;Tfs2(2:7,:);fvs;Tfs2(8,:)];
s_st2 = [fvs;Sfs2(1,:);fvs;Sfs2(2:7,:);fvs;Sfs2(8,:)];
p_st2 = [fvs;Pfs2(1,:);fvs;Pfs2(2:7,:);fvs;Pfs2(8,:)];

t_st3 = [fvs;fvs;Tfs3(1,:);fvs;Tfs3(2:6,:);fvs;Tfs3(7,:)];
s_st3 = [fvs;fvs;Sfs3(1,:);fvs;Sfs3(2:6,:);fvs;Sfs3(7,:)];
p_st3 = [fvs;fvs;Pfs3(1,:);fvs;Pfs3(2:6,:);fvs;Pfs3(7,:)];

t_st4 = [Tfs4(1,:);fvs;Tfs4(2:8,:);fvs;Tfs4(9,:)];
s_st4 = [Sfs4(1,:);fvs;Sfs4(2:8,:);fvs;Sfs4(9,:)];
p_st4 = [Pfs4(1,:);fvs;Pfs4(2:8,:);fvs;Pfs4(9,:)];

t_st5 = [Tfs5(1,:);fvs;Tfs5(2:10,:)];
s_st5 = [Sfs5(1,:);fvs;Sfs5(2:10,:)];
p_st5 = [Pfs5(1,:);fvs;Pfs5(2:10,:)];

t_st6 = [Tfs6(1,:);fvs;Tfs6(2:10,:)];
s_st6 = [Sfs6(1,:);fvs;Sfs6(2:10,:)];
p_st6 = [Pfs6(1,:);fvs;Pfs6(2:10,:)];

%%% add NaN where instrument is missing
Tstacked = stack_vars(t_st1,t_st2,t_st3,t_st4,t_st5,t_st6);
Sstacked = stack_vars(s_st1,s_st2,s_st3,s_st4,s_st5,s_st6);
Pstacked = stack_vars(p_st1,p_st2,p_st3,p_st4,p_st5,p_st6);

save([grdatdir outputfile_stacked '.mat'],...
    'Tstacked','Sstacked','Pstacked','JG','z_stacked')

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

clear Pfs1 Pfs2 Pfs3 Pfs4 Pfs5 Pfs6
clear Sfs1 Sfs2 Sfs3 Sfs4 Sfs5 Sfs6
clear Tfs1 Tfs2 Tfs3 Tfs4 Tfs5 Tfs6
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. GRIDDING

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
    outputfile = ['RTWB_merg_' TS_CLIMATOLOGY_TP ' ' TS_CLIMATOLOGY '_' lastyeardata ];
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
    
    outputfile = ['RTWB_merg_linear_interp' '_'  lastyeardata ];
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

TGfs(pgg>max_depth,:)=NaN;
SGfs(pgg>max_depth,:)=NaN;

disp(['saving: ' outputfile '.mat'])
% Alloacte variables into a structure
RTWB_merg.JG     = JG;
RTWB_merg.Tfs    = Tfs;
RTWB_merg.Sfs    = Sfs;
RTWB_merg.Pfs    = Pfs;
RTWB_merg.PGfs    = pgg';    
RTWB_merg.TGfs   = TGfs;
RTWB_merg.SGfs   = SGfs;

if gridding == 2 % using climatological profiles    
    RTWB_merg.TS_CLIMATOLOGY   = TS_CLIMATOLOGY;
    RTWB_merg.TS_CLIMATOLOGY_TP= TS_CLIMATOLOGY_TP;
    RTWB_merg.TS_CLIMATOLOGY_NAME= TS_CLIMATOLOGY_NAME;
    RTWB_merg.clim_file        = TSclim;
end    
%% 4. Despike time series


stddy_tol  = 10; % 4
std_win    = 3.5; % 3.5 * std of the time series
[nloop]    = 5; % 5
graphics   = 'y';
temp = []; salinity = [];

for i = 1 : length(TGfs(:,1))   % loop through each depth level

    [temp(i,:),dx,ndx] = ddspike(TGfs(i,:),[-std_win*nanstd(TGfs(i,:)),...
        std_win*nanstd(TGfs(i,:))],stddy_tol,[nloop],'y',NaN);
    [salinity(i,:),dx,ndx] = ddspike(SGfs(i,:),[-std_win*nanstd(SGfs(i,:)),...
        std_win*nanstd(SGfs(i,:))],stddy_tol,[nloop],'y',NaN);
end

[m,n] = size(TGfs);
TG_west = NaN * ones(m,n); SG_west = NaN * ones(m,n);

% GDM, 9/4/2013
% Need to interpolate horizontally over nans but...
% Don't want to create fake values above knocked down moorings
% So select a depth below which, interpoation happens

idepth = depthminforhoriz_interp; % multiples of 20
I = find(pgg == idepth);

% copy the top idepth into the new file
% no temporal interpolation
SG_west([1:I-1],:)=SGfs([1:I-1],:);
TG_west([1:I-1],:)=TGfs([1:I-1],:);

%    i = 1; j = 1;
%    for i = 1: 6 % top 100m based on a grid of 20dbar
%        SG_west(j,:) = SGfs(i, :);
%        TG_west(j,:) = TGfs(i, :);
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
    SG_west(i,:) = interp1(JG(is), salinity(i, is), JG);
    % interpolate in time over the missing data
    TG_west(i,:) = interp1(JG(it), temp(i, it), JG);

    %j = j + 1;
end


%%
 % KB, 27/3/2024
    % Added code to remove horizontal interpolation when mooring was knocked down 
    % This was not needed before because data was vertical extrapolated to surface 
    % first which is now done in later step.

for i = 1:length(JG);
    it = find(~isnan(TGfs(:,i)),1,'first');
    is = find(~isnan(SGfs(:,i)),1,'first');

    TG_west(1:it-1,i)=nan;
    SG_west(1:is-1,i)=nan;
end

%% 5. Plotting

figure;
% spurious values created in the interpolation?
hold on; grid on; box on;
plot(SG_west,TG_west,'k')
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
contourf(JG , pgg, TG_west,8); axis ij
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
contourf(JG , pgg, SG_west,8); axis ij
caxis([min(SGfs(:)) max(SGfs(:))])
datetick; ylabel('SALINITY')
title('AFTER DESPIKING AND INTERPOLATION')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
print('-dpng',[grdatdir outputfile '_salinity'])


PDENGfs = sw_pden(SGfs,TGfs,pgg',0);
pden = sw_pden(salinity,temp,pgg',0);   
PDENG_west = sw_pden(SG_west,TG_west,pgg',0);   
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
contourf(JG , pgg,PDENG_west,10); axis ij
  caxis([min(PDENGfs(:)) max(PDENGfs(:))])
datetick; ylabel('POT. DENS.')
title('AFTER DESPIKING AND INTERPOLATION')

 set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)   
print('-dpng',[grdatdir outputfile '_potdens'])

%  MAKE MERGED DAT ONLY FIGURE 

% TEMPERATURE
figure;
ax(1)=subplot(2,1,1)

x=JG;
y=pgg;
z=TG_west;
[c,h]=contourf(x,y,z,100,'LineColor','none');
caxis([ceil(min(min(TG_west))) floor(max(max(TG_west)))]);
C = colorbar;
C.Label.String = 'Conservative Temperature (^{o} C)';
C.Ticks=[floor(min(min(TG_west))):1:ceil(max(max(TG_west)))];
C.TickLength=0.03;
hold on
[c,h]=contour(x , y, z,[2 4 6 8 10 12 14],'LineColor','k','LineWidth',0.25);
clabel(c,h,'LabelSpacing',500,'Color','w')
xlabel(gca,'Year');
ylabel(gca,'Pressure (m)');
[Y,~,~]=datevec(JG);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
axis ij
title('Rockall Trough Western Boundary Gridded Temperature')
cmap=cmocean('Thermal');
colormap(ax(1),cmap);

% SALINITY
ax(2)=subplot(2,1,2)
x=JG;
y=pgg;
z=SG_west;
[cc,hh]=contourf(x,y,z,[35:0.1:35.8],'LineColor','none');
C = colorbar;
cmap=cmocean('Haline');
colormap(ax(2),cmap);
C.Label.String = 'Absolute Salinity (g kg^{-1})';
caxis([35 35.8]);
C.Ticks=[35:0.1:35.8];
C.TickLength=0.03;
hold on
[c,h]=contour(JG , pgg, SGfs,[35 35.2 35.4 35.6 35.8],'LineColor','k','LineWidth',0.25); 
clabel(c,h,'LabelSpacing',500,'Color','w')
xlabel(gca,'Year');
ylabel(gca,'Pressure (m)');
% SORT LABELS
[Y,~,~]=datevec(JG);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
axis ij
title('Rockall Trough Western Boundary Gridded Salinity')

width=35; height=20; FS=12; FN='Arial';
set(gca,'fontsize', FS, 'FontName',FN);
set(gcf,'units','centimeters','position',[5 5 width height])
print('-dpng',[grdatdir outputfile])

%% Duplicate values to the surface 
% Originally writen by Neil Fraser for volume transport code, but moved upstream 
% for simplicity.

% alocate final product
RTWB_merg.TGfs2= TG_west;
RTWB_merg.SGfs2= SG_west;

% % Remove padding NaNs from mooring data
% idx = ~isnan(nanmean(RTWB_merg.SGfs2,1)) ;
% RTWB_merg.Tfs(:,~idx) = [];
% RTWB_merg.Sfs(:,~idx) = [];
% RTWB_merg.Pfs(:,~idx) = [];
% RTWB_merg.TGfs(:,~idx) = [];
% RTWB_merg.SGfs(:,~idx) = [];
% RTWB_merg.TGfs2(:,~idx) = [];
% RTWB_merg.SGfs2(:,~idx) = [];
% RTWB_merg.JG(~idx) = [];

% Grid dimensions
[NZ,NT] = size(RTWB_merg.TGfs2);

% Duplicate SA and CT to surface
% for n = 1:NT
%    k_RTWB = find(~isnan(RTWB_merg.TGfs2(:,n)),1,'first');
%    RTWB_merg.TGfs2(1:k_RTWB-1,n) = RTWB_merg.TGfs2(k_RTWB,n);
%    RTWB_merg.SGfs2(1:k_RTWB-1,n) = RTWB_merg.SGfs2(k_RTWB,n);
%end


% Alloacte variables into a structure

RTWB_merg.comment{1,1}= 'JG -- julian day';
RTWB_merg.comment{2,1}= 'Tfs -- original stacked Temperature data from the deployments';
RTWB_merg.comment{3,1}= 'Sfs -- original stacked salinity data from the deployments ';
RTWB_merg.comment{4,1}= 'Pfs -- original stacked Pressure data from the deployments';    
RTWB_merg.comment{5,1}= 'PGfs -- pressure grid ';
RTWB_merg.comment{6,1}= 'TGfs -- temperature interpolated onto the pressure grid (PGfs)';
RTWB_merg.comment{7,1}= 'SGfs -- salinity interpolated onto the pressure grid (PGfs)';
RTWB_merg.comment{8,1}= 'TGfs2 -- temperature interpolated onto the time grid (JG) after despiking';  
RTWB_merg.comment{9,1}= 'SGfs2 -- salinity interpolated onto the time grid (JG) after despiking';  

save([grdatdir outputfile],'RTWB_merg');
 
%% functions
function Ustacked = stack_vars(U1,U2,U3,U4,U5,U6)
    Ustacked = cat(3,U1,U2,U3,U4,U5,U6);
    Ustacked = sum(Ustacked,3,'omitnan');
    Ustacked(Ustacked==0)=NaN;
end