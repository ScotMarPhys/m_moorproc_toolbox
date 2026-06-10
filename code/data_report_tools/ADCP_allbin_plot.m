function ADCP_allbin_plot(moor,varargin)
%
% Function to load and plot ADCP data at the rodb file format 
% 
%
% Function to load and plot ADCP data at the rodb file format 
%
% function ADCP_allbin_plot('moor','proclvl','procpath')
%
% required inputs:-
%   moor: complete mooring name as string. e.g. 'wb1_1_200420'
%
% optional inputs:-
%   procpath: can specify exact path to proc directory if not using 
%           standard data paths. 
%           e.g. '/Volumes/noc/mpoc/hydro/rpdmoc/rapid/data/moor/proc/'
%   proclvl: can specify level of processing of the data to plot. 
%           e.g. 'proclvl','2': will plot the .use file ; 'proclvl','3' will plot the .microcat and .edt files
%
% functions called:-
%   rodbload, julian, auto_filt
%   from .../exec/moor/tools and .../exec/moor/rodb paths
%   suptitle.m downloaded from Mathworks file exchange website 
%
%
% 07/10/16 - Loic Houpert, SAMS

if nargin <1
    help  plot_stacked
    return
end

a=strmatch('procpath',varargin,'exact');
if a>0
    procpath=char(varargin(a+1));
else
    data_report_tools_dir=which('data_report_tools');
    b=strfind(data_report_tools_dir,'/');
    data_report_tools_dir=data_report_tools_dir(1:b(end));
    procpath=[data_report_tools_dir '../../../moor/proc/']; %DR changed to be relative paths now that data_report_tools are on the network 19/2/12
end

a=strmatch('proclvl',varargin,'exact');
if a>0
    proclvlstr0=char(varargin(a+1));
    proclvl   = str2num(proclvlstr0);
else
    proclvl=2;
    proclvlstr0 = num2str(proclvl);
end

if isunix
    infofile=[procpath,moor,'/',moor,'info.dat'];
elseif ispc
    infofile=[procpath,moor,'\',moor,'info.dat'];
end

% Load vectors of mooring information
% id instrument id, sn serial number, z nominal depth of each instrument
% s_t, e_t, s_d, e_d start and end times and dates
% lat lon mooring position, wd corrected water depth (m)
% mr mooring name
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,...
    'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');


% JULIAN Convert Gregorian date to Julian day.
% JD = JULIAN(YY,MM,DD,HH) or JD = JULIAN([YY,MM,DD,HH]) returns the Julian
% day number of the calendar date specified by year, month, day, and decimal
% hour.
% JD = JULIAN(YY,MM,DD) or JD = JULIAN([YY,MM,DD]) if decimal hour is absent,
% it is assumed to be zero.
% Although the formal definition holds that Julian days start and end at
% noon, here Julian days start and end at midnight. In this convention,
% Julian day 2440000 began at 00:00 hours, May 23, 1968.
jd_start = julian([s_d' hms2h([s_t;0]')]);
jd_end   = julian([e_d' hms2h([e_t;0]')]);

disp(['z : instrument id : serial number'])
for i = 1:length(id)
    disp([z(i),id(i),sn(i)])
end

% Possible ADCP codes - taken from IMP moorings package
iiADCP = find((id>=319) & (id <=328));
vecADCP = sn(iiADCP);

% depths(:,1) = id([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG]);
% depths(:,2) = z([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG]);
depths(:,1) = id(iiADCP);
depths(:,2) = z(iiADCP);
iiiADCP=find(depths(:,1)>=319 & depths(:,1)<=328);
iii=[iiiADCP];


plot_string={};

if proclvl==2
    fileext = '.use';
    extfig  = 'proclvl_2';
elseif proclvl==3
    fileext = '.edt';
    extfig  = 'proclvl_3';    
end 
% -----------------------------------
% START OF READING IN INSTRUMENT DATA
% -----------------------------------
%--------------------------------------
% Now read in ADCP data if required
%--------------------------------------
if iiADCP>0
    j=1;
    
    % loop to read one file at a time

    for i=1:length(vecADCP);
       serialno = vecADCP(i);
       disp('*************************************************************')
       disp(['Reading ADDCP - ',num2str(serialno)])
       disp('*************************************************************')

    
    % first determine how many bin files are to be processed
    % trying to do automatically
    num_bins=dir([procpath,moor,'/adp/',moor,'_',num2str(vecADCP(i)),'_bin*' fileext])
    num_bins=length(num_bins)
    num_bins
    
    namefigadcp1 = [extfig '_' moor '_adcp_' num2str(vecADCP(i)) '_uv_surf'];
    namefigadcp2 = [extfig '_' moor '_adcp_' num2str(vecADCP(i)) '_pgood_surf'];    
    
    for j=1:num_bins % loop for total number of bins
        
        columns = ['YY:MM:DD:HH:Z:T:U:V:W:HDG:PIT:ROL:CS:CD:BEAM1SS:BEAM2SS:BEAM3SS'...
            ':BEAM4SS:BEAM1COR:BEAM2COR:BEAM3COR:BEAM4COR:EV:BEAM1PGP:BEAM2PGP:BEAM3PGP:BEAM4PGP'];
        
        indep  = z(i);
        
        if j<=9
            infile  = [procpath,moor,'/adp/',moor,'_',num2str(vecADCP(i)),'_bin0',num2str(j),fileext];
        else
            infile  = [procpath,moor,'/adp/',moor,'_',num2str(vecADCP(i)),'_bin',num2str(j),fileext];
        end
                
        if exist(infile,'file')==0
            disp(['infile: ' infile ' does not exist.'])

        elseif exist(infile,'file')   > 0 

            [YY,MM,DD,HH,z,t,u,v,w,hdg,pit,rol,spd,direction,Amp1,Amp2,Amp3,Amp4,...
                Beam1Cor,Beam2Cor,Beam3Cor,Beam4Cor,err,PG1,PG2,PG3,PG4] = ...
                rodbload(infile,[columns]);
            
           bad_data=find(t==-9999); t(bad_data)=NaN;
           bad_data=find(z==-9999); z(bad_data)=NaN;
           bad_data=find(u==-9999); u(bad_data)=NaN;
           bad_data=find(v==-9999); v(bad_data)=NaN;
           bad_data=find(w==-9999); w(bad_data)=NaN;
           bad_data=find(hdg==-9999); hdg(bad_data)=NaN;
           bad_data=find(pit==-9999); pit(bad_data)=NaN;
           bad_data=find(rol==-9999); rol(bad_data)=NaN;          
           bad_data=find(err==-9999); err(bad_data)=NaN;       
           bad_data=find(PG1==-9999); PG1(bad_data)=NaN;        
           bad_data=find(PG2==-9999); PG2(bad_data)=NaN;   
           bad_data=find(PG3==-9999); PG3(bad_data)=NaN;        
           bad_data=find(PG4==-9999); PG4(bad_data)=NaN;   
           
            ADCPdata(i).time(j,:) = datenum(YY,MM,DD,HH,0*HH,0*HH); %julian(YY,MM,DD,HH);
            ADCPdata(i).z(j,:) = z;
            ADCPdata(i).t(j,:) = t;          
            ADCPdata(i).u(j,:) = u;
            ADCPdata(i).v(j,:) = v;       
            ADCPdata(i).w(j,:) = w;   
            ADCPdata(i).err(j,:) = err;    
            ADCPdata(i).hdg(j,:) = hdg;       
            ADCPdata(i).pit(j,:) = pit;   
            ADCPdata(i).rol(j,:) = rol;                
            ADCPdata(i).PG1(j,:) = PG1;       
            ADCPdata(i).PG2(j,:) = PG2;    
            ADCPdata(i).PG3(j,:) = PG3;       
            ADCPdata(i).PG4(j,:) = PG4;              

        end
        
        
    end

            ADCPdata(i).name = ['NOR_' num2str(serialno)];
            ADCPdata(i).moor = moor;            
    
    
           j=j+1;
    end
end

% ------------------------------
% Plotting section
% ------------------------------

%----------------------------------------------------
% Graphical parameter
vsblfig = 'on';
zsc = get(0,'MonitorPositions');
scrsz = [1 1 1900 1400]; % zsc(1,:);
figpos = scrsz/60; %[0 0 27 21];
fs1 = 12;
fs2=10;
set(0,'DefaultAxesFontName', 'Helvetica')
set(0,'DefaultAxesFontSize', fs1)
set(0,'DefaultTextFontname', 'Helvetica')
set(0,'DefaultTextFontSize', fs2)
%-----------------------------------------------------
cmapvel = lansey;
    

proc = 1;    
plotylim = [min(-ADCPdata(proc).z(:)) nanmean(-ADCPdata(proc).z(1,:))];

fig=figure('visible',vsblfig,'position',scrsz);
set(fig,'PaperUnits','centimeters','PaperOrientation','portrait',... 
				    'Paperposition',figpos)
subplot(2,1,1)
surf(ADCPdata(proc).time, -ADCPdata(proc).z,ADCPdata(proc).u);
set(gca,'ylim', plotylim)
shading flat
view([0 90])
caxis([-50 50])
title('U (cm.s-1)')
datetick
colormap(cmapvel)
colorbar

subplot(2,1,2)
surf(ADCPdata(proc).time, -ADCPdata(proc).z,ADCPdata(proc).v);
set(gca,'ylim', plotylim)
shading flat
view([0 90])
caxis([-50 50])
title('V (cm.s-1)')
datetick
colormap(cmapvel)
colorbar

print('-dpng',namefigadcp1)

    
fig=figure('visible',vsblfig,'position',scrsz);
set(fig,'PaperUnits','centimeters','PaperOrientation','portrait',... 
				    'Paperposition',figpos)


subplot(3,1,1)
surf(ADCPdata(proc).time, -ADCPdata(proc).z,ADCPdata(proc).PG4+ADCPdata(proc).PG4);
set(gca,'ylim', plotylim)
shading flat
view([0 90])
caxis([0 100])
title('PG4 - pgood 4 beam solutions')
datetick
colormap(gca,parula)
colorbar
    
subplot(3,1,2)
surf(ADCPdata(proc).time, -ADCPdata(proc).z,ADCPdata(proc).PG4+ADCPdata(proc).PG1);
set(gca,'ylim', plotylim)
shading flat
view([0 90])
caxis([0 100])
title('PG1 - pgood 3 beam solutions')
datetick
colormap(gca,parula)
colorbar

subplot(3,1,3)
surf(ADCPdata(proc).time, -ADCPdata(proc).z,ADCPdata(proc).err);
set(gca,'ylim', plotylim)
shading flat
view([0 90])
caxis([0 5])
title('Error Velocity (cm.s-1)')
datetick
colormap(gca,parula)
colorbar

print('-dpng',namefigadcp2)

    
end
