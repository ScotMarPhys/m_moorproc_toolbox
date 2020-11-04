% velocity_grid_adcp.m
% Interpolate current meter data (u,v,press)
%  1. 2-day low-pass ; 2. intepolate onto a 12-hour grid ; 3. interpolate
%  vertically using a variety of methods ; 4. create mean u,v
% Jan 2017 - Adapted for OSNAP moorings (by L. Houpert)

clear all ; close all ;
prognam = 'velocity_grid_adcp' ;

moor = 'rtadcp1_01_2014' ;
jg_start = datenum(2014,7,19,00,00,00);  jg_end = datenum(2015,6,20,00,00,00);

dnumi              = jg_start: 0.5: jg_end; % full time series using 2 samples per day
pgrid = 0:20:2000; % depths in 20dbar bins
gap_max      = 10; % no more than 10 days of nan in a row

basedir      = '/media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap';
%basedir      = '/home/sa02lh/Data/Dropbox/Work/Postdoc_OSNAP/OSNAP_mooring/backup_mdrive';
p_hydrogrid.moor   = moor;   
mooringpath  = [basedir '/data/moor/proc/' moor ];
out_path     = [basedir '/data/moor/proc/velocity_grid/'];

do_plot = 0;

% --- get moring information from infofile
infofile =[mooringpath '/' moor 'info.dat'];
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

vec=find((id>=319) & (id <=328)) % Possible ADCP codes - taken from IMP moorings package
sn=sn(vec);
z = z(vec);

num_bins=dir([moor,'/adp/*',num2str(sn) '*.raw']);
num_bins=length(num_bins);

columns = ['YY:MM:DD:HH:Z:T:U:V:W:HDG:PIT:ROL:CS:CD:BEAM1SS:BEAM2SS:BEAM3SS'...
    ':BEAM4SS:BEAM1COR:BEAM2COR:BEAM3COR:BEAM4COR:EV:BEAM1PGP:BEAM2PGP:BEAM3PGP:BEAM4PGP'];
% last good data bin 37
% for proc = 1: 1 : num_bins % loop for total number of bins

ufi = [];
vfi = [];
wfi = [];
pfi = [];
ufii  = [];
vfii  = [];
wfii  = [];
cc = 0;
for proc=1:32 ;
    if (proc<10);
        infile  = [mooringpath,'/adp/',moor,'_',num2str(sn),'_bin0',sprintf('%d',proc),'.edt'];
    else
        infile  = [mooringpath,'/adp/',moor,'_',num2str(sn),'_bin',sprintf('%d',proc),'.edt'];
    end
    
    [YY,MM,DD,HH,z,t,u,v,w,heading,pitch,roll,spd,direction,Amp1,Amp2,Amp3,Amp4,...
        Beam1Cor,Beam2Cor,Beam3Cor,Beam4Cor,err,PG1,PG2,PG3,PG4] = ...
        rodbload(infile,[columns]);
    dnum = datenum(YY,MM,DD,HH,0,0);
    
    if isempty(u) % no instrument data (ex: rtwb2_01_2014)
        continue
    end
    cc = cc+1;   

    bad_data=find(t==-9999); t(bad_data)=nan;
    bad_data=find(z==-9999); z(bad_data)=nan;
    bad_data=find(u==-9999); u(bad_data)=nan;
    bad_data=find(v==-9999); v(bad_data)=nan;
    bad_data=find(w==-9999); w(bad_data)=nan;
    
    p=z;  
    
    unan_sum = sum(isnan(u));
    vnan_sum = sum(isnan(v));
    wnan_sum = sum(isnan(w));
    pnan_sum = sum(isnan(p));
    % filtering parameter
    sr  = median(diff(dnum)); % sampling interval
    co = sr/((1/sr)*2); % 2-day low-pass
    
    % 2-day low pass
    unnan = find(~isnan(u));
    vnnan = find(~isnan(v));
    wnnan = find(~isnan(w));    
    pnnan = find(~isnan(p));
    uf0(unnan)  = auto_filt(u(unnan),sr,co,'low',4);
    vf0(unnan)  = auto_filt(v(unnan),sr,co,'low',4);
    wf0(unnan)  = auto_filt(w(unnan),sr,co,'low',4);      
    pf0(unnan)  = auto_filt(p(unnan),sr,co,'low',4);
    uf = interp1(dnum(unnan),uf0(unnan),dnum,'linear');
    vf = interp1(dnum(unnan),vf0(unnan),dnum,'linear');
    wf = interp1(dnum(unnan),wf0(unnan),dnum,'linear');        
    pf = interp1(dnum(unnan),pf0(unnan),dnum,'linear');
    if unan_sum*sr > gap_max
            gapI        = gap_mark(u,gap_max,1/sr);
            uf(gapI) = NaN;
    end
    if vnan_sum*sr > gap_max
            gapI        = gap_mark(v,gap_max,1/sr);
            vf(gapI) = NaN;
    end    
    if wnan_sum*sr > gap_max
            gapI        = gap_mark(w,gap_max,1/sr);
            wf(gapI) = NaN;
    end       
    if pnan_sum*sr > gap_max
            gapI        = gap_mark(p,gap_max,1/sr);
            pf(gapI) = NaN;
    end  
    
    % put data onto a 12-hour grid
    ufi(cc,:) = interp1(dnum,uf,dnumi,'linear');
    vfi(cc,:) = interp1(dnum,vf,dnumi,'linear');
    wfi(cc,:) = interp1(dnum,wf,dnumi,'linear');    
    pfi(cc,:) = interp1(dnum,pf,dnumi,'linear');
    
    if do_plot ;
        % Plot individual records
        figure(proc);clf;
        [p1]=subplot(3,1,1);hold on;grid on;set(gca,'layer','top');
        plot(dnum,u,'k');
        plot(dnum,uf,'y','linewidth',2);
        plot(dnumi,ufi(proc,:),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);
        datetick('x',24)
        ylim([-25 25]);
        ylabel('u [cm/s]');
        title(['Bin : ',sprintf('%d',proc),' : Depth : ',sprintf('%5.1f',nanmean(p)), 'dbar']);
        
        [p2]=subplot(3,1,2);hold on;grid on;set(gca,'layer','top');
        plot(dnum,v,'k');
        plot(dnum,vf,'y','linewidth',2);
        plot(dnumi,vfi(proc,:),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5)
        datetick('x',24)
        ylim([-50 50]);
        ylabel('v [cm/s]');
        
        [p3]=subplot(3,1,3);hold on;grid on;set(gca,'layer','top');
        plot(dnum,p,'k');
        plot(dnum,pf,'y','linewidth',2);
        plot(dnumi,pfi(proc,:),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5)
        datetick('x',24)
        ylim([100 750]);
        ylabel('depth [m]');
        set(gca,'YDIR','reverse');
        linkaxes([p1 p2 p3],'x');
    end
    
end


%% Now onto vertical interpolation of the subsampled, low-pass data 
% matices: dnumi, ufi(6,length(dnumi)), vfi, pfi
ufii_pchip  = nan(length(pgrid),length(dnumi));
vfii_pchip  = nan(length(pgrid),length(dnumi));
wfii_pchip  = nan(length(pgrid),length(dnumi));
ufii_linear = nan(length(pgrid),length(dnumi));
vfii_linear = nan(length(pgrid),length(dnumi));
wfii_linear = nan(length(pgrid),length(dnumi));
ufii_akima  = nan(length(pgrid),length(dnumi));
vfii_akima  = nan(length(pgrid),length(dnumi));
wfii_akima  = nan(length(pgrid),length(dnumi));


maxdepth = max(nanmean(pfi,2));
% interpolation
for i=1:length(dnumi) ;
    iok = find(~isnan(ufi(:,i)));
    if length(iok)>2    
    ufii_pchip(:,i) = interp1(pfi(iok,i),ufi(iok,i),pgrid,'pchip','extrap') ;
    vfii_pchip(:,i) = interp1(pfi(iok,i),vfi(iok,i),pgrid,'pchip','extrap') ;
    wfii_pchip(:,i) = interp1(pfi(iok,i),wfi(iok,i),pgrid,'pchip','extrap') ;
    % mask extrapolated data below the deepest measurement (== bottom)
    inodata = find(pgrid>maxdepth);
    ufii_pchip(inodata,i) = nan;
    vfii_pchip(inodata,i) = nan;    
    wfii_pchip(inodata,i) = nan;      
    end
end

% This is a really bad method
% for i=1:length(dnumi) ;
%     ufii_spline(:,i) = interp1(pfi(:,i),ufi(:,i),pgrid,'spline','extrap') ;
%     vfii_spline(:,i) = interp1(pfi(:,i),vfi(:,i),pgrid,'spline','extrap') ;
% end

for i=1:length(dnumi) ;
    iok = find(~isnan(ufi(:,i)));
    if length(iok)>2    
    ufii_linear(:,i) = interp1(pfi(iok,i),ufi(iok,i),pgrid,'linear','extrap') ;
    vfii_linear(:,i) = interp1(pfi(iok,i),vfi(iok,i),pgrid,'linear','extrap') ;
    wfii_linear(:,i) = interp1(pfi(iok,i),wfi(iok,i),pgrid,'linear','extrap') ;    
    % mask extrapolated data below the deepest measurement (== bottom)
    inodata = find(pgrid>maxdepth);
    ufii_linear(inodata,i) = nan;
    vfii_linear(inodata,i) = nan;    
    wfii_linear(inodata,i) = nan;    
    end
end

for i=1:length(dnumi) ;
    iok = find(~isnan(ufi(:,i)));
    if length(iok)>2    
    ufii_akima(:,i) = akima(pfi(iok,i),ufi(iok,i),pgrid) ;
    vfii_akima(:,i) = akima(pfi(iok,i),vfi(iok,i),pgrid) ;
    wfii_akima(:,i) = akima(pfi(iok,i),wfi(iok,i),pgrid) ;
    % mask extrapolated data below the deepest measurement (== bottom)
    inodata = find(pgrid>maxdepth);
    ufii_akima(inodata,i) = nan;
    vfii_akima(inodata,i) = nan;    
    wfii_akima(inodata,i) = nan;    
    end
end

%     
% %% Make a mean velocity from all the methods
% for i= 1:size(ufii_pchip,1);
%     for j=1:size(ufii_pchip,2);
%         ufii(i,j) = nanmean([ufii_pchip(i,j) ufii_linear(i,j) ufii_akima(i,j)]);
%         vfii(i,j) = nanmean([vfii_pchip(i,j) vfii_linear(i,j) vfii_akima(i,j)]);
%         wfii(i,j) = nanmean([wfii_pchip(i,j) wfii_linear(i,j) wfii_akima(i,j)]);       
%     end
% end

%% Make a mean velocity from all the methods
ufii = ufii_linear;
vfii = vfii_linear;
wfii = wfii_linear; 


%% This just saves all the original data. Most of if is crap from the last bin loaded
eval(['save ' out_path moor '_velocity_grid.mat']);


