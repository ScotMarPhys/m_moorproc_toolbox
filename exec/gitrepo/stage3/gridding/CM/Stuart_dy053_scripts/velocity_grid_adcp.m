% velocity_grid_adcp.m
% Interpolate current meter data (u,v,press)
%  1. 2-day low-pass ; 2. intepolate onto a 12-hour grid ; 3. interpolate
%  vertically using a variety of methods ; 4. create mean u,v

clear all ; close all ;
prognam = 'velocity_grid_adcp' ;

moor = 'rtadcp1_01_2014' ;

cd '/Users/sa02sc/Documents/cruises/DY053/DY053_Cruise_Report/Stuart/RT_Fluxes/adcp/'

% --- get moring information from infofile
infofile =[moor '/' moor 'info.dat'];


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
    for proc=1:32 ;
    if (proc<10);
        infile  = [moor,'/adp/',moor,'_',num2str(sn),'_bin0',sprintf('%d',proc),'.edt'];
    else
        infile  = [moor,'/adp/',moor,'_',num2str(sn),'_bin',sprintf('%d',proc),'.edt'];
    end
    
    [YY,MM,DD,HH,z,t,u,v,w,heading,pitch,roll,spd,direction,Amp1,Amp2,Amp3,Amp4,...
        Beam1Cor,Beam2Cor,Beam3Cor,Beam4Cor,err,PG1,PG2,PG3,PG4] = ...
        rodbload(infile,[columns]);
    dnum = datenum(YY,MM,DD,HH,0,0);
    p=z;
    
    % Blitz the NaNs
    warning('off')
    u = interp1(dnum,u,dnum,'spline');
    v = interp1(dnum,v,dnum,'spline');
    p = interp1(dnum,p,dnum,'spline');
    
    % filtering parameter
    sr  = median(diff(dnum)); % sampling interval
    co = sr/((1/sr)*2); % 2-day low-pass
    
    % 2-day low pass
    uf  = auto_filt(u,sr,co,'low',4);
    vf  = auto_filt(v,sr,co,'low',4);
    pf  = auto_filt(p,sr,co,'low',4);
    
    % put data onto a 12-hour grid
%     dnumi = 735799:0.5:736135; % This is for nortec times
      dnumi = 735903:0.5:736134;
    ufi(proc,:) = interp1(dnum,uf,dnumi,'linear');
    vfi(proc,:) = interp1(dnum,vf,dnumi,'linear');
    pfi(proc,:) = interp1(dnum,pf,dnumi,'linear');
    
    do_plot = 0 ;
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


