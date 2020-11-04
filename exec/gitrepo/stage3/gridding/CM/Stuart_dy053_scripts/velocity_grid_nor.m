% velocity_grid_nor.m
% Interpolate current meter data (u,v,press)
%  1. 2-day low-pass ; 2. intepolate onto a 12-hour grid ; 3. interpolate
%  vertically using a variety of methods ; 4. create mean u,v

clear all ; close all ;
prognam = 'velocity_grid_nor' ;

% moor = 'rteb1_01_2014' ;
moor = 'rtwb1_01_2014' ;

% --- get moring information from infofile
infofile =[moor '/' moor 'info.dat'];


[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

vec=find((id==368|id==370)); % Nortek id number

sn=sn(vec); z = z(vec);

do_plot = 1 ;
for proc = 1 : length(vec);
%         for proc = 1 : 1 ;
    columns = 'YY:MM:DD:HH:T:P:U:V:W:HDG:PIT:ROL:USS:VSS:WSS:IPOW:CS:CD';
    indep  = z(proc);
    infile  = [moor '/nor/',moor,'_',num2str(sn(proc)),'.use'];

    [YY,MM,DD,HH,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = ...
        rodbload(infile,[columns]);
    dnum = datenum(YY,MM,DD,HH,0,0);
    %          disp(['s/n : ',sprintf('%5d',sn(proc)),' : start :,',num2str(dnum(1)),' : end :',num2str(dnum(end))]);
    %     end
    warning('off')
    % Blitz the NaNs
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
    dnumi = 735799:0.5:736135;
    ufi(proc,:) = interp1(dnum,uf,dnumi,'linear');
    vfi(proc,:) = interp1(dnum,vf,dnumi,'linear');
    pfi(proc,:) = interp1(dnum,pf,dnumi,'linear');
    
 if do_plot ;   
% Plot individual records
    figure(proc);clf;
    [p1]=subplot(3,1,1);hold on;grid on;set(gca,'layer','top');
    plot(dnum,u,'k');
    plot(dnum,uf,'y','linewidth',2);
    plot(dnumi,ufi(proc,:),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);
    datetick('x',24)
    ylim([-50 50]);
    ylabel('u [cm/s]');
    title(['s/n : ',sprintf('%5d',sn(proc)),' : Depth : ',sprintf('%5.1f',nanmean(p)), 'dbar']);

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
    ylim([min(p) max(p)]);
    ylabel('press [dbar]');

    linkaxes([p1 p2 p3],'x');
 end

end

%% Now onto vertical interpolation of the subsampled, low-pass data 
% matices: dnumi, ufi(6,length(dnumi)), vfi, pfi

% Define pressure grid
if(moor == 'rteb1_01_2014')
    pgrid = 0:10:1800 ;
else
    pgrid=0:10:1600 ;
end

% interpolation
for i=1:length(dnumi) ;
    ufii_pchip(:,i) = interp1(pfi(:,i),ufi(:,i),pgrid,'pchip','extrap') ;
    vfii_pchip(:,i) = interp1(pfi(:,i),vfi(:,i),pgrid,'pchip','extrap') ;
end

% This is a really bad method
% for i=1:length(dnumi) ;
%     ufii_spline(:,i) = interp1(pfi(:,i),ufi(:,i),pgrid,'spline','extrap') ;
%     vfii_spline(:,i) = interp1(pfi(:,i),vfi(:,i),pgrid,'spline','extrap') ;
% end

for i=1:length(dnumi) ;
    ufii_linear(:,i) = interp1(pfi(:,i),ufi(:,i),pgrid,'linear','extrap') ;
    vfii_linear(:,i) = interp1(pfi(:,i),vfi(:,i),pgrid,'linear','extrap') ;
end

for i=1:length(dnumi) ;
    ufii_akima(:,i) = akima(pfi(:,i),ufi(:,i),pgrid) ;
    vfii_akima(:,i) = akima(pfi(:,i),vfi(:,i),pgrid) ;
end

%% Make a mean velocity from all the methods
for i= 1:size(ufii_pchip,1);
    for j=1:size(ufii_pchip,2);
        ufii(i,j) = nanmean([ufii_pchip(i,j) ufii_linear(i,j) ufii_akima(i,j)]);
        vfii(i,j) = nanmean([vfii_pchip(i,j) vfii_linear(i,j) vfii_akima(i,j)]);
    end
end

%% This just saves all the original data. Most of if is crap from the last bin loaded
eval(['save(''',moor,'_velocity_grid.mat'')']);

