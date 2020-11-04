% velocity_grid_nor.m
% Interpolate current meter data (u,v,press)
%  1. 2-day low-pass ; 2. intepolate onto a 12-hour grid ; 3. interpolate
%  vertically using a variety of methods ; 4. create mean u,v

clear all ; close all ;
prognam = 'nor_stats.m' ;

% moor = 'rteb1_01_2014' ;
% moor = 'rtwb1_01_2014' ;
moor = 'rtwb2_01_2014' ;

% --- get moring information from infofile
infofile =[moor '/' moor 'info.dat'];


[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

vec=find((id==368|id==370)); % Nortek id number

sn=sn(vec); z = z(vec);


for proc = length(sn) ;
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
    
end
ufi=ufi(proc,:);
vfi=vfi(proc,:);
    

%%

figure(1);clf;
subplot(2,1,1),hold on ; grid on;
plot(dnum,u,'k');
plot(dnumi,ufi,'y','linewidth',2)
datetick('x',24)
title(['s\n : ',sprintf('%5.0f',sn(proc)),' at ',sprintf('%4.0f',z(proc)),' m'])
ylabel('u:zonal');
ylim([-20 20]);

subplot(2,1,2),hold on ; grid on;
plot(dnum,v,'k');
plot(dnumi,vfi,'y','linewidth',2)
datetick('x',24)
title([sprintf('%6.3f',lat),'°N           ',sprintf('%6.3f',lon),'°W'])
ylabel('v:meridional');
ylim([-20 20]);
eval(['print -djpeg f1_uv_time_',prognam,'_',moor,'_.jpg']);
%%
nlag=40;
[acor,lag]=xcorr(vfi ,vfi,nlag,'coeff');
% integrate from zero to first zero crossing
for i=nlag+1:nlag*2;
    d=acor(i+1)-acor(i);
    if(d>0);
        Ilag=i;
    end
end

    
    
% [~,Ilag] = min(abs(acor - 0));
Ilag = 53;
[~,Ilag0] = min(abs(acor - 1));

area=2*trapz(lag(Ilag0:Ilag),acor(Ilag0:Ilag));

figure(2);clf;hold on ;grid on;
plot(lag,acor,'k')
plot(lag(Ilag),acor(Ilag),'ko');

disp(['Timescale of independence : ',sprintf('%6.1f',area/2),' days']); % this is the number of 12-hour lags
mnu=nanmean(ufi);
sdu=nanstd(ufi);
seu=sdu./(sqrt(length(ufi)/area/2));
disp(['Mean u [cm/s] : ',sprintf('%3.1f',mnu),' : SD : ',sprintf('%3.1f',sdu),' : SE : ',sprintf('%3.1f',seu)])

mnv=nanmean(vfi);
sdv=nanstd(vfi);
sev=sdv./(sqrt(length(vfi)/area/2));
disp(['Mean v [cm/s] : ',sprintf('%3.1f',mnv),' : SD : ',sprintf('%3.1f',sdv),' : SE : ',sprintf('%3.1f',sev)])








