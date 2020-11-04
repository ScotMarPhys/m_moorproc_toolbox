% velocity_grid_adcp_analysis.m
% Read gridded adcp file produced by velocity_grid_adcp.m

clear all ; close all ;
prognam = 'velocity_grid_adcp_analysis' ;

% moor = 'rtadcp1_01_2014' ;

cd '/Users/sa02sc/Documents/cruises/DY053/DY053_Cruise_Report/Stuart/RT_Fluxes/adcp/'

%  Most of if is crap from the last bin loaded
load('rtadcp1_01_2014_velocity_grid.mat')
%%
if do_plot ;
    figure(1);clf; figure(12);clf;
    ic = 1;        nsub = 6;
    for i=1:6:36 ;

        % Extract the velocity at thepressure level nearest the mean instrument
        % pressure
        [~,Ip] = min(abs(pgrid - nanmean(pfi(i,:))));
%         disp(['Inst mean press : ',sprintf('%5.1f',nanmean(pfi(i,:))),' pgrid pressure match at : ',num2str(pgrid(Ip))]);
        
        % timeseries at equivalent pressures of splined data and input data
        figure(10);
        subplot(nsub,1,ic);hold on;grid on;set(gca,'layer','top');
        plot(dnumi,ufi(i,:),'k','linewidth',4);
        plot(dnumi,ufii(Ip,:),'r','linewidth',2);
        datetick('x',24)
        ylim([-25 25]);
        ylabel('u [cm/s]');
        title(['Bin : ',sprintf('%2d',proc),' : Depth : ',sprintf('%5.1f',nanmean(pgrid(i))), 'dbar']);
        
        figure(11);
        subplot(nsub,1,ic);hold on;grid on;set(gca,'layer','top');
        plot(dnumi,vfi(i,:),'k','linewidth',4);
        plot(dnumi,vfii(Ip,:),'r','linewidth',2);
        
        datetick('x',24)
        ylim([-50 50]);
        ylabel('v [cm/s]');
        title(['Bin : ',sprintf('%2d',proc),' : Depth : ',sprintf('%5.1f',nanmean(pgrid(i))), 'dbar']);
        
        ic=ic+1
        
    end
end
%% Now what do vertical profiles look like
jt=1:20:length(ufii);
figure(12);clf;hold on ;grid on;set(gca,'layer','top');
% plot(ufii_pchip(:,jt),pgrid,'r','linewidth',2);
% plot(ufii_spline(:,jt),pgrid,'g','linewidth',2);
plot(ufii_linear(:,jt),pgrid,'b','linewidth',2);
% plot(ufii_akima(:,jt),pgrid,'c','linewidth',2);

plot(ufi(:,jt),pfi(:,jt),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);

% title(['pchip (r)  ','linear (b)  ','akima (c)  ','data (k)']);
set(gca,'YDIR','reverse');
xlim([-30 30]);
eval(['print -djpeg f12_vert_prof_',prognam,'_',moor,'_.jpg']);

%% 
figure(13);clf;hold on ; grid on ; set(gca,'layer','top');
ci=-50:10:50;
colormap(lbmap(length(ci)-1,'RedBlue')); 
colormap(flipud(colormap));
[c,h] = contourf(dnumi,pgrid,ufii,ci);
set(get(h,'Children'),'LineStyle','none');
contour(dnumi,pgrid,ufii,[0 0],'k','linewidth',1)
%     axis([350.4,350.8,0,600]);
ylim([0 750]);
caxis([min(ci) max(ci)]);
set(gca,'YDIR','reverse');
datetick('x',24)
colorbarf(c,h);
xlabel('date');ylabel('press [db]');title([moor,' : zonal [cm/s]']);

eval(['print -djpeg f13_u_',prognam,'_',moor,'.jpg']);

figure(14);clf;hold on ; grid on ; set(gca,'layer','top');
ci=-50:10:50;
colormap(lbmap(length(ci)-1,'RedBlue')); 
colormap(flipud(colormap));
[c,h] = contourf(dnumi,pgrid,vfii,ci);
set(get(h,'Children'),'LineStyle','none');
contour(dnumi,pgrid,vfii,[0 0],'k','linewidth',1)
%     axis([350.4,350.8,0,600]);
ylim([0 750]);
caxis([min(ci) max(ci)]);
set(gca,'YDIR','reverse');
datetick('x',24)
colorbarf(c,h);
xlabel('date');ylabel('press [db]');title([moor,' : meridional [cm/s]']);
eval(['print -djpeg f14_v_',prognam,'_',moor,'.jpg']);

%% Lets see a vertial mean profile

Irot = 20;
Irot=24; % This value from the chart analysis 
[u_across,v_along] = uvrot(ufii,vfii,Irot);


figure(15);clf;hold on ; grid on;
plot(nanmean(ufii,2),pgrid,'k','linewidth',2);
plot(nanmean(vfii,2),pgrid,'r','linewidth',2);

plot(nanmean(u_across,2),pgrid,'k+','linewidth',3);
plot(nanmean(v_along,2),pgrid,'r+','linewidth',3);

% plot(nanmean(ufi,2),nanmean(pfi,2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
% plot(nanmean(vfi,2),nanmean(pfi,2),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8);
xlabel('velocity [cm/s]');ylabel('press [db]');title([moor,' : Irot= ',sprintf('%3.1f',Irot),'°']);

set(gca,'YDIR','reverse')
legend('zonal','meridional','across','along','location','southeast')
eval(['print -djpeg f15_uv_mean_profile_',prognam,'_',moor,'.jpg']);

%%

ufii_mn = nanmean(ufii,2);vfii_mn = nanmean(vfii,2);
%%
% figure(155);clf;hold on ;grid on;
% 
% quiver(0,0,ufii_mn(:),vfii_mn(:))
% % axis('equal')
%% Timeseries

figure(16);clf;hold on;grid on;
pref=260;[~,Ip] = min(abs(pgrid - pref));plot(dnumi,vfii(Ip,:),'k','linewidth',2);
pref=300;[~,Ip] = min(abs(pgrid - pref));plot(dnumi,vfii(Ip,:),'r','linewidth',2);
pref=500;[~,Ip] = min(abs(pgrid - pref));plot(dnumi,vfii(Ip,:),'g','linewidth',2);
pref=620;[~,Ip] = min(abs(pgrid - pref));plot(dnumi,vfii(Ip,:),'b','linewidth',2);
plot([dnumi(1) dnumi(end)],[0 0],'k','linewidth',3);

datetick('x',24);
legend('260','300','500','620');
xlabel('date');ylabel('meridional velocity [cm/s]');title('RTADCP1');
% title([sprintf('%5.1f',pgrid(Ip)),' dbar'])

eval(['print -djpeg f16_v_timeseries_',prognam,'_',moor,'.jpg']);

%% find rotation angle to minimise u_across

pref=620;[~,Ip] = min(abs(pgrid - pref));
u_zonal = ufii(Ip,:); 
v_merid = vfii(Ip,:);

% Irot = 23.9 ; % This angle by geometry on ETOPO2 chart by hand
% [u_across,v_along] = uvrot(u_zonal,v_merid,Irot);
% disp(['Rotation : ',sprintf('%4.1f',Irot),'°']);
% sprintf(['Mean across slope current ',sprintf('%3.1f',nanmean(u_across)),' cm/s'])



disp('Find angle of rotation that minimises across slope flow');
disp('angle across along');

izcount=1;
for iz = 1:length(pgrid);
    
icount=1;
Irot=16:0.1:21;
for i=1:length(Irot);
clear u_across v_along ;
    [u_across,v_along] = uvrot(ufii(iz,:),vfii(iz,:),Irot(i));
    u_across_mn(icount) = nanmean(u_across) ;
% disp([sprintf('%4.1f',Irot),' : ',sprintf('%4.2f',nanmean(u_across))]);
icount = icount+1 ;
end

    [~,Irot_minimum] = min(abs(u_across_mn));
    Irotz(izcount)=Irot(Irot_minimum);
    izcount = izcount + 1 ;
end


figure(17);clf;hold on;grid on;
plot(Irotz,pgrid,'ko-','linewidth',2)
set(gca,'YDIR','reverse');
ylabel('depth');xlabel('angle');title('angle to minimise mean across velocity');
eval(['print -djpeg f17_rotation_angle_',prognam,'_',moor,'.jpg']);
%% 18.4° angle minimises across slope flow at the top of the boundary layer

pref=380;[~,Ip] = min(abs(pgrid - pref));
u_zonal = ufii(Ip,:); 
v_merid = vfii(Ip,:);

Irot = 18.4 ; 
[u_across,v_along] = uvrot(u_zonal,v_merid,Irot);

disp(['Mean zonal flow : ',sprintf('%4.2f',nanmean(u_zonal)),' cm/s']);
disp(['Rotate by : ',sprintf('%3.1f',Irot),'°']);
disp(['Mean across flow : ',sprintf('%4.2f',nanmean(u_across)),' cm/s']);

figure(18);clf
subplot(2,1,1);hold on ;grid on ;
plot(dnumi,u_zonal,'r','linewidth',2)
plot(dnumi,v_merid,'k','linewidth',2)
datetick('x')
ylim([-20 50]);
ylabel('velocity');xlabel('date');title('depth = 380m');
legend('zonal','meridional');

subplot(2,1,2);hold on ;grid on ;
plot(dnumi,u_across,'r','linewidth',2)
plot(dnumi,v_along,'k','linewidth',2)
datetick('x')
ylim([-20 50]);
ylabel('velocity');xlabel('date');title('Rotated by 18.4°');
legend('u\_across','v\_along');

eval(['print -djpeg f18_rotated_velocities_',prognam,'_',moor,'.jpg']);

