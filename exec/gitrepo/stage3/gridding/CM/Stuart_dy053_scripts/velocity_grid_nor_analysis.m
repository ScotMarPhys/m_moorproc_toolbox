% velocity_grid_nor_analysis.m


clear all ; close all ;
prognam = 'velocity_grid_nor_analysis' ;

moor = 'rteb1_01_2014' ;
% moor = 'rtwb1_01_2014' ;

% --- get moring information from infofile
infofile =[moor '/' moor 'info.dat'];

%% This just saves all the original data. Most of if is crap from the last bin loaded
eval(['load(''',moor,'_velocity_grid.mat'')']);

%%
if do_plot ;
    figure(10);clf; figure(11);clf;
    for i=1:length(sn) ;
        % Extract the velocity at thepressure level nearest the mean instrument
        % pressure
        [~,Ip] = min(abs(pgrid - nanmean(pfi(i,:))));
        disp(['Inst mean press : ',sprintf('%5.1f',nanmean(pfi(i,:))),' pgrid pressure match at : ',num2str(pgrid(Ip))]);
        
        % timeseries at equivalent pressures of splined data and input data
        figure(10);
        subplot(length(sn),1,i);hold on;grid on;set(gca,'layer','top');
        plot(dnumi,ufi(i,:),'k','linewidth',4);
        plot(dnumi,ufii_pchip(Ip,:),'r','linewidth',2);
        datetick('x',24)
        ylim([-20 20]);
        ylabel('u [cm/s]');
        title(['s/n : ',sprintf('%5d',sn(proc)),' : Depth : ',sprintf('%5.1f',nanmean(p)), 'dbar']);
        
        figure(11);
        subplot(length(sn),1,i);hold on;grid on;set(gca,'layer','top');
        plot(dnumi,vfi(i,:),'k','linewidth',4);
        plot(dnumi,vfii_pchip(Ip,:),'r','linewidth',2);
        
        datetick('x',24)
        ylim([-20 20]);
        ylabel('v [cm/s]');
        title(['s/n : ',sprintf('%5d',sn(proc)),' : Depth : ',sprintf('%5.1f',nanmean(p)), 'dbar']);
        
        %     linkaxes([p1 p2],'x');
        
    end
end
%% Now what do vertical profiles look like
jt=1:60:length(ufii_akima);
figure(12);clf;hold on ;grid on;set(gca,'layer','top');
plot(ufii_pchip(:,jt),pgrid,'r','linewidth',2);
% plot(ufii_spline(:,jt),pgrid,'g','linewidth',2);
plot(ufii_linear(:,jt),pgrid,'b','linewidth',2);
plot(ufii_akima(:,jt),pgrid,'c','linewidth',2);

plot(ufi(:,jt),pfi(:,jt),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);

title(['pchip (r)  ','linear (b)  ','akima (c)  ','data (k)']);
set(gca,'YDIR','reverse');

eval(['print -djpeg f12_vert_prof_',prognam,'_',moor,'_.jpg']);

%% 
figure(13);clf;hold on ; grid on ; set(gca,'layer','top');
ci=-50:10:50;
colormap(lbmap(length(ci)-1,'RedBlue')); 
colormap(flipud(colormap));
[c,h] = contourf(dnumi,pgrid,ufii,ci);
set(get(h,'Children'),'LineStyle','none');
contour(dnumi,pgrid,ufii_pchip,[0 0],'k','linewidth',1)
%     axis([350.4,350.8,0,600]);
ylim([0 1800]);
caxis([min(ci) max(ci)]);
set(gca,'YDIR','reverse');
datetick('x',24)
colorbarf(c,h);
xlabel('date');ylabel('press [db]');title([moor,' : u [cm/s]']);

eval(['print -djpeg f13_u_',prognam,'_',moor,'.jpg']);

figure(14);clf;hold on ; grid on ; set(gca,'layer','top');
ci=-50:10:50;
colormap(lbmap(length(ci)-1,'RedBlue')); 
colormap(flipud(colormap));
[c,h] = contourf(dnumi,pgrid,vfii,ci);
set(get(h,'Children'),'LineStyle','none');
contour(dnumi,pgrid,vfii_pchip,[0 0],'k','linewidth',1)
%     axis([350.4,350.8,0,600]);
ylim([0 1800]);
caxis([min(ci) max(ci)]);
set(gca,'YDIR','reverse');
datetick('x',24)
colorbarf(c,h);
xlabel('date');ylabel('press [db]');title([moor,' : v [cm/s]']);
eval(['print -djpeg f14_v_',prognam,'_',moor,'.jpg']);

%% Lets see a vertial mean profile
figure(15);clf;hold on ; grid on;
plot(nanmean(ufii,2),pgrid,'k','linewidth',2);
plot(nanmean(vfii,2),pgrid,'r','linewidth',2);

plot(nanmean(ufi,2),nanmean(pfi,2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
plot(nanmean(vfi,2),nanmean(pfi,2),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8);
xlabel('velocity [cm/s]');ylabel('press [db]');title([moor]);

set(gca,'YDIR','reverse')
legend('u','v','location','southeast')
eval(['print -djpeg f15_uv_mean_profile_',prognam,'_',moor,'.jpg']);

