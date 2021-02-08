% velocity_grid_nor_analysis.m

clear all ; close all ;
prognam = 'velocity_grid_nor' ;


for iyear=1:5
switch(iyear)
    case(1)
        year = '_01_2014';
    case(2)
        year = '_02_2015';    
    case(3)
        year = '_03_2016';   
    case(4)
        year = '_04_2017'; 
    case(5)
        year = '_05_2018';
end

 for imoor=1:3
     switch(imoor)
         case(1)
             moorselect = 'rtwb1' ;
         case(2)
             moorselect = 'rtwb2' ;
         case(3)
             moorselect = 'rteb1' ;      
     end
    

moor = [moorselect year];

basedir      = '/media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap';
basedir      = getenv('OSNAP');
gridpath  = [basedir '/data/moor/proc/velocity_grid/'];
plotpath  = [gridpath 'PLOTS/'];

%% This just saves all the original data. Most of if is crap from the last bin loaded
% eval(['load(''',gridpath,moor,'_velocity_grid.mat'')']);
load([gridpath moor '_velocity_grid.mat'])
%%
if do_plot ;
    figure(10);clf; figure(11);clf;
    for i=1:length(pfi(:,1)) ;
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

% eval(['print -djpeg ' plotpath 'f12_vert_prof_',prognam,'_',moor,'_.jpg']);
print(figure(12),'-djpeg', [plotpath 'f12_vert_prof_',prognam,'_',moor,'_.jpg']);

%% 
figure(13);clf;hold on ; grid on ; set(gca,'layer','top');
ci=-50:10:50;
colormap(lbmap(length(ci)-1,'RedBlue')); 
colormap(flipud(colormap));
imagesc(dnumi,pgrid,ufii)
[c,h] = contourf(dnumi,pgrid,ufii,ci);
set(get(h,'Children'),'LineStyle','none');
contour(dnumi,pgrid,ufii_pchip,[0 0],'k','linewidth',1)
%     axis([350.4,350.8,0,600]);
ylim([0 2000]);
caxis([min(ci) max(ci)]);
set(gca,'YDIR','reverse');
datetick('x',24)
colorbar
xlabel('date');ylabel('press [db]');
title([moor(1:4) ' ' moor(7:8) ' ' moor(10:end) ' : u [cm/s]']);
print(figure(13),'-djpeg', [plotpath 'f13_u_',prognam,'_',moor,'.jpg']);

figure(14);clf;hold on ; grid on ; set(gca,'layer','top');
ci=-50:10:50;
colormap(lbmap(length(ci)-1,'RedBlue')); 
colormap(flipud(colormap));
[c,h] = contourf(dnumi,pgrid,vfii,ci);
set(get(h,'Children'),'LineStyle','none');
contour(dnumi,pgrid,vfii_pchip,[0 0],'k','linewidth',1)
%     axis([350.4,350.8,0,600]);
ylim([0 2000]);
caxis([min(ci) max(ci)]);
set(gca,'YDIR','reverse');
datetick('x',24)
colorbar
xlabel('date');ylabel('press [db]');
title([moor(1:4) ' ' moor(7:8) ' ' moor(10:end) ' : v [cm/s]']);
print(figure(14),'-djpeg',[plotpath 'f14_v_',prognam,'_',moor,'.jpg']);

%% Lets see a vertial mean profile
figure(15);clf;hold on ; grid on;
plot(nanmean(ufii,2),pgrid,'k','linewidth',2);
plot(nanmean(vfii,2),pgrid,'r','linewidth',2);

plot(nanmean(ufi,2),nanmean(pfi,2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
plot(nanmean(vfi,2),nanmean(pfi,2),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8);
xlabel('velocity [cm/s]');ylabel('press [db]');

title([moor(1:4) ' ' moor(7:8) ' ' moor(10:end)]);

set(gca,'YDIR','reverse')
legend('u','v','location','southeast')
% eval(['print -djpeg ' plotpath 'f15_uv_mean_profile_',prognam,'_',moor,'.jpg']);
print(figure(15),'-djpeg',[plotpath 'f15_uv_mean_profile_',prognam,'_',moor,'.jpg']);

close all
 end
end
