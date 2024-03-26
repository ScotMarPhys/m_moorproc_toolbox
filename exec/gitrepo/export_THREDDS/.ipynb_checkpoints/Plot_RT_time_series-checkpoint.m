% PLOT_RT_TIME_SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars; close('all');
indir   = 'C:\Users\sa01ld\m_moorproc_toolbox';

datafile = fullfile(indir,'data/moor/THREDDS/Rockall-Trough-Mooring-Time-Series-2022.nc');
finfo=ncinfo(datafile);

mooring='EAST';
% mooring='WEST';

% index of west mooring for current meter plots, either 1 or 2
 w_nm='2';

% cruise dates
KN221=datenum(2014,7,18,00,00,00);
PE399=datenum(2015,6,21,00,00,00);
DY053=datenum(2016,7,01,00,00,00);
dy078=datenum(2017,5,12,00,00,00);
ar30=datenum(2018,7,10,00,00,00);
dy120=datenum(2020,10,18,00,00,00);
jc238=datenum(2022,07,25,00,00,00);


%% 1. produce merged plot with cruise occupations
% figure;
% %%%%%%%%%%%%%%%%%% TEMPERATURE %%%%%%%%%%%%%%%%%%%
% x=ncread(datafile,'TIME')+datenum(1950,1,1,0,0,0);
% y=ncread(datafile,'PRES');
% z=gsw_pt_from_CT(ncread(datafile,['SG_' mooring]),ncread(datafile,['TG_' mooring])-273.15);
% ax(1)=subplot(3,1,1)
% [c,h]=contourf(x,y,z,100,'LineColor','none');
% caxis([ceil(min(min(z))) floor(max(max(z)))]);
% C = colorbar;
% C.Label.String = '\theta (^{o} C)';
% C.Ticks=[floor(min(min(z))):2:ceil(max(max(z)))];
% C.TickLength=0.135;
% ylabel(gca,'Depth (m)');
% [Y,~,~]=datevec(x);
% date1=datenum(min(Y),1,1);
% date2=datenum(max(Y),1,1);
% datetick('x','Keeplimits');
% axis ij
% % title('Rockall Trough Eastern Boundary Gridded Temperature')
% cmap=cmocean('Thermal');
% colormap(cmap);
% text(KN221,0,'KN221','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
% text(PE399,0,'PE399','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
% text(DY053,0,'DY053','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
% text(dy078,0,'DY078','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
% text(dy120,0,'DY120','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
% text(ar30,0,'AR30','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
% hold on
% plot([KN221 KN221],[0 1800],'--k','LineWidth',2);
% plot([PE399 PE399],[0 1800],'--k','LineWidth',2);
% plot([DY053 DY053],[0 1800],'--k','LineWidth',2);
% plot([dy078 dy078],[0 1800],'--k','LineWidth',2);
% plot([dy120 dy120],[0 1800],'--k','LineWidth',2);
% plot([ar30 ar30],[0 1800],'--k','LineWidth',2);
% 
% %%%%%%%%%%%%%%%%%%% SALINITY %%%%%%%%%%%%%%%%%%%
% ax(2)=subplot(3,1,2);
% x=ncread(datafile,'TIME')+datenum(1950,1,1,0,0,0);
% y=ncread(datafile,'PRES');
% z=ncread(datafile,['SG_' mooring]);
% [c,h]=contourf(x,y,z,12,'LineColor','none');
% caxis([round(min(min(z)),1) round(max(max(z)),1)]);
% C = colorbar;
% cmap=cmocean('-Haline');
% colormap(ax(2),cmap);
% C.Label.String = 'S_{a} (g kg^{-1})';
% C.Ticks=[round(min(min(z)),1):0.2:round(max(max(z)),1)];
% C.TickLength=0.135;
% hold on
% % [c,h]=contour(JG , pgg, SGfs,[35.2 35.4 35.6],'LineColor','k','LineWidth',0.25); 
% % clabel(c,h,'LabelSpacing',500,'Color','w')
% ylabel(gca,'Depth (m)');
% % SORT LABELS
% [Y,~,~]=datevec(x);
% date1=datenum(min(Y),1,1);
% date2=datenum(max(Y),1,1);
% datetick('x','Keeplimits');
% axis ij
% hold on
% plot([KN221 KN221],[0 1800],'--k','LineWidth',2);
% plot([PE399 PE399],[0 1800],'--k','LineWidth',2);
% plot([DY053 DY053],[0 1800],'--k','LineWidth',2);
% plot([dy078 dy078],[0 1800],'--k','LineWidth',2);
% plot([dy120 dy120],[0 1800],'--k','LineWidth',2);
% plot([ar30 ar30],[0 1800],'--k','LineWidth',2);
% 
% %%%%%%%%%%%%%%%%%%% DENSITY %%%%%%%%%%%%%%%%%%%
% ax(3)=subplot(3,1,3);
% x=ncread(datafile,'TIME')+datenum(1950,1,1,0,0,0);
% y=ncread(datafile,'PRES');
% z=gsw_rho(ncread(datafile,['SG_' mooring]),ncread(datafile,['TG_' mooring])-273.15,0)-1000;
% [c,h]=contourf(x,y,z,12,'LineColor','none');
% caxis([round(min(min(z)),1) round(max(max(z)),1)]);
% C = colorbar;
% cmap=cmocean('Dense');
% colormap(ax(3),cmap);
% C.Label.String = '\sigma (kg/m^{3})';
% C.Ticks=[round(min(min(z)),1):0.4:round(max(max(z)),1)];
% C.TickLength=0.135;
% hold on
% % [c,h]=contour(JG , pgg, SGfs,[35.2 35.4 35.6],'LineColor','k','LineWidth',0.25); 
% % clabel(c,h,'LabelSpacing',500,'Color','w')
% ylabel(gca,'Depth (m)');
% % SORT LABELS
% [Y,~,~]=datevec(x);
% date1=datenum(min(Y),1,1);
% date2=datenum(max(Y),1,1);
% datetick('x','Keeplimits');
% axis ij
% hold on
% plot([KN221 KN221],[0 1800],'--k','LineWidth',2);
% plot([PE399 PE399],[0 1800],'--k','LineWidth',2);
% plot([DY053 DY053],[0 1800],'--k','LineWidth',2);
% plot([dy078 dy078],[0 1800],'--k','LineWidth',2);
% plot([dy120 dy120],[0 1800],'--k','LineWidth',2);
% plot([ar30 ar30],[0 1800],'--k','LineWidth',2);
% 
% 
% width=35; height=20; FS=20; FN='Helvetica';
% set(ax(1:3),'fontsize', FS, 'FontName',FN);
% set(gcf,'units','centimeters','position',[5 5 width height])
% 
% print('-dpng',['Figures/' mooring '_TSRho'])
% print('-dpng','C:\Users\sa01ld\Desktop\OSNAP\ScotMarPhys.OSNAP-Mooring-Processing.io\img\rtwb_all')

%% 2.plot gridded T-S data
figure;
KN221=datenum(2014,7,18,00,00,00);
PE399=datenum(2015,6,21,00,00,00);
DY053=datenum(2016,7,01,00,00,00);
dy078=datenum(2017,5,12,00,00,00);
ar30=datenum(2018,7,10,00,00,00);
dy120=datenum(2020,10,18,00,00,00);

%%%%%%%%%%%%%%%%%%% TEMPERATURE %%%%%%%%%%%%%%%%%%%
x=ncread(datafile,'TIME')+datenum(1950,1,1,0,0,0);
y=ncread(datafile,'PRES');
z=gsw_pt_from_CT(ncread(datafile,['SG_' mooring]),ncread(datafile,['TG_' mooring])-273.15);
ax(1)=subplot(3,1,1);
[c,h]=contourf(x,y,z,100,'LineColor','none');
caxis([ceil(min(min(z))) floor(max(max(z)))]);
C = colorbar;
C.Label.String = '\theta (^{o} C)';
C.Ticks=[floor(min(min(z))):2:ceil(max(max(z)))];
C.TickLength=0.135;
ylabel(gca,'Depth (m)');
[Y,~,~]=datevec(x);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
axis ij
% title('Rockall Trough Eastern Boundary Gridded Temperature')
cmap=cmocean('Thermal');
colormap(cmap);
text(KN221,0,'KN221','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(PE399,0,'PE399','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(DY053,0,'DY053','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(dy078,0,'DY078','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(dy120,0,'DY120','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(ar30,0,'AR30','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(jc238,0,'JC238','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
hold on
plot([KN221 KN221],[0 1800],'--k','LineWidth',2);
plot([PE399 PE399],[0 1800],'--k','LineWidth',2);
plot([DY053 DY053],[0 1800],'--k','LineWidth',2);
plot([dy078 dy078],[0 1800],'--k','LineWidth',2);
plot([dy120 dy120],[0 1800],'--k','LineWidth',2);
plot([ar30 ar30],[0 1800],'--k','LineWidth',2);
plot([jc238 jc238],[0 1800],'--k','LineWidth',2);


%%%%%%%%%%%%%%%%%%% SALINITY %%%%%%%%%%%%%%%%%%%
ax(2)=subplot(3,1,2);
x=ncread(datafile,'TIME')+datenum(1950,1,1,0,0,0);
y=ncread(datafile,'PRES');
z=ncread(datafile,['SG_' mooring]);
[c,h]=contourf(x,y,z,12,'LineColor','none');
caxis([round(min(min(z)),1) round(max(max(z)),1)]);
C = colorbar;
cmap=cmocean('-Haline');
colormap(ax(2),cmap);
C.Label.String = 'S_{a} (g kg^{-1})';
C.Ticks=[round(min(min(z)),1):0.2:round(max(max(z)),1)];
C.TickLength=0.135;
hold on
% [c,h]=contour(JG , pgg, SGfs,[35.2 35.4 35.6],'LineColor','k','LineWidth',0.25); 
% clabel(c,h,'LabelSpacing',500,'Color','w')
ylabel(gca,'Depth (m)');
% SORT LABELS
[Y,~,~]=datevec(x);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
axis ij
hold on
plot([KN221 KN221],[0 1800],'--k','LineWidth',2);
plot([PE399 PE399],[0 1800],'--k','LineWidth',2);
plot([DY053 DY053],[0 1800],'--k','LineWidth',2);
plot([dy078 dy078],[0 1800],'--k','LineWidth',2);
plot([dy120 dy120],[0 1800],'--k','LineWidth',2);
plot([ar30 ar30],[0 1800],'--k','LineWidth',2);


%%%%%%%%%%%%%%%%%%% DENSITY %%%%%%%%%%%%%%%%%%%
ax(3)=subplot(3,1,3);
x=ncread(datafile,'TIME')+datenum(1950,1,1,0,0,0);
y=ncread(datafile,'PRES');
z=gsw_rho(ncread(datafile,['SG_' mooring]),ncread(datafile,['TG_' mooring])-273.15,0)-1000;
[c,h]=contourf(x,y,z,12,'LineColor','none');
caxis([round(min(min(z)),1) round(max(max(z)),1)]);
C = colorbar;
cmap=cmocean('Dense');
colormap(ax(3),cmap);
C.Label.String = '\sigma (kg/m^{3})';
C.Ticks=[round(min(min(z)),1):0.4:round(max(max(z)),1)];
C.TickLength=0.135;
hold on
% [c,h]=contour(JG , pgg, SGfs,[35.2 35.4 35.6],'LineColor','k','LineWidth',0.25); 
% clabel(c,h,'LabelSpacing',500,'Color','w')
ylabel(gca,'Depth (m)');
% SORT LABELS
[Y,~,~]=datevec(x);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
axis ij
hold on
plot([KN221 KN221],[0 1800],'--k','LineWidth',2);
plot([PE399 PE399],[0 1800],'--k','LineWidth',2);
plot([DY053 DY053],[0 1800],'--k','LineWidth',2);
plot([dy078 dy078],[0 1800],'--k','LineWidth',2);
plot([dy120 dy120],[0 1800],'--k','LineWidth',2);
plot([ar30 ar30],[0 1800],'--k','LineWidth',2);


width=35; height=20; FS=20; FN='Helvetica';
set(ax(1:3),'fontsize', FS, 'FontName',FN);
set(gcf,'units','centimeters','position',[5 5 width height])
print('-dpng',['C:\Users\sa01ld\Desktop\OSNAP\ScotMarPhys.OSNAP-Mooring-Processing.io\img\' mooring '_TS'])


%% 3 plot currents
if strcmpi(mooring,'WEST')
    mooring=[mooring '_' w_nm];
end

figure;
%%%%%%%%%%%%%%%%%%% U %%%%%%%%%%%%%%%%%%%
x=ncread(datafile,'TIME')+datenum(1950,1,1,0,0,0);
y=ncread(datafile,'PRES');
z=ncread(datafile,['U_' mooring])/100;
figure;
[c,h]=contourf(x,y,z,100,'LineColor','none');
caxis([min(min(z)) max(max(z))]);
C = colorbar;cmap=cmocean('Balance','pivot',0);colormap(cmap);
C.Label.String = 'm s^{-1}';
C.Ticks=[floor(min(min(z))):0.1:ceil(max(max(z)))];
C.TickLength=0.06;
ylabel(gca,'Depth (m)');
[Y,~,~]=datevec(x);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
text(KN221,0,'KN221','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(PE399,0,'PE399','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(DY053,0,'DY053','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(dy078,0,'DY078','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(dy120,0,'DY120','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
text(ar30,0,'AR30','Rotation',45,'fontsize', 16, 'FontName','Helvetica');
hold on
plot([KN221 KN221],[0 1800],'--k','LineWidth',2);
plot([PE399 PE399],[0 1800],'--k','LineWidth',2);
plot([DY053 DY053],[0 1800],'--k','LineWidth',2);
plot([dy078 dy078],[0 1800],'--k','LineWidth',2);
plot([dy120 dy120],[0 1800],'--k','LineWidth',2);
plot([ar30 ar30],[0 1800],'--k','LineWidth',2);
axis ij
width=35; height=10; FS=14; FN='Helvetica';
set(gca,'fontsize', FS, 'FontName',FN);
set(gcf,'units','centimeters','position',[5 5 width height])
print('-dpng',['Figures/' mooring '_NOR_V'])
print('-dpng',['C:\Users\sa01ld\Desktop\OSNAP\ScotMarPhys.OSNAP-Mooring-Processing.io\img\' mooring '_NOR_V']);

%%%%%%%%%%%%%%%%%%% V %%%%%%%%%%%%%%%%%%%
x=ncread(datafile,'TIME')+datenum(1950,1,1,0,0,0);
y=ncread(datafile,'PRES');
z=ncread(datafile,['V_' mooring])/100;
figure;[c,h]=contourf(x,y,z,100,'LineColor','none');
caxis([min(min(z)) max(max(z))]);
C = colorbar;cmap=cmocean('Balance','pivot',0);colormap(cmap);
C.Label.String = 'm s^{-1}';
C.Ticks=[floor(min(min(z))):0.1:ceil(max(max(z)))];
C.TickLength=0.06;
ylabel(gca,'Depth (m)');
[Y,~,~]=datevec(x);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
hold on
plot([KN221 KN221],[0 1800],'--k','LineWidth',2);
plot([PE399 PE399],[0 1800],'--k','LineWidth',2);
plot([DY053 DY053],[0 1800],'--k','LineWidth',2);
plot([dy078 dy078],[0 1800],'--k','LineWidth',2);
plot([dy120 dy120],[0 1800],'--k','LineWidth',2);
plot([ar30 ar30],[0 1800],'--k','LineWidth',2);
axis ij
width=35; height=10; FS=14; FN='Helvetica';
set(gca,'fontsize', FS, 'FontName',FN);
set(gcf,'units','centimeters','position',[5 5 width height])
print('-dpng',['Figures/' mooring '_NOR_U'])
print('-dpng',['C:\Users\sa01ld\Desktop\OSNAP\ScotMarPhys.OSNAP-Mooring-Processing.io\img\' mooring '_NOR_U']);

%%%%%%%%%%%%%%%%%%% W %%%%%%%%%%%%%%%%%%%

x=ncread(datafile,'TIME')+datenum(1950,1,1,0,0,0);
y=ncread(datafile,'PRES');
z=ncread(datafile,['W_' mooring])/100;
figure;
[c,h]=contourf(x,y,z,100,'LineColor','none');
caxis([min(min(z)) max(max(z))]);
C = colorbar;cmap=cmocean('Balance','pivot',0);colormap(cmap);
C.Label.String = 'm s^{-1}';
C.Ticks=[floor(min(min(z))):0.01:ceil(max(max(z)))];
C.TickLength=0.06;
ylabel(gca,'Depth (m)');
[Y,~,~]=datevec(x);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
hold on
plot([KN221 KN221],[0 1800],'--k','LineWidth',2);
plot([PE399 PE399],[0 1800],'--k','LineWidth',2);
plot([DY053 DY053],[0 1800],'--k','LineWidth',2);
plot([dy078 dy078],[0 1800],'--k','LineWidth',2);
plot([dy120 dy120],[0 1800],'--k','LineWidth',2);
plot([ar30 ar30],[0 1800],'--k','LineWidth',2);
axis ij
width=35; height=10; FS=14; FN='Helvetica';
set(gca,'fontsize', FS, 'FontName',FN);
set(gcf,'units','centimeters','position',[5 5 width height])
print('-dpng',['Figures/' mooring '_NOR_W'])
print('-dpng',['C:\Users\sa01ld\Desktop\OSNAP\ScotMarPhys.OSNAP-Mooring-Processing.io\img\' mooring '_NOR_W']);

clear