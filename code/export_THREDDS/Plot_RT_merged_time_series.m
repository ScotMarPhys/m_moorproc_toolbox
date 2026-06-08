%% Input data
tsdir                       = [basedir '/osnap/data/moor/proc/hydro_grid_merged/'];
vldir                       = [basedir '/osnap/data/moor/proc/velocity_grid_merged/'];%[pathosnap '/data/moor/proc/velocity_grid_merged/'];
thredds                     = [basedir '/osnap/data/moor/THREDDS/']
%% T S DATA 
% Load western array data
load([tsdir 'RTWB_merg_linear_interp_2022.mat']);
% rename var and make NaN 99999
TG_WEST                 =RTWB_merg.TGfs2;
% TG_WEST(isnan(TG_WEST)) =99999;

SG_WEST                 =RTWB_merg.SGfs2;
pressure                =RTWB_merg.PGfs(:,1); 
time                    =RTWB_merg.JG; 
    
%Load eastern array data
load([tsdir 'RTEB_merg_linear_interp_2022.mat']);
% rename vars
TG_EAST                 =RTEB_merg.TGfs2;
% TG_EAST(isnan(TG_EAST)) =99999;
SG_EAST                 =RTEB_merg.SGfs2; 

%  VELOCITY
% western boundary 1
ffile               ='RTWB1_merg_linear_interp_201407_202407.mat';
load([vldir ffile]);
% rename vars
U_WEST_1             =RTWB1_merg_CM.UGfs2; 
V_WEST_1             =RTWB1_merg_CM.VGfs2; 

% western boundary 2
ffile               ='RTWB2_merg_linear_interp_201407_202407.mat';
load([vldir ffile]);
% rename vars
U_WEST_2             =RTWB2_merg_CM.UGfs2;
V_WEST_2             =RTWB2_merg_CM.VGfs2;
            
% eastern boundary 
ffile               ='RTEB1_merg_linear_interp_201407_202407.mat';
load([vldir ffile]);
% rename vars
U_EAST             =RTEB_merg_CM.UGfs2;
V_EAST             =RTEB_merg_CM.VGfs2;


% MAKE MERGED DATA ONLY FIGURE
% Temperature
% cruise dates
KN221=datenum(2014,7,18,00,00,00);
PE399=datenum(2015,6,21,00,00,00);
DY053=datenum(2016,7,01,00,00,00);
dy078=datenum(2017,5,12,00,00,00);
ar30=datenum(2018,7,10,00,00,00);
dy120=datenum(2020,10,18,00,00,00);
jc238=datenum(2022,07,25,00,00,00);
DY181=datenum(2024,07,12,00,00,00);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rteb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.04, 0.7, 0.9]);

ax(1)=subplot(3,1,1);
%[left bottom width height]
ax(1).Position=[0.05 0.65 0.9 0.28];
x=time;y=pressure;z=TG_EAST;
[c,h]=contourf(x,y,z,36,'LineColor','none');
caxis([ceil(min(min(TG_EAST))) floor(max(max(TG_EAST)))]);
cmap=cmocean('Thermal',45);
colormap(ax(1),cmap(9:end,:));
hold on
[c,h]=contour(x , y, z,[2 4 6 8 10 12 14], '--','LineColor','k');
h.LineWidth = 0.0001;
clabel(c,h,'LabelSpacing',1000,'Color','w')
% xlabel(gca,'Year');
ylabel(gca,'Depth (m)');
[Y,~,~]=datevec(time);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
set(gca,'XTickLabel',[]);
C = colorbar;
C.Label.String = 'Temperature (^{o} C)';
C.Ticks=[floor(min(min(TG_EAST))):1:ceil(max(max(TG_EAST)))];
C.TickLength=0.03;
axis ij

% title('Rockall Trough Eastern Boundary Gridded Temperature')
text(KN221,0,'KN221','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(PE399,0,'PE399','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(DY053,0,'DY053','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(dy078,0,'DY078','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(dy120,0,'DY120','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(ar30,0,'AR30','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(jc238,0,'JC238','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(DY181,0,'DY181','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
hold on
plot([KN221 KN221],[0 1800],'--w','LineWidth',1);
plot([PE399 PE399],[0 1800],'--w','LineWidth',1);
plot([DY053 DY053],[0 1800],'--w','LineWidth',1);
plot([dy078 dy078],[0 1800],'--w','LineWidth',1);
plot([dy120 dy120],[0 1800],'--w','LineWidth',1);
plot([ar30 ar30],[0 1800],'--w','LineWidth',1);
plot([jc238 jc238],[0 1800],'--w','LineWidth',1);
plot([DY181 DY181],[0 1800],'--w','LineWidth',1);

% Salinity

ax(2)=subplot(3,1,2)
%[left bottom width height]
ax(2).Position=[0.05 0.35 0.9 0.28];
x=time;y=pressure;z=SG_EAST;
[c,h]=contourf(x,y,z,12,'LineColor','none');
caxis([round(min(min(z)),1) round(max(max(z)),1)]);
C = colorbar;
cmap=cmocean('Haline');
colormap(ax(2),cmap);
C.Label.String = 'Salinity (g kg^{-1})';
C.Ticks=[round(min(min(z)),1):0.1:round(max(max(z)),1)];
C.TickLength=0.03;
hold on
[c,h]=contour(time , pressure, SG_EAST,[35.2 35.4 35.6],'LineColor','k','LineWidth',0.25); 
clabel(c,h,'LabelSpacing',1000,'Color','w')
% xlabel(gca,'Year');
ylabel(gca,'Depth (m)');
% SORT LABELS
[Y,~,~]=datevec(time);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
set(gca,'XTickLabel',[]);
axis ij
% title('Rockall Trough Eastern Boundary Gridded Salinity')
hold on
plot([KN221 KN221],[0 1800],'--w','LineWidth',1);
plot([PE399 PE399],[0 1800],'--w','LineWidth',1);
plot([DY053 DY053],[0 1800],'--w','LineWidth',1);
plot([dy078 dy078],[0 1800],'--w','LineWidth',1);
plot([dy120 dy120],[0 1800],'--w','LineWidth',1);
plot([ar30 ar30],[0 1800],'--w','LineWidth',1);
plot([jc238 jc238],[0 1800],'--w','LineWidth',1);
plot([DY181 DY181],[0 1800],'--w','LineWidth',1);

% Density

ax(3)=subplot(3,1,3);
%[left bottom width height]
ax(3).Position=[0.05 0.05 0.9 0.28];
x=time;y=pressure;z=gsw_rho(SG_EAST,TG_EAST,0); 
[c,h]=contourf(x,y,z,12,'LineColor','none');
caxis([round(min(min(z)),1) round(max(max(z)),1)]);
C = colorbar;
cmap=cmocean('Dense');
colormap(ax(3),cmap);
C = colorbar;
C.Label.String = 'Potential density (kg/m^{3})';
C.Ticks=[round(min(min(z)),1):0.2:round(max(max(z)),1)];
C.TickLength=0.03;
hold on
[c,h]=contour(time , pressure, z,[1027.8 1027.5 1027.2],'LineColor','k','LineWidth',0.25); 
clabel(c,h,'LabelSpacing',1000,'Color','w')
xlabel(gca,'Year');
ylabel(gca,'Depth (m)');
% SORT LABELS
[Y,~,~]=datevec(time);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
axis ij
% title('Rockall Trough Eastern Boundary Gridded Density')
hold on
plot([KN221 KN221],[0 1800],'--w','LineWidth',1);
plot([PE399 PE399],[0 1800],'--w','LineWidth',1);
plot([DY053 DY053],[0 1800],'--w','LineWidth',1);
plot([dy078 dy078],[0 1800],'--w','LineWidth',1);
plot([dy120 dy120],[0 1800],'--w','LineWidth',1);
plot([ar30 ar30],[0 1800],'--w','LineWidth',1);
plot([jc238 jc238],[0 1800],'--w','LineWidth',1);
plot([DY181 DY181],[0 1800],'--w','LineWidth',1);

% save figure
width=35; height=20; FS=12; FN='Arial';
for hh=1:numel(ax)
    set(ax(hh),'fontsize', FS, 'FontName',FN);
end
print('-dpng',[thredds 'RTEB_MERGED_MCAT'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rtwb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.04, 0.7, 0.9]);

ax(1)=subplot(3,1,1);
%[left bottom width height]
ax(1).Position=[0.05 0.65 0.9 0.28];
x=time;y=pressure;z=TG_WEST;
[c,h]=contourf(x,y,z,36,'LineColor','none');
caxis([ceil(min(min(TG_WEST))) floor(max(max(TG_WEST)))]);
cmap=cmocean('Thermal',45);
colormap(ax(1),cmap(9:end,:));
hold on
[c,h]=contour(x , y, z,[2 4 6 8 10 12 14], '--','LineColor','k');
h.LineWidth = 0.0001;
clabel(c,h,'LabelSpacing',1000,'Color','w')
% xlabel(gca,'Year');
ylabel(gca,'Depth (m)');
[Y,~,~]=datevec(time);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
set(gca,'XTickLabel',[]);
C = colorbar;
C.Label.String = 'Temperature (^{o} C)';
C.Ticks=[floor(min(min(TG_WEST))):1:ceil(max(max(TG_WEST)))];
C.TickLength=0.03;
axis ij

% title('Rockall Trough Eastern Boundary Gridded Temperature')
text(KN221,0,'KN221','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(PE399,0,'PE399','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(DY053,0,'DY053','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(dy078,0,'DY078','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(dy120,0,'DY120','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(ar30,0,'AR30','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(jc238,0,'JC238','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
text(DY181,0,'DY181','Rotation',45,'fontsize', 12, 'FontName','Helvetica');
hold on
plot([KN221 KN221],[0 1800],'--w','LineWidth',1);
plot([PE399 PE399],[0 1800],'--w','LineWidth',1);
plot([DY053 DY053],[0 1800],'--w','LineWidth',1);
plot([dy078 dy078],[0 1800],'--w','LineWidth',1);
plot([dy120 dy120],[0 1800],'--w','LineWidth',1);
plot([ar30 ar30],[0 1800],'--w','LineWidth',1);
plot([jc238 jc238],[0 1800],'--w','LineWidth',1);
plot([DY181 DY181],[0 1800],'--w','LineWidth',1);

% Salinity

ax(2)=subplot(3,1,2)
%[left bottom width height]
ax(2).Position=[0.05 0.35 0.9 0.28];
x=time;y=pressure;z=SG_WEST;
[c,h]=contourf(x,y,z,12,'LineColor','none');
caxis([round(min(min(z)),1) round(max(max(z)),1)]);
C = colorbar;
cmap=cmocean('Haline');
colormap(ax(2),cmap);
C.Label.String = 'Salinity (g kg^{-1})';
C.Ticks=[round(min(min(z)),1):0.1:round(max(max(z)),1)];
C.TickLength=0.03;
hold on
[c,h]=contour(time , pressure, SG_WEST,[35.2 35.4 35.6],'LineColor','k','LineWidth',0.25); 
clabel(c,h,'LabelSpacing',1000,'Color','w')
% xlabel(gca,'Year');
ylabel(gca,'Depth (m)');
% SORT LABELS
[Y,~,~]=datevec(time);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
set(gca,'XTickLabel',[]);
axis ij
% title('Rockall Trough Eastern Boundary Gridded Salinity')
hold on
plot([KN221 KN221],[0 1800],'--w','LineWidth',1);
plot([PE399 PE399],[0 1800],'--w','LineWidth',1);
plot([DY053 DY053],[0 1800],'--w','LineWidth',1);
plot([dy078 dy078],[0 1800],'--w','LineWidth',1);
plot([dy120 dy120],[0 1800],'--w','LineWidth',1);
plot([ar30 ar30],[0 1800],'--w','LineWidth',1);
plot([jc238 jc238],[0 1800],'--w','LineWidth',1);
plot([DY181 DY181],[0 1800],'--w','LineWidth',1);

% Density

ax(3)=subplot(3,1,3);
%[left bottom width height]
ax(3).Position=[0.05 0.05 0.9 0.28];
x=time;y=pressure;z=gsw_rho(SG_WEST,TG_WEST,0); 
[c,h]=contourf(x,y,z,12,'LineColor','none');
caxis([round(min(min(z)),1) round(max(max(z)),1)]);
C = colorbar;
cmap=cmocean('Dense');
colormap(ax(3),cmap);
C = colorbar;
C.Label.String = 'Potential density (kg/m^{3})';
C.Ticks=[round(min(min(z)),1):0.2:round(max(max(z)),1)];
C.TickLength=0.03;
hold on
[c,h]=contour(time , pressure, z,[1027.8 1027.5 1027.2],'LineColor','k','LineWidth',0.25); 
clabel(c,h,'LabelSpacing',1000,'Color','w')
xlabel(gca,'Year');
ylabel(gca,'Depth (m)');
% SORT LABELS
[Y,~,~]=datevec(time);
date1=datenum(min(Y),1,1);
date2=datenum(max(Y),1,1);
datetick('x','Keeplimits');
axis ij
% title('Rockall Trough Eastern Boundary Gridded Density')
hold on
plot([KN221 KN221],[0 1800],'--w','LineWidth',1);
plot([PE399 PE399],[0 1800],'--w','LineWidth',1);
plot([DY053 DY053],[0 1800],'--w','LineWidth',1);
plot([dy078 dy078],[0 1800],'--w','LineWidth',1);
plot([dy120 dy120],[0 1800],'--w','LineWidth',1);
plot([ar30 ar30],[0 1800],'--w','LineWidth',1);
plot([jc238 jc238],[0 1800],'--w','LineWidth',1);
plot([DY181 DY181],[0 1800],'--w','LineWidth',1);

% save figure
width=35; height=20; FS=12; FN='Arial';
for hh=1:numel(ax)
    set(ax(hh),'fontsize', FS, 'FontName',FN);
end
print('-dpng',[thredds 'RTWB_MERGED_MCAT'])

