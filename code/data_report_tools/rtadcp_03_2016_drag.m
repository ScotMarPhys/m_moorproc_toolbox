% rtadcp1_03_2016_drag.m
% Plot cruise track and cableout details during the drag operation to
% recover this mooring
% DY078, Stuart

% close all ; clear all ;
cd ~/osnap/exec/dy078/data_report_tools ;

[d,h]=mload('~/cruise/data/nav/posmvpos/pos_dy078_01.nc','/') ;
dnum = datenum(2017,1,1,0,0,d.time);
dnum_start = datenum(2017,5,13,13,30,0);
[~,is] = min(abs(dnum-dnum_start))

dnum_end = datenum(2017,5,13,17,00,0);
[~,ie] = min(abs(dnum-dnum_end))

%% load cabledata
mdatapuptechsas(17,133,130000,17,133,190000,'','winch','winch_rtadcp_03_2016_drag','-');

[dw,hw] =  mload('winch_rtadcp_03_2016_drag.nc','/');
dnum_w= datenum(2017,5,13,0,0,dw.time);

%%
% ADCP seabed position
lat_ad = 57.0983 ; lon_ad =	-9.3551 ;
% find ticks at 15min intervals starting from is
is15 = is:900:ie;

figure(1);clf;hold on;grid on;set(gca,'Layer','top');
plot(d.long(is:ie),d.lat(is:ie),'k','Linewidth',2);
plot(d.long(is15),d.lat(is15),'ro','Linewidth',2,...
    'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',7); % 15 minute ticks
plot(lon_ad,lat_ad,'k*','Linewidth',2,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10); % ADCP Seabed position
text(d.long(is15(3))+0.0005,d.lat(is15(3)),'1400','FontSize',12);
text(d.long(is15(7)),d.lat(is15(7))+0.0009,'1500','FontSize',12);
text(d.long(is15(11))+0.0008,d.lat(is15(11)),'1600','FontSize',12);
text(d.long(is15(15))+0.0008,d.lat(is15(15)),'1700','FontSize',12);

scale=1/cos((pi()/180)*57.1);
daspect([scale 1 1]);

xlabel('Longitude [°W]');ylabel('Latitude [°N]');title('Drag course and times');

print -djpeg f1_rtadcp1_03_2016_drag_track.jpg
%%
figure(2);clf;
subplot(2,1,1);hold on;grid on;set(gca,'Layer','top');
plot(dnum_w,dw.cablout,'k','linewidth',2);
set(gca,'YDIR','reverse');
datetick('x');
ylim([0 4000]);
title('Deep tow cable');
xlabel('time [hrs]');ylabel('cable out [m]');

subplot(2,1,2);hold on;grid on;set(gca,'Layer','top');
plot(dnum_w,dw.tension,'k','linewidth',2);
set(gca,'YDIR','reverse');
datetick('x');
ylim([0 3]);
title('Deep tow cable');
xlabel('time [hrs]');ylabel('Tension [tons]');

print -djpeg f2_rtadcp1_03_2016_drag_cable.jpg



