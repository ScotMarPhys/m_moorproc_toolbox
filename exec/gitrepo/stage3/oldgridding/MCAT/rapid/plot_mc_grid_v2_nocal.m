% Plot 2015 uncalibrated microcat grid WB array (WB1+WB2)
% Grid = rtwb_02_2015_grid_v2_nocal.mat
% esdu, Jul 16

close all; clear;

moor = 'rtwb_02_2015';

addpath C:\Code\Matlab\Various\nansuite
addpath C:\CRUISES\DY053\Moorings\matlab\tools
filenm=['C:\CRUISES\DY053\Moorings\processed\hydro_grid_v2_nocal\' moor '_grid_v2_nocal.mat'];
load(filenm);

% Convert mstar time to matlab time
time_str=gregorian(jd);
mtime=datenum(time_str);

% Calculate sigma-theta
sigt=nan(size(TGfs));
for j=1:length(TGfs)
    sigt(:,j)=sw_pden(SGfs(:,j),TGfs(:,j),p_grid,0);
end


%%
scrsz = get(0,'ScreenSize');
fig=figure('Position',[1 1 scrsz(3) scrsz(4)]);

p1=subplot(3,1,1);
imagesc(mtime,p_grid,TGfs)
hold on
%plot(mtime,Pfs(:,:),'.k')
datetick('x',12); ylim([0 1800]);
set(gca,'fontsize',8)
caxis([3.5 11.5]);
colorbar
title(['Temperature OSNAP ' moor ' - not calibrated'],'interpreter','none');
ylabel('Pressure (db)')

p2=subplot(3,1,2);
imagesc(mtime,p_grid,SGfs)
hold on
%plot(mtime,Pfs(:,:),'.k')
datetick('x',12); ylim([0 1800]);
set(gca,'fontsize',8)
caxis([34.9 35.4]);
colorbar
title(['Salinity  OSNAP ' moor ' - not calibrated'],'interpreter','none');
ylabel('Pressure (db)')

p3=subplot(3,1,3);
imagesc(mtime,p_grid,sigt)
hold on
%plot(mtime,Pfs(:,:),'.k')
datetick('x',12); ylim([0 1800]);
set(gca,'fontsize',8)
caxis([1026.85 1027.85]);
colorbar
title(['Density  OSNAP ' moor ' - not calibrated'],'interpreter','none');
ylabel('Pressure (db)')

linkaxes([p1 p2 p3],'xy')

saveas(fig,['hydrogrid_timeseries_nocal_' moor '.tif']);
close