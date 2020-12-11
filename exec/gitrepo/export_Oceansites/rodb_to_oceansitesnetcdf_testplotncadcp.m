
namefigadcp1 = [netcdffilepath '_uv'];
namefigadcp2 = [netcdffilepath '_err'];

ADCPncdata.time = repmat(ADCPncdata.TIME',length(ADCPncdata.BINDEPTH),1);

% ------------------------------
% Plotting section
% ------------------------------

%----------------------------------------------------
% Graphical parameter
vsblfig = 'off';
zsc = get(0,'MonitorPositions');
scrsz = [1 1 1900 1400]; % zsc(1,:);
figpos = scrsz/60; %[0 0 27 21];
fs1 = 12;
fs2=10;
set(0,'DefaultAxesFontName', 'Helvetica')
set(0,'DefaultAxesFontSize', fs1)
set(0,'DefaultTextFontname', 'Helvetica')
set(0,'DefaultTextFontSize', fs2)
%-----------------------------------------------------
cmapvel = lansey;

proc = 1;    
plotylim = [min(-ADCPncdata.PRES(:)) max(-ADCPncdata.PRES(:))];

fig=figure('visible',vsblfig,'position',scrsz);
set(fig,'PaperUnits','centimeters','PaperOrientation','portrait',... 
				    'Paperposition',figpos)
subplot(2,1,1)
pcolor(ADCPncdata.time, -ADCPncdata.PRES,ADCPncdata.UCUR);
set(gca,'ylim', plotylim)
shading flat
view([0 90])
caxis([-0.5 0.5])
title([moor ' U (m.s-1)'])
datetick
colormap(cmapvel)
colorbar

subplot(2,1,2)
pcolor(ADCPncdata.time, -ADCPncdata.PRES,ADCPncdata.VCUR);
set(gca,'ylim', plotylim)
shading flat
view([0 90])
caxis([-0.5 0.5])
title([moor ' V (m.s-1)'])
datetick
colormap(cmapvel)
colorbar

print('-dpng',namefigadcp1)

    
fig=figure('visible',vsblfig,'position',scrsz);
set(fig,'PaperUnits','centimeters','PaperOrientation','portrait',... 
				    'Paperposition',figpos)


% subplot(3,1,1)
% surf(ADCPncdata.time, -ADCPncdata.PRES,ADCPncdata.PG4+ADCPncdata.PG4);
% set(gca,'ylim', plotylim)
% shading flat
% view([0 90])
% caxis([0 100])
% title('PG4 - pgood 4 beam solutions')
% datetick
% colormap(gca,parula)
% colorbar
%     
% subplot(3,1,2)
% surf(ADCPncdata.time, -ADCPncdata.PRES,ADCPncdata.PG4+ADCPncdata.PG1);
% set(gca,'ylim', plotylim)
% shading flat
% view([0 90])
% caxis([0 100])
% title('PG1 - pgood 3 beam solutions')
% datetick
% colormap(gca,parula)
% colorbar

subplot(3,1,3)
surf(ADCPncdata.time, -ADCPncdata.PRES,ADCPncdata.ECUR);
set(gca,'ylim', plotylim)
shading flat
view([0 90])
caxis([0 0.05])
title('Error Velocity (m.s-1)')
datetick
colormap(gca,parula)
colorbar

print('-dpng',namefigadcp2)

    