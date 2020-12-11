
namefig1 = [netcdffilepath '_ptc'];
mcatncdata.time = mcatncdata.TIME; %repmat(nortekncdata.TIME',length(nortekncdata.DEPTH),1);

% ------------------------------
% Plotting section
% ------------------------------
bad_data=find(mcatncdata.PRES==99999); mcatncdata.PRES(bad_data)=nan;
bad_data=find(mcatncdata.TEMP==99999); mcatncdata.TEMP(bad_data)=nan;
bad_data=find(mcatncdata.CNDC==99999); mcatncdata.CNDC(bad_data)=nan;
bad_data=find(mcatncdata.PSAL==99999); mcatncdata.PSAL(bad_data)=nan;

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

fig=figure('visible',vsblfig,'position',scrsz);
set(fig,'PaperUnits','centimeters','PaperOrientation','portrait',... 
				    'Paperposition',figpos)
subplot(3,1,1)
plot(mcatncdata.time, mcatncdata.PRES);
title([moor ' P (db)'])
datetick

subplot(3,1,2)
plot(mcatncdata.time, mcatncdata.TEMP);
title([moor ' T (C)'])
datetick

subplot(3,1,3)
plot(mcatncdata.time, mcatncdata.PSAL);
title([moor ' SAL'])
datetick

print('-dpng',namefig1)

  