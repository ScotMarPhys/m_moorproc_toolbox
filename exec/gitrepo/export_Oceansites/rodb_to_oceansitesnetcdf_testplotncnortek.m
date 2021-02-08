
namefignor1 = [netcdffilepath '_uv'];
namefignor2 = [netcdffilepath '_pw'];
nortekncdata.time = nortekncdata.TIME; %repmat(nortekncdata.TIME',length(nortekncdata.DEPTH),1);

% ------------------------------
% Plotting section
% ------------------------------
bad_data=find(nortekncdata.UCUR==99999); nortekncdata.UCUR(bad_data)=nan;
bad_data=find(nortekncdata.VCUR==99999); nortekncdata.VCUR(bad_data)=nan;
bad_data=find(nortekncdata.WCUR==99999); nortekncdata.WCUR(bad_data)=nan;
bad_data=find(nortekncdata.PRES==99999); nortekncdata.PRES(bad_data)=nan;

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


fig=figure('visible',vsblfig,'position',scrsz);
set(fig,'PaperUnits','centimeters','PaperOrientation','portrait',... 
				    'Paperposition',figpos)
subplot(2,1,1)
plot(nortekncdata.time, nortekncdata.UCUR);
title([moor ' U (m.s-1)'])
datetick

subplot(2,1,2)
plot(nortekncdata.time, nortekncdata.VCUR);
title([moor ' V (m.s-1)'])
datetick

print('-dpng',namefignor1)

    
fig=figure('visible',vsblfig,'position',scrsz);
set(fig,'PaperUnits','centimeters','PaperOrientation','portrait',... 
				    'Paperposition',figpos)

subplot(2,1,1)
plot(nortekncdata.time, nortekncdata.WCUR);
title([moor ' W (m.s-1)'])
datetick

subplot(2,1,2)
plot(nortekncdata.time, nortekncdata.PRES);
title('pres (db)')
datetick


print('-dpng',namefignor2)

    