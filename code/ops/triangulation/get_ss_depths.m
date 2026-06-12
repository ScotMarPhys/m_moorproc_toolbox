function [depths] = get_ss_depths(lats,lons)

 lats = [27.8152,27.8281,27.8453,27.8755,27.8970,27.9272,27.9659,28.0089,28.0219,28.0348];
 lons = [ 346.6720,346.6087,346.5307,346.4089,346.3017,346.1847,346.0141,345.8484,345.7900,345.7461]-360;

 dis = sw_dist(lats,lons);
 dis = [0 dis];
 cdis =cumsum(dis);
 
if min(lats) > 90 | min(lats) < -90
	lats = dm2dd(lats);
	lons = dm2dd(lons);
end

if length(lats) > 1
	dlat = max(lats)-min(lats);
	dlon = max(lons)-min(lons);
	dd = max([0.1*dlat 0.1*dlon 0.02]);
else
	dd = 0.02
end

east = max(lons)+dd;west = min(lons)-dd;
north = max(lats)+dd;south = min(lats)-dd;

[elevations,elat,elon]=mygrid_sand([south,north,west,east],1);
elon = elon - 360;

%x =elat'*ones(size(elon));
%y = ones(size(elat'))*elon;

[yy,xx] = meshgrid(elat,elon);
depths = interp2(yy,xx,elevations',lats,lons);

spd = 7.5;
tsetup = 0.25;
tme = dis/spd + tsetup + 2*abs(depths)/3600;
ctme = cumsum(tme);

for i = 1:length(lats)
	fprintf(1,'%6.2f N %6.2f W Depth %6.1f m \n',lats(i),lons(i),depths(i))
end
fprintf(1,' \n \n')

for i = 1:length(lats)
%	fprintf(1,'%7.2f N %7.2f W Depth %6.1f m \n',dd2dm(lats(i)),dd2dm(lons(i)),depths(i))
	fprintf(1,'%3i %7.2f N %7.2f W Depth %8.1f m Dist %5.1f nm  %5.1f nm %5.1f hours %5.1f hours \n', ...
	        i,dd2dm(lats(i)),dd2dm(lons(i)),depths(i),dis(i),cdis(i),tme(i),ctme(i))
end

m_proj('lambert','lon',[west, east],'lat',[south, north]);
m_grid('box','on','color','k','linewidth',[1],'fontsize',[12]);
hold on
m_pcolor(elon,elat,elevations); shading flat
h = m_plot(lons,lats,'k+','LineWidth',2.5,'MarkerSize',12);

