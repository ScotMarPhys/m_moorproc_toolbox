function plot_positions(dti1,dti2)
% Find the positions and depths at given times
m_setup;
this_cruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
rootdir = ['/noc/users/pstar/rpdmoc/rapid/data/exec/' this_cruise];
addpath([rootdir '/mfiles/misc']);

nav_stream = 'posmvpos';
dep_stream = 'ea600m';
no_fixes = 2;

tme(1) = datenum(dti1);
tme(2) = datenum(dti2);
for i=1:2
	tmt(i) = tme(i) - MEXEC_G.Mtechsas_torg;
end

% Now get data from Techsas
ddtt = 3/24;
pos = mtload(nav_stream,datevec(tme(1)-ddtt),datevec(tme(end)+ddtt),'time lat long gndspeed gndcourse','q');
pos.lon = pos.long;
pos.tme = pos.time + MEXEC_G.Mtechsas_torg;
lat = interp1(pos.time,pos.lat,tmt);
lon = interp1(pos.time,pos.lon,tmt);
dep = mtload(dep_stream,datevec(tme(1)-ddtt), datevec(tme(end)+ddtt),'time depthm','q');    
dep.snd = dep.depthm;
iabs = dep.snd == 0; 
wd2 = interp1(dep.time(~iabs),dep.snd(~iabs),tme);

% Work corrected water depths
for i=1:no_fixes
    corr_struct = 	mcarter(lat(1),lon(1),wd2(i));
    wd_corr(i) = corr_struct.cordep;
end    

% Average ships speed
disa2b = sw_dist(lat,lon);
itx = pos.tme > tme(1) & pos.tme < tme(2);
latp = [lat(1) pos.lat(itx) lat(2)];
lonp = [lon(1) pos.lon(itx) lon(2)];
tima2b = (24*(tme(2)-tme(1)));
dispath = sum(sw_dist(latp,lonp));
spda2b = disa2b/(24*(tme(2)-tme(1)));
spdpath = dispath/(24*(tme(2)-tme(1)));
iout = 1;
fprintf(iout,'\n Time   Lat   Lon  Uncorr Depth Corr depth \n');
for i = 1:no_fixes
  fprintf(iout,' %s  %8.4f %8.4f %7.1f % 7.1f \n', ...
          datestr(tme(i)),lat(i),lon(i),wd2(i),wd_corr(i));
  fprintf(iout,'        %8.2f N   %8.2f E \n',dd2dm(lat(i)),dd2dm(lon(i)));
end

fprintf(1,'Elapsed time = %7.1f hours \n',tima2b);
fprintf(1,'A to B distance = %7.1f Nm \n',disa2b);
fprintf(1,'Path distance = %7.1f Nm \n',dispath);
fprintf(1,'A to B average speed = %7.1f kts \n',spda2b);
fprintf(1,'Path average speed = %7.1f kts \n',spdpath);

tdif = tme(2)-tme(1);
aaa = 1+floor(0.001+(tdif*24/3)/8);
xtic = [floor(24*tme(1))/24:aaa*3/24:ceil(24*tme(2))/24];
%datestr(xtic)
figure
subplot(3,1,1)
plot(pos.tme,pos.gndspeed);
hold on
y1 = ylim; t1 = tme(1);
for i = 1:no_fixes
  t1 = tme(i);
  plot([t1 t1],y1,'r');
end
xlm = [tme(1)-1/24 tme(2)+1/24];
set(gca,'XTick',xtic); xlim(xlm) ;
datetick('x','dd-HH', 'keeplimits','keepticks');
xlabel('Time'),ylabel('Speed');grid minor
subplot(3,1,2)
plot(pos.tme,pos.lon);
hold on
y1 = ylim; t1 = tme(1);
for i = 1:no_fixes
  t1 = tme(i);
  plot([t1 t1],y1,'r');
end
set(gca,'XTick',xtic); xlim(xlm);
datetick('x','dd-HH', 'keeplimits','keepticks');
xlabel('Time'),ylabel('Longitude');grid minor
subplot(3,1,3)
plot(pos.tme,pos.lat);
hold on
y1 = ylim; t1 = tme(1);
for i = 1:no_fixes
  t1 = tme(i);
  plot([t1 t1],y1,'r');
end
set(gca,'XTick',xtic); xlim(xlm);
datetick('x','dd-HH', 'keeplimits','keepticks');
xlabel('Time'),ylabel('Latitude');grid minor


figure
plot(pos.lon,pos.lat,'b');
hold on
grid on
hb =plot(lon,lat,'r+','MarkerSize',12);
rlon = max(pos.lon)-min(pos.lon);
rlat = max(pos.lat)-min(pos.lat);
mlon = 0.5*(max(pos.lon)+min(pos.lon));
mlat = 0.5*(max(pos.lat)+min(pos.lat));
rx = max(rlon,rlat);
xlm = mlon + 0.5*[-rx rx];
ylm = mlat + 0.5*[-rx rx];
xlim(xlm);
ylim(ylm);