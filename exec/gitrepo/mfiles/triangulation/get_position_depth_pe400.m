% function get_position_depth_pe399
% Find the positions and depths at given times, modified by Loic Houpert
% for the cruise PE399 26/06/2015

m_setup;
this_cruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
rootdir = ['/home/mstar/osnap/exec/' this_cruise];
addpath([rootdir '/mfiles/misc']);
nav_stream = 'posmvpos'; %penav
dep_stream = 'ea600m';

disp('Please enter date of at least two location: ')
yy = input('Year: ');
MM = input('Month (MM): ');
dd = input('Day (dd): ');

disp('Now enter times (UTC): ')
ient = 1; i = 0;

while ient == 1
    i = i+1;
    disp(['Enter time of position number ' num2str(i-1)])
    stime = input('Time (hh:mm) or (hh:mm:ss) or press return to finish: ','s');
    if strmatch(stime,'');
      ient = 0;
    else
      ctime(i) = {stime};
      hh(i) = str2num(stime(1:2));
      mm = str2num(stime(4:5));
      if length(stime) == 8
        ss = str2num(stime(7:8));
      else
        ss = 0;
      end  
      tme(i) = datenum([yy MM dd hh(i) mm ss]) ;
    end
end


% how many points
no_fixes=length(tme);

% in case we go over midnight
for i = 2:no_fixes
  if tme(i) < tme(i-1)
    tme(i) = tme(i)+1;
	fprintf(1,'Assuming monotonic time adding one day!')
  end
end


% Techsas has different time origin
tmm = tme - MEXEC_G.Mtechsas_torg;
dtt = 0.1*(tme(end)-tme(1));

% Now get data from Techsas
pos = mtload(nav_stream,datevec(tme(1)-dtt), ...
                        datevec(tme(end)+dtt));
pos.lon = pos.long;
pos.tme = pos.time + MEXEC_G.Mtechsas_torg;
lat = interp1(pos.tme,pos.lat,tme);
lon = interp1(pos.tme,pos.lon,tme);

dep = mtload(dep_stream,datevec(tme(1)-dtt), ...
                        datevec(tme(end)+dtt));    
dep.snd = dep.depthm; % uncomment for dy053
dep.tme = dep.time + MEXEC_G.Mtechsas_torg;
iabs = dep.snd == 0; 
wd2 = interp1(dep.tme(~iabs),dep.snd(~iabs),tme);

% Work corrected water depths
for i=1:no_fixes
    corr_struct = 	mcarter(lat(1),lon(1),wd2(i));
    wd_corr(i) = corr_struct.cordep;
end    

iout = 1;
fprintf(iout,'Time   Lat   Lon  Uncorr Depth Corr depth \n');
for i = 1:no_fixes
  fprintf(iout,' %s  %8.4f %8.4f %7.1f % 7.1f \n', ...
          char(ctime(i)),lat(i),lon(i),wd2(i),wd_corr(i));
  fprintf(iout,'        %8.2f N   %8.2f E \n',dd2dm(lat(i)),dd2dm(lon(i)));
end

figure
subplot(3,1,1)
plot(dep.tme(~iabs),dep.snd(~iabs));
hold on
y1 = ylim; 
for i = 1:no_fixes
  t1 = tme(i);
  plot([t1 t1],y1,'r');
end
datetick;xlabel('Time'),ylabel('Depth');grid on
subplot(3,1,2)
plot(pos.tme,pos.lon);
hold on
y1 = ylim; t1 = tme(1);
for i = 1:no_fixes
  t1 = tme(i);
  plot([t1 t1],y1,'r');
end
datetick;xlabel('Time'),ylabel('Longitude');grid on
subplot(3,1,3)
plot(pos.tme,pos.lat);
hold on
y1 = ylim; t1 = tme(1);
for i = 1:no_fixes
  t1 = tme(i);
  plot([t1 t1],y1,'r');
end
datetick;xlabel('Time'),ylabel('Latitude');grid on

datevec(tme)
