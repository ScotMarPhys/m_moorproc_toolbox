function Anchor5
% Find the anchor position from 3 or more positions and ranges
% This version user enters only time information in file  <name>_times.txt
% Results written
% and position is found from navigation files on Techsas.
% Format of text file shoudl be 
%    YYYY MM DD HH MM SS range  
% First line is anchor drop (enter zero for range) and followign lines are traingulation times
% Also solves for traingulationa position.
% Saves restults to /noc/users/pstar/rpdmoc/rapid/data/moor/raw/d382/moor_positions/
% This version uses calls to m_map toolbox to use m_fdist and m_idist to calculate
% distances on spheroid.
% DAS October 2012 Updated for JC103 May 2014

% Where to get navigation info from
%nav_stream = 'gpsfugro'  % D382
nav_stream = 'posmvpos';
dep_stream = 'ea600m';

m_setup
this_cruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
rootdir = ['/noc/users/pstar/rpdmoc/rapid/data/exec/' this_cruise];

% To access other mfiles
addpath([rootdir '/mfiles/misc']);

dirout = ['/noc/users/pstar/rpdmoc/rapid/data/moor/raw/' this_cruise '/moor_positions/'];
%dirout = '/noc/users/pstar/rpdmoc/rapid/data/moor/raw/d382/moor_positions/';

% info from user
fprintf(1,'\n Enter mooring name (e.g. ''ebh3'')). Times will then be read from <name>_times.txt: \n');
fprintf(1,' - in the folder %s \n',dirout)
fprintf(1,' - with one range per line in format  YYYY MM DD HH MM SS range. \n')
fprintf(1, ' Outpur will be saved to <name>_triangle.txt in the smae directory. \n\n');
loc_name = input('Enter mooring name (e.g. ''ebh3''): ');

if strmatch(loc_name, '')
  return
else
  if exist(dirout,'dir')
    filein = [dirout loc_name '_times.txt'];
    fileout = [dirout loc_name '_triangle.txt'];
    disp(['Reading times from ' filein]);
    if exist(filein,'file')
      indata = load(filein);
    else
      disp('input file does not exst')
      return
    end
    disp(['Opening output file ' fileout]);
    disp('If file exists contents will be overwritten')
    iout = fopen(fileout,'w+');
  else
    disp(['Directory ' dirout ' does not exist']);
    return
  end
end

rht=input('Enter approximate height of release above seabed in metres: ');
td =input('Enter approximate transducer depth in metres: ');

% Read in data from text file assume first point is anchor drop
for i = 1:size(indata,1)
  dvec = indata(i,1:6);
  tme(i) = datenum(dvec);
  if i == 1
    range(1) = -999;
  else
    range(i) = indata(i,7);
  end
  if size(indata,2) == 8
    dptvl(i) = indata(i,8);
  else 
	dptvl(i) = -999;
  end
end

% Now get data from Techsas
pos = mtload(nav_stream,datevec(tme(1)-0.04), ...
                        datevec(tme(end)+0.01),'time lat long','q');
pos.lon = pos.long;
pos.time = pos.time + MEXEC_G.Mtechsas_torg;

% interpolate for positions then depths
lat = interp1(pos.time,pos.lat,tme);
lon = interp1(pos.time,pos.lon,tme);

% information on depth
if dptvl(1) < 0
  dep = mtload(dep_stream,datevec(tme(1)-0.01), ...
                      datevec(tme(1)+0.01),'time depthm','q');   
  dep.time = dep.time + MEXEC_G.Mtechsas_torg;
% Need to check the depths
  iabs = dep.depthm == 0; 
  wd2 = interp1(dep.time(~iabs),dep.depthm(~iabs),tme);
  f1 = figure;
  plot(dep.time(~iabs),dep.depthm(~iabs));
  hold on
  y1 = ylim; t1 = tme(1);
  plot([t1 t1],y1,'r');
  datetick;xlabel('Time'),ylabel('Depth');
  fprintf(1,'Uncorrected water depth at anchor drop %5.1f m \n',wd2(1));
  ich = input('Enter 1 to accept or 2 to change values: ');
  if ich == 2
    wd=input('Enter uncorrected water depth in metres: ');
  else
    wd = wd2(1);
  end
  close(f1)
% The corrected water depth
  corr_struct = 	mcarter(lat(1),lon(1),wd);
  wd_corr = corr_struct.cordep;
else
	wd0 = dptvl(1);
	%wd_corr = carter(latM,lonM,wd);
	corr_struct = 	mcarter(lat(1),lon(1),wd0);
	wd_corr = corr_struct.cordep;
	wd = 2*wd0-wd_corr;
	wd_corr = wd0;
end


% Separate the anchor drop position
latD = lat(1);
lonD = lon(1);
tmeD = tme(1);
range = range(2:end);
lat = lat(2:end);
lon = lon(2:end);
tme = tme(2:end);
no_fixes=length(range);

% If no fixes then just give anchor drop position
if no_fixes == 0
	fprintf(1,'Anchor drop at: %s %8.4f %8.4f Corr. water depth: %6.1f \n', ...
	              datestr(tmeD,31),latD,lonD,wd_corr);
	fprintf(1,'%8.2f N %8.2f W \n',dd2dm(latD),dd2dm(lonD))

	fprintf(iout,'%s  %s \n',loc_name,datestr(tmeD,1));
	fprintf(iout,'Anchor drop at: %s %8.4f %8.4f Corr. water depth: %6.1f \n', ...
	              datestr(tmeD,31),latD,lonD,wd_corr);
	fprintf(iout,'%8.4f %8.4f \n',dd2dm(latD),dd2dm(lonD));
	latdeg = floor(latD);
	londeg = floor(-lonD);
	latmin = 60*(latD-latdeg);
	lonmin = 60*(-lonD-londeg);

	title4 = sprintf('Latitude %i %5.2f N, Longitude %i %5.2f W',latdeg,latmin,londeg,lonmin);
	fprintf(iout,'%s',title4);
	return
end

% Work out effective horizontal range
% Not sure about sound speed correction - i thin option 2 is correct
irn = 3;
if max(range) > 8000
	% Otherwise mcarter fails
	irn = 3
end
for i=1:no_fixes
% Horizontal range from ship to position vertically above release
  if irn == 1
    rangeh(i)=sqrt((range(i))^2 - (wd_corr-rht-td)^2); 
  elseif irn == 2
    corr_struct = 	mcarter(lat(1),lon(1),range(i));
    rangeC = corr_struct.cordep;
	rangeh(i)=sqrt((rangeC)^2 - (wd_corr-rht-td)^2); 
  elseif irn == 3
    rangeh(i)=sqrt((range(i))^2 - (wd-rht-td)^2); 
  end
end    

% Open a new plot
figure
% Decide on boudaries for plot
mpd = 111.2*1000;
axylim = 1.05;
south = latD -axylim*max(rangeh)/mpd;
north = latD +axylim*max(rangeh)/mpd;
west =  lonD -axylim*max(rangeh)/(mpd*cos(latD*pi/180));
east =  lonD +axylim*max(rangeh)/(mpd*cos(latD*pi/180));
m_proj('lambert','lon',[west, east],'lat',[south, north]);
% Following is if use m_map to plot - but this does not allow zooming so have disabled
% m_grid('box','on','color','k','linewidth',[1],'fontsize',[14]);
hold on
% plot ships track
plot(pos.lon,pos.lat,'k--')

xdegr = -180:2:180;
for i=1:no_fixes
  [clon,clat,az3] = m_fdist(lon(i),lat(i),xdegr,rangeh(i));
  clon = clon-360;
  if i < 6 
    plot(lon(i),lat(i),'b+');
  else
	plot(lon(i),lat(i),'r+')
  end
  if i < 6 
    plot(clon,clat,'b');
  else
	plot(clon,clat,'r');
  end
end

titletext1=['Triangulation Survey for: ',loc_name];
titletext2=sprintf( ...
    'Corrected water depth: %5.0f m. Release Height: %3.0f m. Transducer depth: %3.0f', ...
    wd_corr,rht,td);
title({titletext1;titletext2});
xlabel('Longitude'); ylabel('Latitude')


% Calculate triangulation point by solving least squares problem
% Nb take differences of equations to make linear

lon0 = mean(lon);
lat0 = mean(lat);
lat0 = latD;
lon0 = lonD;

[long1,latg1,x,y,reser0] = solve_anchor(lon0,lat0,lon,lat,rangeh,no_fixes);
[APlon,APlat,x,y,reser] = solve_anchor(long1,latg1,lon,lat,rangeh,no_fixes);

if max(reser0-reser) > 0.1
	fprintf(1,'Something not right')
    keyboard
end

% plot estimate and anchor drop point
h1 = plot(APlon,APlat,'mo');
text((east-west)*0.1+west,(north-south)*0.90+south,['Anchor drop'],'color','g');
plot(lonD,latD,'g+');
plot(lonD,latD,'go');
grid on
%hh = plot(tlon,tlat,'MarkerSize',20,'LineWidth',1.5);

% determine anchor position from figure(1)
disp('Estimated position:' )
fprintf(1,'Latitude %8.4f Longitude %8.4f \n',APlat,APlon)
fprintf(1,'Max residual error: %6.1f m and RMS residual error: %6.1f m\n',max(reser),sqrt(sum(reser.^2)/no_fixes))
ich = input('Enter 1 to accept or 2 to change values: ');
if ich == 2
  disp('Use figure to determine anchor seabed position:')
  APlat = input('Latitude = ');
  APlon = input('Longitude = ');
  sol(1) = mpd*(APlon-lon0);
  sol(2) = mpd*(APlat-lat0);
  [dis,az12,az21] = m_idist(lon0,lat0,APlon,APlat);
  az = az12*pi/180;
  sol(1) = dis*sin(az);
  sol(2) = dis*cos(az);
  reser  = abs( sqrt((x-sol(1)).^2+(y-sol(2)).^2)-rangeh);
  fprintf(1,'Max residual error: %6.1f m \n',max(reser))
  set(h1,'Visible','off')
end

% User comments?
us_com = input('Enter any comments to be saved here:','s');

% Degrees and minutes
latdeg = floor(APlat);
londeg = floor(-APlon);
latmin = 60*(APlat-latdeg);
lonmin = 60*(-APlon-londeg);

title4 = sprintf('Latitude %i %5.2f N, Longitude %i %5.2f W',latdeg,latmin,londeg,lonmin);
%fallback1 = sw_dist([APlat lat(1)],[APlon lon(1)],'nm');
%fallback2 = sw_dist([APlat lat(1)],[APlon lon(1)],'km');

plot(APlon,APlat,'r+');
plot(APlon,APlat,'ro');
plot(lonD,latD,'g+');
plot(lonD,latD,'go');
titletext3=['Red = anchor seabed position. ',num2str(APlat),'N ',num2str(APlon),'W.'];
%titletext4=['Fall back = ' num2str(fallback2) ' km = ' num2str(fallback1) ' nm.'];
title({titletext1;titletext2;titletext3;title4});

% Calcualte distance from drop to moored positions
fallback=sw_dist([latD APlat],[lonD APlon],'km')*1000;
fallback=sprintf('%5.0f',fallback);

% Text for plot
text((east-west)*0.1+west,(north-south)*0.95+south,['Fallback = ' fallback 'm'],'color','k')
text((east-west)*0.1+west,(north-south)*0.90+south,['Anchor drop'],'color','g');
text((east-west)*0.1+west,(north-south)*0.85+south,['Anchor location'],'color','r')

% Finish plot
xlim([west east]);
ylim([south north]);

fprintf(iout,'%s  %s \n',loc_name,datestr(tmeD,1));
fprintf(iout,'Anchor drop at: %s %8.4f %8.4f Corr. water depth: %6.1f \n', ...
              datestr(tmeD,31),latD,lonD,wd_corr);
fprintf(iout,'Date      Time       Lat   Lon  Slant range  Horiz range Residaul \n');
for i = 1:no_fixes
  fprintf(iout,'%s  %8.4f %8.4f %7.0f %7.0f %7.0f \n', ...
        datestr(tme(i),31),lat(i),lon(i),range(i),rangeh(i),reser(i));
end
fprintf(iout,'Triangulated position \n');
fprintf(iout,'%8.4f %8.4f \n',APlat,APlon);
fprintf(iout,'%s \n',title4);
fprintf(iout,'Fallback  %s m \n',fallback);
fprintf(iout,'Comments: %s \n',us_com);

% Plot to file is saving
if strmatch(loc_name, '')
  return
else
  plotout = [dirout loc_name '_triangle'];
  print('-depsc', plotout)
  eval(['!cat ' fileout])
end

% Check bathymetry from swath if available - might want to rerun again after this
if strcmp(loc_name(1:2),'wb') & ~strcmp(loc_name(1:3),'wb4');
    topo_map = 'gr_kn182_plot.mat';
elseif strcmp(loc_name(1:2),'ma')
	loc_sh = loc_name(1:4);
    if strmatch(loc_sh,char([{'mar3'}])) > 0
    	topo_map = 'mar34.mat';
    elseif strmatch(loc_sh,char([{'mar2'},{'mar1'}])) > 0
	    topo_map = 'mar12.mat';
    elseif strmatch(loc_sh,char([{'mar0'}])) > 0
	    topo_map = 'mar0_JC064_swath.mat';
	end
else
	topo_map = 'none';
end
if ~strcmp(topo_map,'none')
  rdpath = ['/noc/users/pstar/rpdmoc/rapid/data/exec/' this_cruise '/mfiles/rapid_widgit_v2/data/'];
  load([rdpath topo_map]);
  newdep = interp2(lon,lat,depth,APlon,APlat);
  fprintf(1,'Corrected depth at triangulated positoin is %7.1f \n',newdep);
  fprintf(1,'triangualtion calc used %7.1f \n',wd_corr);
else
    fprintf(1,'No swath data available depth used was %7.1f \n',wd_corr);
end

% Check what bathymetry was when ship passed over the anchor point
ixt = pos.time > tme(1)-0.5 & pos.time < tme(1) +0.5;
lat_sh = pos.lat(ixt);lon_sh = pos.lon(ixt); tme_sh = pos.time(ixt);
dis_sh = 1000*111.2*sqrt( (lat_sh-APlat).^2+cos(APlat*pi/180)^2*(lon_sh-APlon).^2);
im = find(dis_sh == min(dis_sh));
if length(im) > 1
	im = im(1)
end
tp_sh = tme_sh(im);

if dptvl(1) < 0
	  wd_sh = interp1(dep.time(~iabs),dep.depthm(~iabs),tp_sh);
	  corr_struct = mcarter(lat(1),lon(1),wd);
	  wd_corr_sh = corr_struct.cordep;
else
	  wd_corr_sh = 0;
end

fprintf(1,'\n Closest ship track to anchor was %6.1f m \n',dis_sh(im))
fprintf(1,' at time %s \n',datestr(tp_sh))
fprintf(1,' Where depth was %6.1f m \n',wd_corr_sh)  