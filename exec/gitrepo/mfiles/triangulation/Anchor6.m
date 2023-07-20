function Anchor6
% This version allows for chanigng depth of target
%
% Find the release position from 3 or more positions and ranges
% This version user enters only time information in file  <moorname>_times.txt
% created by the user and stored in:
%   /noc/users/pstar/rpdmoc/rapid/data/moor/raw/<cruise>/moor_positions/
%
% Results written and position is found from navigation files on Techsas.
%
% Format of text file should be 
%    YYYY MM DD HH MM SS range depth(if applicable)
%
% First line is anchor drop:
%   - enter zero for range
%   - and NEGATIVE of uncorrected water depth at time of anchor drop
%   - if no depth data are available, enter a POSITIVE number as best
%   estimate of depth
%
% following lines are traingulation times (enter zero for depth)
%
% e.g.
%
% ebh4l6_times.txt:
% 2015 10 28 10 28 20 0 -1014
% 2015 10 28 10 46 10 1157 0
% 2015 10 28 11 00 30 1549 0
% 2015 10 28 11 11 30 1547 0
%
% Also solves for traingulation position.
%
% Saves restults to /noc/users/pstar/rpdmoc/rapid/data/moor/raw/<cruise>/moor_positions/
%
% This version uses calls to m_map toolbox to use m_fdist and m_idist to calculate
% distances on spheroid.
%
% DAS October 2012 Updated for JC103 May 2014
%

% Techsas stream to get navigation info from:

global MEXEC_G MOORPROC_G

mcruise = MOORPROC_G.cruise;
opt1 = 'ship'; opt2 = 'datasys_best'; get_cropt %default_navstream
nav_stream = default_navstream;
dirout = fullfile(MOORPROC_G.moordatadir,'raw',mcruise,'moor_positions');

% User input
fprintf(1,'\n Enter mooring name (e.g. ebh3)). Times will then be read from <name>_times.txt: \n');
fprintf(1,' - in the folder %s \n',dirout)
fprintf(1,' - with one range per line in format  YYYY MM DD HH MM SS range. \n')
fprintf(1, ' Output will be saved to <name>_triangle.txt in the same directory. \n\n');
loc_name = input('Enter mooring name (e.g. ebh3): ','s');

if isempty(loc_name)
  return
else
  if exist(dirout,'dir')
    filein = fullfile(dirout, [loc_name '_times.txt']);
    fileout = fullfile(dirout, [loc_name '_triangle.txt']);
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

td  = input('Enter ship transducer depth in metres: ');
dep0 = input('Enter approximate depth(m) at time of first range: ');
was = input('Enter rate of ascent (m/hr): ');

% Read in data from text file assume first point is anchor drop
tme = nan(1,size(indata,1));
range = tme;
dptvl = tme;
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
switch MEXEC_G.Mshipdatasystem
    case 'techsas'
        pos = mtload(nav_stream,datevec(tme(1)-0.04), datevec(tme(end)+0.01),'time lat long');
        pos.lon = pos.long;
        pos.time = pos.time + uway_torg;
    case 'rvdas'
        pos = mrload(nav_stream,datevec(tme(1)-0.04), datevec(tme(end)+0.01),'time latitude longitude');
        pos.lat = pos.latitude;
        pos.lon = pos.longitude;
        pos.time = pos.dnum;
    case 'scs'
        timvar = 'time'; latvar = 'GPS-Furuno-GGA-lat'; lonvar = 'GPS-Furuno-GGA-lon';
        pos = msload(nav_stream,datevec(tme(1)-0.04), datevec(tme(end)+0.01),[timvar ' ' latvar ' ' lonvar]);
        latvar = replace(latvar,'-','_'); lonvar = replace(lonvar,'-','_');
        pos.lon = pos.(lonvar);
        pos.lat = pos.(latvar);
        pos.time = pos.(timvar) + uway_torg;
end

% interpolate for positions then depths
lat = interp1(pos.time,pos.lat,tme);
lon = interp1(pos.time,pos.lon,tme);


% Separate the anchor drop position
latD = lat(1);
lonD = lon(1);
tmeD = tme(1);
range = range(2:end);
lat = lat(2:end);
lon = lon(2:end);
tme = tme(2:end);
no_fixes=length(range);
if isempty(was)
    dep = dep0
else
    dep = dep0 - (tme-tme(1))*24*was
end

% Work out effective horizontal range
% Not sure about sound speed correction - i think option 2 is correct
 
rangeh = sqrt(range.^2 - (dep-td).^2);

% Open a new plot
figure
% Decide on boudaries for plot
mpd = 111.2*1000;
axylim = 1.5;
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
  [clon,clat,~] = m_fdist(lon(i),lat(i),xdegr,rangeh(i));
  clon = clon-360;
  if i < 4 
    plot(lon(i),lat(i),'b+');
  else
	plot(lon(i),lat(i),'r+')
  end
  if i < 4 
    plot(clon,clat,'b');
  else
	plot(clon,clat,'r');
  end
end

titletext1=['Triangulation Survey for: ',loc_name];
titletext2=sprintf('Releases depth: %5.0f m. rate of ascent: %5.1f (m/hour)', dep0,was);
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

% Plot to file is saving
if strmatch(loc_name, '')
  return
else
  plotout = [dirout loc_name '_triangle'];
  print('-depsc', plotout)
  eval(['!cat ' fileout])
end

% Now data that we save ni the output file
fprintf(iout,'%s  %s \n',loc_name,datestr(tmeD,1));
fprintf(iout,'Anchor drop at: %s %8.4f %8.4f \n', ...
              datestr(tmeD,31),latD,lonD);
fprintf(iout,'Date      Time       Lat   Lon  Slant range  Horiz range Residaul \n');
for i = 1:no_fixes
  fprintf(iout,'%s  %8.4f %8.4f %7.0f %7.0f %7.0f \n', ...
        datestr(tme(i),31),lat(i),lon(i),range(i),rangeh(i),reser(i));
end
fprintf(iout,'Trilaterated position \n');
fprintf(iout,'%8.4f %8.4f \n',APlat,APlon);
fprintf(iout,'%s \n',title4);
fprintf(iout,'Fallback  %s m \n',fallback);
fprintf(iout,'Comments: %s \n',us_com);
