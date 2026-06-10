function rapid_usbl_lite(varargin)
% This script creates plot showing navigation information.
% *********** PLEASE CHECK POSITIONS ARE CORRECT *************
% The positions are stored in the file:
%  - moor_pos.dat
% See README file for instuctions on displaying the plot on the web.
% Copies final plots to the webserver directory.  E.g.
%  - /srv/www/vhosts/jc103
%
% Plot is to file only and there is no screen plot.
% Updates when advance 5% of distance to waypoint (minimum 10s)
% The script will continue for one day or until user interupts
% Primarily for use during deployment of moorings.
% 
% There are variable number of inputs:  this_moor,plot_size,cont_int
% Usually only need the first one
%   - this_moor - is a string and can have values:
%   	- 'ship'  Will create plot cetred on ship's position
%		-  any of the mooring names listed in moor_pos.dat
% 		- 'check' Will list all mooring names, positions and nominal depths
%   - plot_size - is a real number and is size of plot in nautical miles
%   - cont_int - is contour interval usually do not need to set this
%   - web_dis - default is 1 to copu files to webserver and have no matlab display
%                 but if 0 then will display in Matlab only and will plot jsut once
%
% Calls following scripts:  rapid_web_map.m, dirav.m, dm2dd.m, dd2dm2.m dinterp1.m
%
% Nb the calculation of ETA assumes that hte ship is heading directly 
% towards the waypoint
%
% Update DAS April 2014 on JC103

% Default values of input variables
this_moor = 'check';
plot_size = 0.5;		% Default choice
cont_int = -999;   	% < 0 means automatic choice
web_dis = 1;

% Optional values from user
if nargin >= 1;
	this_moor = char(varargin{1});
	if isempty(this_moor) this_moor = 'ship'; end
end
if nargin >= 2;
	plot_size = varargin{2};
	if isempty(plot_size) plot_size = 8; end
end
if nargin >= 3;
	cont_int = varargin{3};
	if isempty(cont_int) cont_int = -999; end
end
if nargin >= 4;
	web_dis = varargin{4};
	if isempty(web_dis) web_dis = 1; end
end


% Some preliminaries
twait = 60;   % Recalculated later
nicwx = 1440; % minutes in a day
m_setup;      % m_setup needed as we will access data from Techsas

% Find the root
root_at_sea = ['/local/users/pstar/rpdmoc/rapid/data/exec/' MEXEC_G.MSCRIPT_CRUISE_STRING];
root_at_noc = ['/noc/users/pstar/rpdmoc/rapid/data/exec/' MEXEC_G.MSCRIPT_CRUISE_STRING];

if exist(root_at_sea,'dir')
	rootdir = root_at_sea;
elseif exist(root_at_noc,'dir')
	rootdir = root_at_noc;
else
	fprintf(1,'Where are we?')
	return
end

% To access other mfiles that are needed
addpath([rootdir '/mfiles/misc']);

% Will save plots here so do not need to replot bathymetry each time
work_dir = [rootdir '/mfiles/rapid_widgit_lite/working/'];
web_dir  = [rootdir '/mfiles/rapid_widgit_lite/web_server_home/'];

% Check Techsas data streams are oknav_stream = 'posmvpos';
nav.stream = 'posmvpos';
nav.vars = 'time lat long gndspeed  gndcourse';
test_tech_stream(nav.stream,nav.vars);
usbl.stream = 'usbl';
usbl.vars = 'time lat long alt';
test_tech_stream(usbl.stream,usbl.vars);
dep.stream = 'sim';
dep.vars = 'time depthm';
test_tech_stream(dep.stream,dep.vars);
gyr.stream = 'gyro_s';
gyr.vars = 'time heading';
test_tech_stream(gyr.stream,gyr.vars);
log_C.stream = 'log_chf';
if mtresolve_stream(log_C.stream) ~= log_C.stream
	log_C.vars = 'time speedfa speedps';
	test_tech_stream(log_C.stream,log_C.vars);
	no_log_C = 0;
else
	no_log_C = 1;
end
log_S.stream = 'log_skip';
log_S.vars = 'time forewaterspeed portwaterspeed';
test_tech_stream(log_S.stream,log_S.vars);
win.stream = 'winch';
win.vars = 'time cablout rate';
test_tech_stream(win.stream,win.vars);

% Read in info from file
fprintf(1,' ****** Check location information is correct and up to date ******** \n')
fprintf(1,' ******  -- Reading data from moor_pos.dat -- ******* \n ')
fomp = fopen('moor_pos.dat');
aabb = textscan(fomp,'%s%f%f%f');
mooring.names = aabb{1};
mooring.lat = dm2dd(aabb{2});
mooring.lon = dm2dd(aabb{3});
mooring.depth = aabb{4};

% Some other details of the plots
kmpd = 111.2; % Km per degree latitude
plsz = plot_size;    % Square plot with size in nm

% Find which mooring we are looking at
moor_list = char(mooring.names);
if isempty(this_moor)
  disp('Need to provide a location descriptor e.g."ship" or a mooring')
  return
else
  imoor = strmatch(this_moor,moor_list);
% Strmatch can give more than one match (e.g.'ebh4' will match to 'ebh4l' so will check later
end

% Either a mistake, using ship's position or have a match
if isempty(imoor)
  if strmatch(this_moor, 'ship') == 1
% Use ship's location
    disp('Using ship location')
	xtemp = mtload(nav_stream, now, now-0.01, 'time lat long');
        mlat = xtemp.lat(end);
	mlon = xtemp.long(end);
        mdepth = NaN;
%        plsz = 10;
  elseif strmatch(this_moor, 'user') == 1
% User enters location
    disp('User to enter location')
    mlat = input('Please enter decimal latitude:');
    mlon = input('Please enter decimal longitude:')
    mdepth = NaN;
    plsz = 10;  
  elseif strmatch(this_moor, 'check') == 1
% List mooring names and positions ao that they can be checked
    fprintf(1,'Checking instrument positions \n')
    nms = length(mooring.names);
    fprintf(1,'Name    Latitude  Longitude  Approx Depth  \n')
    for i = 1:nms
      fprintf(1,'%7s  %8.4f  %8.4f    %6.1f ', ...
             char(mooring.names(i)),mooring.lat(i),mooring.lon(i),mooring.depth(i));
      [latdd, latmm] = dd2dm2(mooring.lat(i));
      [londd, lonmm] = dd2dm2(mooring.lon(i));
      fprintf(1,' %4.0f %5.2f N %4.0f %5.2f E  \n',latdd,latmm,londd,lonmm)
    end
    return
  else
    disp('Sorry no matching mooring name')
    return  
  end
else
% It is possible to have two matches so check which one
  imoor = find(strcmp(this_moor,mooring.names));
% Use this mooring position
  disp(['Plotting for mooring ' this_moor])
  mdepth = mooring.depth(imoor);
  mlat = mooring.lat(imoor);
  mlon = mooring.lon(imoor);
end

% Set values for map boundary
dlat = plsz/60;
dlon = 1.4*dlat*cos(pi*mlat/180);
north = mlat + dlat/2;
south = mlat - dlat/2;
west = mlon - dlon/2;
east = mlon + dlon/2;

% **** Start mina loop now ****
for icwx = 1:nicwx
%  t0 = now+1.3/1440;  % Clock on Oceanuson D382 was slow by 90s 
%  t0 = now;           % then the time was corrected
%  and then it becamse fast!
% SO better to see last data written to Techsas
  [ddd,uuu] = mtlast(nav.stream);
% Allow 10s lag as all stream don't update same time
  t0 = ddd.time + MEXEC_G.uway_torg-10/86400;

% Is there a plot file already?
  if strmatch(this_moor, 'ship') == 1
    xtemp = mtload(nav.stream, now, now-0.012, 'time lat long');
    mlat = xtemp.lat(end);
    mlon = xtemp.long(end);
    mdepth = NaN;
    north = mlat + dlat/2;
    south = mlat - dlat/2;
    west = mlon - 0.5*dlon/2;
    east = mlon + 1.5*dlon/2;
    rapid_usbl_map(web_dis,work_dir,this_moor,mlat,mlon,mdepth,north,south,east,west,cont_int);
  elseif exist([work_dir this_moor '_map.fig'],'file') & icwx > 1
    open([work_dir this_moor '_map.fig']);
  else
    rapid_usbl_map(web_dis,work_dir,this_moor,mlat,mlon,mdepth,north,south,east,west,cont_int);
  end

  m_proj('lambert','lon',[west, east],'lat',[south, north]);

% Neeed to get some navigation data etc.  Get GPS first
  tfreq = 300;
  tback = 0.0834;
  xdata = mtload(nav.stream, t0-tback, t0, nav.vars);
% What time is it
  tX = xdata.time(end) + MEXEC_G.uway_torg;
% and where are we
  latdeg = floor(mlat);
  londeg = floor(-mlon);
  latmin = 60*(mlat-latdeg);
  lonmin = 60*(-mlon-londeg);
  postx = sprintf(': %i %5.2f N   %i %5.2f W ',latdeg,latmin,londeg,lonmin);
  pltitxt = sprintf(' %s%s \n %s ',this_moor,postx,datestr(tX));
  disp(['Starting new plot now ' pltitxt])
  
% Add details to existing plot
  subplot(1,2,1)
  hold on;
  title(pltitxt,'FontSize',14);
  tmx =  downsample(xdata.time,tfreq);
  tme  =  MEXEC_G.uway_torg + tmx;
  lats =  downsample(xdata.lat,tfreq);
  lons =  downsample(xdata.long,tfreq);
  gspd =  downsample(xdata.gndspeed,tfreq);
  gcse =  downsample(xdata.gndcourse,tfreq);

  if ~no_log_C
	  ydata = mtload(log_C.stream, t0-tback, t0, log_C.vars);
	  wsfa =  ndinterp1(ydata,'speedfa',tmx);
	  wsps =  ndinterp1(ydata,'speedps',tmx);
  else
	 ycdata = mtload(log_S.stream, t0-tback, t0, log_S.vars);
	 wsfa =  ndinterp1(ycdata,'forewaterspeed',tmx);
	 wsps =  ndinterp1(ycdata,'portwaterspeed',tmx);
  end
 
  zdata = mtload(dep.stream, t0-tback, t0, dep.vars);
  wdep =  ndinterp1(zdata,'depthm',tmx);    % Some NaN because tries to exrapolate

  adata = mtload(gyr.stream, t0-tback, t0, gyr.vars);
  gyro =  ndinterp1(adata,'heading',tmx);

  wdata = mtload(win.stream, t0-tback, t0, win.vars);
  win_cab =  ndinterp1(wdata,'cablout',tmx); win_cab = win_cab(end);
  win_rat =  ndinterp1(wdata,'rate',tmx); win_rat = win_rat(end);
  
% USBL data
  udata = mtload(usbl.stream, t0-tback, t0, usbl.vars);
  lons2 = ndinterp1(udata,'long',tmx);
  lats2 = ndinterp1(udata,'lat',tmx);
  zzi = ndinterp1(udata,'alt',tmx);
  for i = 1:length(tmx)-1
	  iuse = udata.time > tmx(i) & udata.time < tmx(i+1);
	  lons2(i) = nanmedian(udata.long(iuse));
	  lats2(i) = nanmedian(udata.lat(iuse));	
  end  
% Need to do some averageing
  navsec = 0.5*twait;
  navsec = min(navsec,600);  % maximum 10 minutes
  navsec = max(navsec,60);    % Minimum one minute 
  fprintf(1,'Averageing ship speed for %5.1f seconds \n',navsec)
  sog = nanmean(xdata.gndspeed(end-navsec:end));
  cog = dirav(xdata.gndcourse(end-navsec:end));
  dagyro = dirav(adata.heading(end-navsec:end));
 
  x1 = lons(end);
  y1 = lats(end);
  
  ixb = find(~isnan(udata.long));
  x2 = nanmedian(udata.long(ixb(end-18):ixb(end)));
  y2 = nanmedian(udata.lat(ixb(end-18):ixb(end)));

  ixz = find(~isnan(zzi)) ; 
  z2 = zzi(ixz(end));
  dz2 = (zzi(ixz(end-1))-zzi(ixz(end)))/(ixz(end)-ixz(end-1));

  if strmatch(this_moor,'ship') == 1
    etakm = NaN;
    etam = NaN;
    etaNM = NaN;
    etaGMT = NaN;
  else
    etakm = m_lldist([x1,mlon],[y1,mlat]);
    kmX = m_lldist([x1,mlon],[mlat,mlat]);
    kmY = m_lldist([mlon,mlon],[y1,mlat]);
    dir2m = atan2(sign(mlon-x1)*kmX,sign(mlat-y1)*kmY);  % A fudge ok for small disances
	
    etakm2 = m_lldist([x2,mlon],[y2,mlat]);
    kmX2 = m_lldist([x2,mlon],[mlat,mlat]);
    kmY2 = m_lldist([mlon,mlon],[y2,mlat]);
    dir2m2 = atan2(sign(mlon-x2)*kmX2,sign(mlat-y2)*kmY2);  % A fudge ok for small disances
	dir2m2 = dir2m2*180/acos(-1);

    kmX3 = m_lldist([x2,x1],[y1,y1]);
    kmY3 = m_lldist([x1,x1],[y2,y1]);
    dir2m3 = atan2(sign(x1-x2)*kmX3,sign(y1-y2)*kmY3);  % A fudge ok for small disances
	dir2m3 = dir2m3*180/acos(-1);
	
    ho_rn = m_lldist([x2,x1],[y2,y1]);	
	sl_rn = sqrt(1e6*ho_rn^2+z2^2);
    etaNM = etakm*60/111.2;
    etam = 60*etaNM/sog;
%   etam = etam*cos(dir2m-pi*cog/180);  
    time_now=now;
    etaGMT=datestr(time_now+etam/60/24,'HH:MM   dd/mm/yy');
  end

% Scale vectors so in reasonable size range
  vec_scl = (0.005+0.05/sog);
  dx = vec_scl * (1/cos(lats(1)*pi/180))*sog*sin(cog*pi/180)/60;
  dy = vec_scl * sog*cos(cog*pi/180)/60;
 
  xx2 = x1+dx;
  yy2 = y1+dy;
 
  wspd = sqrt(wsfa.^2+wsps.^2);
  wdir = atan2(wsps,wsfa)+gyro*pi/180;

  if sum(~isnan(wdep)) >= 1
%   corr_depth = carter(xdata.lat(end),xdata.long(end),zdata.depthm(end));
    corr_struct = mcarter(lats(end),lons(end),wdep(end));
    corr_depth = corr_struct.cordep;
  else
    corr_depth = NaN;
  end

% Scale vectors so in reasonable size range for the log
vec_scl = (0.2+0.25/nanmean(wspd(end-1:end)));
dx2 = vec_scl * (1/cos(lats(1)*pi/180))*wspd.*sin(wdir)/60;
dy2 = vec_scl * wspd.*cos(wdir)/60;
  
% Use gyro for shipas heading
%vec_scl = (east-west)/20;
%dx2 = vec_scl * (1/cos(lats(1)*pi/180))*sin(dgy*pi/180)/60;
%dy2 = vec_scl * cos(dgy*pi/180)/60;

  x3 = x1+nanmean(dx2(end-1:end));
  y3 = y1+nanmean(dy2(end-1:end));
  
  % Just gyro heading
  vec_scl = (east-west)/20;
  dx3 = vec_scl * (1/cos(lats(1)*pi/180))*sin(dagyro*pi/180);
  dy3 = vec_scl * cos(dagyro*pi/180);
  
  x4 = x1 + dx3;
  y4 = y1 + dy3;
  
  subplot(1,2,1)
  h1  =  m_plot(lons,lats,'r+-');
  h1b =  m_plot(lons(end),lats(end),'kx','LineWidth',1,'MarkerSize',9);
  iu2u = ~isnan(lons2);
  h2 =  m_plot(lons2(iu2u),lats2(iu2u),'b+-');
%  h2b =  m_plot(lons2(end),lats2(end),'kx','LineWidth',1,'MarkerSize',9);
  h2b =  m_plot(x2,y2,'kx','LineWidth',1,'MarkerSize',9);
  h2 =  m_plot([x1 xx2], [y1 yy2],'g-','LineWidth',1.5);
%  h3 =  m_plot([x1 x1+dx2], [y1 y1+dy2],'k-','LineWidth',1.5);
  h4 =  m_plot([x1 x4], [y1 y4],'k-','LineWidth',1.5);
  
  subplot(1,2,2)
  set(gca,'XTickLabel',[])
  set(gca,'YTickLabel',[])
  xx00 = 0.0;
  str1 = sprintf('Speed O.G.  %5.1f knots  ',sog);
  ha = text(xx00,0.925,str1);
  set(ha,'BackgroundColor',[.75 .9 .75],'FontSize',14)
  if cog < 0
    cog = cog+360;
  end
  str6 = sprintf('Course O.G.  %5.1f deg  ',cog);
  hb = text(xx00,0.86,str6);
  set(hb,'BackgroundColor',[.75 .9 .75],'FontSize',14);

  agyro = adata.heading(end);
  str8 = sprintf('Gyro heading %6.1f deg. ',agyro);
  he = text(xx00,0.795,str8);
  set(he,'BackgroundColor',[.75 .9 .75],'FontSize',14)

  str4 = sprintf('Bearing beacon to target %4.0f deg',dir2m2);
  he = text(xx00,0.7,str4);
  set(he,'BackgroundColor',[.85 .85 .85],'FontSize',14)

  str4 = sprintf('Bearing beacon to ship   %4.0f deg',dir2m3);
  he = text(xx00,0.63,str4);
  set(he,'BackgroundColor',[.85 .85 .85],'FontSize',14)
    
  str7 = sprintf('Approx dist. ship to waypt. %5.0f m ',1000*etakm);
  hg = text(xx00,0.54,str7);
  set(hg,'BackgroundColor',[.9 .8 .8],'FontSize',14);

  str17 = sprintf('Approx dist. bcn 2 waypt.  %5.0f m ',1000*etakm2);
  hg = text(xx00,0.47,str17);
  set(hg,'BackgroundColor',[.9 .8 .8],'FontSize',14);

  str6 = sprintf('Slant range ship 2 bcn %5.0f m ',sl_rn);
  hb = text(xx00,0.40,str6);
  set(hb,'BackgroundColor',[.75 .9 .75],'FontSize',14);
%
str66 = sprintf('Winch cable out    %4.0f m ',win_cab);
hb = text(xx00,0.33,str66);
set(hb,'BackgroundColor',[.75 .9 .65],'FontSize',14);

str68 = sprintf('Winch cable rate   %5.1f m/min ',win_rat);
hb = text(xx00,0.26,str68);
set(hb,'BackgroundColor',[.75 .9 .65],'FontSize',14);
 
  str7 = sprintf('Depth of beacon   %5.0f m ',-z2);
  hg = text(xx00,0.18,str7);
  set(hg,'BackgroundColor',[.8 .8 .9],'FontSize',14);

  str17 = sprintf('Rate of descent  %6.1f m/min ',dz2/5);
  hg = text(xx00,0.1,str17);
  set(hg,'BackgroundColor',[.8 .8 .9],'FontSize',14);
  
  str17 = sprintf('Corr depth. %7.0f m ',corr_depth);
  hg = text(xx00,0.02,str17);
  set(hg,'BackgroundColor',[.8 .8 .8],'FontSize',14);
  



% Print a plot
  set(gca,'Visible','off')
  set(gcf,'PaperUnits','centimeters','PaperPosition',[-2 0 25 16 ])
  if web_dis == 1
  	 eval(['print -dpng ' work_dir 'latest_plot']);
     close(gcf);
     eval(['!\cp  ' work_dir 'latest_plot.png ' web_dir]);
    twait = 10;
	
    fprintf(1,'Waiting for %i seconds before running agian \n',floor(twait))
    pause(twait);
    fprintf(1,' \n')
% *** End main loop ***
  else
	return
  end
  end
