function rapid_widgit_lite(varargin)
% This script creates plot showing navigation information.
% *********** PLEASE CHECK POSITIONS ARE CORRECT *************
% The positions are stored in the file:
%  - moor_pos.dat
% See README file for instuctions on displaying the plot on the web.
% Copies final plots to the local directory web_server_home which should
% be linked to the correct folder (/var/www/html on Koaekea)
%
% Normally plot is to file only and there is no screen plot.
% Updates when advance 5% of distance to waypoint (minimum 10s)
% The script will continue for one day or until user interupts
% Primarily for use during deployment of moorings.
%
% There are a number of optional inputs:
%   - this_moor - is a string and can have values:
%   	- 'ship'  Will create plot cetred on ship's position
%		-  any of the mooring names listed in moor_pos.dat
% 		- 'check' Will list all mooring names, positions and nominal depths
%   - falback - fallback in metres.  The distance by which the ship will overrun the target
%   - plot_size - is a real number and is size of plot in nautical miles
%   - cont_int - is contour interval usually do not need to set this
%   - web_dis - default is 1 to copy files to webserver and have no matlab display
%               but if 0 then will display in Matlab only and will plot jsut once
%
% Calls following scripts:  rapid_web_map.m, dirav.m, dm2dd.m, dd2dm2.m dinterp1.m
%
% Nb the calculation of ETA assumes that hte ship is heading directly
% towards the waypoint
%
% Update DAS October 2018 on JC174
% When using matlab 2015 need this
warning('off','snctools:nc_varget:fillValueMismatch')
if exist('m_common.m','file')
    m_common
else
    m_setup
end

% Some preliminaries
twait = 60;    % Recalculated later - time between plots
nicwx = 1440;  % minutes in a day
kmpd  = 111.2; % km per degree latitude

root_at_sea = ['/local/users/pstar/projects/rpdmoc/rapid/data/exec/' MEXEC_G.MSCRIPT_CRUISE_STRING];
root_at_noc = ['/noc/users/pstar/rpdmoc/rapid/data/exec/' MEXEC_G.MSCRIPT_CRUISE_STRING];

if exist(root_at_sea,'dir')
    rootdir = root_at_sea;
elseif exist(root_at_noc,'dir')
    rootdir = root_at_noc;
else
    fprintf(1,'Where are we?')
    return
end

addpath([rootdir '/mfiles/misc']);

% Now look for positions in the data directory
pos_direct = [fileparts(fileparts(rootdir)) '/moor/raw/' MEXEC_G.MSCRIPT_CRUISE_STRING '/moor_positions/tables/'];
moor_targets = [pos_direct 'moor_pos_target.dat'];
moor_inwater = [pos_direct 'moor_pos_pre_cruise.dat'];

% Make sure techsas streams are linked
switch MEXEC_G.Mshipdatasystem
    case 'techsas'
		[s,w] = unix('techsas_linkscript');
		if s ~- 0
    		disp(w)
		end
	end

% Check Techsas data streams are ok
switch MEXEC_G.Mshipdatasystem
    case 'techsas'
		t_origin = MEXEC_G.uway_torg;
        nav.stream = 'posmvpos';
        nav.vars = 'time lat long gndspeed  gndcourse';
        dep.stream = 'sim';
        dep.vars = 'time depthm';
        gyr.stream = 'gyro_s';
        gyr.vars = 'time heading';
        log_S.stream = 'log_skip';
        log_S.vars = 'time forewaterspeed portwaterspeed';
    case 'rvdas'
        nav.stream = 'pospmv';
        nav.vars = 'time latitude longitude altitude';
        % nav streams separated in rvdas so need to load 2 streams and
        % merger later
        nav2.stream = 'vtgcnav';
        nav2.vars = 'time courseoverground speedknots';
        dep.stream = 'singleb';
        dep.vars = 'time waterdepth';
        gyr.stream = 'hdtgyro';
        gyr.vars = 'time heading';
        log_S.stream = 'logskip';
        log_S.vars = 'time speed_forward_raw speed_stbd_raw';
	case 'scs'
		ms_raw_to_sed(57);
		update_allmat;
		t_origin = 0;
        nav.stream = 'posfur';
        nav.vars   = 'time GPS-Furuno-GGA-lat GPS-Furuno-GGA-lon';
        dep.stream = 'knudsen';
        dep.vars   = 'time SingleBeam-Knudsen-PKEL99-uncdepth';
        gyr.stream = 'gyro';
        gyr.vars   = 'time Gyro1-HDT-heading';
        log_S.stream = 'dopplerlog';
        log_S.vars = 'time SpeedLog-Furuno-VBW-speedfa SpeedLog-Furuno-VBW-speedps';
end
test_tech_stream(nav.stream,nav.vars);
test_tech_stream(dep.stream,dep.vars);
test_tech_stream(gyr.stream,gyr.vars);
test_tech_stream(log_S.stream,log_S.vars);
%% log_C.stream = 'log_chf';
%% log_C.vars = 'time speedfa speedps';
%% try
%%     test_tech_stream(log_C.stream,log_C.vars);
%%     no_log_C = 0;
%% catch
%%     no_log_C = 1;
%% end
 no_log_C = 1;
 
% Default values of input variables
this_moor = 'check';
fall_backm = 0;      % Fallback in metres
plot_size  = 8;		% Default choice
cont_int = -999;   	% < 0 means automatic choice
web_dis = 1;

% Optional values from user
if nargin >= 1;
	this_moor = char(varargin{1});
	if isempty(this_moor) this_moor = 'ship'; end
end
if nargin >= 2;
	fall_backm = varargin{2};
	if isempty(fall_backm) fall_backm = 0; end
end
fprintf(1,'\n *** Fallback = %5.0fm *** \n\n',fall_backm)
if nargin >= 3;
	plot_size = varargin{3};
	if isempty(plot_size) plot_size = 8; end
end
if nargin >= 4;
	cont_int = varargin{4};
	if isempty(cont_int) cont_int = -999; end
end
if nargin >= 5;
	web_dis = varargin{5};
	if isempty(web_dis) web_dis = 1; end
end

% Will save plots here so do not need to replot bathymetry each time
work_dir = [rootdir '/mfiles/rapid_widgit_lite/working/'];
web_dir  = [rootdir '/mfiles/rapid_widgit_lite/web_server_home/'];

% Read in info from file
fprintf(1,' ****** Check location information is correct and up to date ******** \n')
fprintf(1,' ****** >> Reading data from  %s ******* \n ',moor_targets)
fprintf(1,' ******     %s     ******')

% Open file with target positions
fomp = fopen(moor_targets);
aabb = textscan(fomp,'%s%f%f%f');
mooring.names = aabb{1};
mooring.lat = dm2dd(aabb{2});
mooring.lon = dm2dd(aabb{3});
mooring.depth = aabb{4};
fclose(fomp);

% Some other details of the plots
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
switch MEXEC_G.Mshipdatasystem
case 'techsas'
xtemp = mtlast(nav.stream);
case 'rvdas'
xtemp = mrlast(nav.stream);
case 'scs'
xtemp = mslast(nav.stream);
xtemp.lat  =xtemp.GPS_Furuno_GGA_lat;
xtemp.long =xtemp.GPS_Furuno_GGA_lon;
end
    mlat = xtemp.lat;
	mlon = xtemp.long;
    mdepth = NaN;
  elseif strmatch(this_moor, 'user') == 1
% User enters location
    disp('User to enter location')
    mlat = input('Please enter decimal latitude:');
    mlon = input('Please enter decimal longitude:');
    mdepth = NaN;
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

% Now read in positions of all moorings in the water
fomp = fopen(moor_targets);
aabb = textscan(fomp,'%s%f%f%f');
mooring.names = aabb{1};
mooring.lat = dm2dd(aabb{2});
mooring.lon = dm2dd(aabb{3});
mooring.depth = aabb{4};

% Set values for map boundary
dlat  = plsz/60;
dlon  = dlat*cos(pi*mlat/180);
north = mlat + dlat/2;
south = mlat - dlat/2;
west  = mlon - dlon/2;
east  = mlon + dlon/2;

% **** Start main loop now ****
for icwx = 1:nicwx
switch MEXEC_G.Mshipdatasystem
case 'techsas'
  [ddd,uuu] = mtlast(nav.stream);
case 'rvdas'
[ddd,nnn,uuu] = mrlast(nav.stream);
    ddd.time=ddd.utctime;
case 'scs'
  [ddd] = mslast(nav.stream);
end
% Allow 10s lag as all stream don't update same time
  t0 = ddd.time + t_origin-10/86400;

% Is there a plot file already?
  if strmatch(this_moor, 'ship') == 1
	  switch MEXEC_G.Mshipdatasystem
	  case 'techsas'
	  xtemp = mtlast(nav.stream);
	  case 'rvdas'
	  xtemp = mrlast(nav.stream);
	  case 'scs'
	  xtemp = mslast(nav.stream);
	  xtemp.lat  = xtemp.GPS_Furuno_GGA_lat;
	  xtemp.long = xtemp.GPS_Furuno_GGA_lon;
	  end
    mlat = xtemp.lat(end);
    mlon = xtemp.long(end);
    mdepth = NaN;
    north = mlat + dlat/2;
    south = mlat - dlat/2;
    west = mlon - dlon/2;
    east = mlon + dlon/2;
    rapid_web_map(web_dis,work_dir,this_moor,mlat,mlon,mdepth,north,south,east,west,cont_int);
  elseif exist([work_dir this_moor '_map.fig'],'file') & icwx > 1
    open([work_dir this_moor '_map.fig']);
  else
    rapid_web_map(web_dis,work_dir,this_moor,mlat,mlon,mdepth,north,south,east,west,cont_int);
  end
  
% Plot location of mooring
  m_proj('lambert','lon',[west, east],'lat',[south, north]);
  imp = mooring.lon >= west & mooring.lon <= east & mooring.lat >= south & mooring.lat <= north;
  if ~isempty(imp)
  	m_plot(mooring.lon(imp),mooring.lat(imp),'+','Color',[0.5 0.5 0.5],'MarkerSize',12,'LineWidth',1)
  	m_plot(mlon,mlat,'k+','MarkerSize',12,'LineWidth',1.5)
	if exist('drlat','var')
  	   m_plot(drlon,drlat,'mx','MarkerSize',12,'LineWidth',2)
    end
  end

% Neeed to get some navigation data etc.  Get GPS first
  tfreq = 300;
  switch MEXEC_G.Mshipdatasystem
    case 'techsas'
        xdata = mtload(nav.stream, t0-0.0417, t0, nav.vars);
    case 'rvdas'
        xdata = mrload(nav.stream, t0-0.0417, t0, nav.vars);
        xdata.time=xdata.dnum;
        xdata.lat=xdata.latitude;
        xdata.long=xdata.longitude;
        % load 2nd nav stream for speed/course
        xdata2 = mrload(nav2.stream, t0-0.0417, t0, nav2.vars);
        % interpolate onto same timegrid
        xdata.gndspeed=interp1(xdata2.dnum,xdata2.speedknots,xdata.time);
        xdata.gndcourse=interp1(xdata2.dnum,xdata2.course,xdata.time);
    case 'scs'
        xdata = msload(nav.stream, t0-0.0417, t0, nav.vars);
        xdata.lat  = xdata.GPS_Furuno_GGA_lat;
        xdata.long = xdata.GPS_Furuno_GGA_lon;
        xdata.gndspeed = NaN*xdata.lat;
        xdata.gndcourse= NaN*xdata.lat;
  end

% What time is it
  tX = xdata.time(end) + t_origin;
% and where are we
  latdeg = floor(mlat);
  londeg = floor(-mlon);
  latmin = 60*(mlat-latdeg);
  lonmin = 60*(-mlon-londeg);
  postx = sprintf(': %i %5.2f N   %i %5.2f W ',latdeg,latmin,londeg,lonmin);
  if fall_backm > 0
	  fbtxt = sprintf(', F.B.=%4im',fall_backm);
  else
	  fbtxt = '';
  end
  pltitxt = sprintf(' %s%s \n %s%s ',this_moor,postx,datestr(tX),fbtxt);
  disp(['Starting new plot now ' pltitxt])
  
% Add details to existing plot
  subplot(1,2,1);
  hold on;
  title(pltitxt,'FontSize',14);
  tmx =  downsample(xdata.time,tfreq);
  tme  =  t_origin + tmx;
  lats =  downsample(xdata.lat,tfreq);
  lons =  downsample(xdata.long,tfreq);
  gspd =  downsample(xdata.gndspeed,tfreq);
  gcse =  downsample(xdata.gndcourse,tfreq);

% depth data
  switch MEXEC_G.Mshipdatasystem
    case 'techsas'
        zdata = mtload(dep.stream, t0-0.0417, t0,dep.vars);
    case 'rvdas'
        zdata = mrload(dep.stream, t0-0.0417, t0,dep.vars);
        zdata.time=zdata.dnum;
        zdata.depthm=zdata.waterdepth;
    case 'scs'
        zdata = msload(dep.stream, t0-0.0417, t0,dep.vars);
        zdata.depthm = zdata.SingleBeam_Knudsen_PKEL99_uncdepth;
  end
  wdep = ndinterp1(zdata,'depthm',tmx);
  if isempty(zdata)
	  zdata.time = xdata.time;
	  zdata.depthm = NaN*zdata.time;
  end
  
% wsfa =  downsample(ydata.speedfa,300);
% wsps =  downsample(ydata.speedps,300);
% zdata = mtload('ea600m', t0-0.0417, t0, 'time snd');
% wdep =  downsample(zdata.snd,30);
% adata = mtload('gyro_s', t0-0.0417, t0, 'time heading');
% gyro =  downsample(adata.heading,300);

% Need to do some averageing
  navsec = 0.5*twait;
  navsec = min(navsec,600);  % maximum 10 minutes
  navsec = max(navsec,60);    % Minimum one minute 
% Gyro data
  switch MEXEC_G.Mshipdatasystem
    case 'techsas'
        gdata = mtload(gyr.stream, t0-navsec/86400, t0, gyr.vars);
    case 'rvdas'
        gdata = mrload(gyr.stream, t0-navsec/86400, t0, gyr.vars);
        gdata.time=gdata.dnum;
    case 'scs'
        gdata = msload(gyr.stream, t0-navsec/86400, t0, gyr.vars);
        gdata.heading = gdata.Gyro1_HDT_heading;
  end
  gyro =  ndinterp1(gdata,'heading',tmx);

  fprintf(1,'Averaging ship speed for %5.1f seconds \n',navsec)
  sog = nanmean(xdata.gndspeed(end-navsec:end));
  cog = dirav(xdata.gndcourse(end-navsec:end));
  if isnan(sog)
	  [sog,cog] = get_sog_cog(xdata);
  end
  dagyro = dirav(gyro);

  switch MEXEC_G.Mshipdatasystem
    case 'techsas'
        ycdata = mtload(log_S.stream, t0-navsec/86400, t0, log_S.vars);
    case 'rvdas'
        ycdata = mrload(log_S.stream, t0-navsec/86400, t0, log_S.vars);
        ycdata.forewaterspeed=ycdata.speed_forward_raw;
        ycdata.portwaterspeed=-ycdata.speed_stbd_raw; % I think this need the opposite sign
    case 'scs'
        ycdata = msload(log_S.stream, t0-navsec/86400, t0, log_S.vars);
        ycdata.forewaterspeed = ycdata.SpeedLog_Furuno_VBW_speedfa;
        ycdata.portwaterspeed = ycdata.SpeedLog_Furuno_VBW_speedps; 
  end
  wsco_f =  nanmean(ycdata.forewaterspeed);
  wsco_p =  nanmean(ycdata.portwaterspeed);
  wsco = sqrt(wsco_f.^2+wsco_p.^2);

  if ~no_log_C
    switch MEXEC_G.Mshipdatasystem
        case 'techsas'
  	        ydata = mtload(log_C.stream, t0-navsec/86400, t0, log_C.vars);
        case 'rvdas'
            ydata = mrload(log_C.stream, t0-navsec/86400, t0, log_C.vars);
    end
  	wsfa =  nanmean(ydata.speedfa);
  	wsps =  nanmean(ydata.speedps);
  else
	  wsfa = wsco_f;
	  wsps = wsco_p;
  end


  cpfa = sqrt(wsco.^2/(wsfa.^2+wsps.^2));
  wsfa = cpfa*wsfa;
  wsps = cpfa*wsps;
  wsfa = wsco_f;
  wsps = wsco_p;

  x1 = lons(end);
  y1 = lats(end);
 
  if strmatch(this_moor,'ship') == 1
    etakm = NaN;
    etam = NaN;
    etaNM = NaN;
    etaGMT = NaN;
  else
	if fall_backm > 0
	    [etakm,dir2m,dir2s] = m_idist(x1,y1,mlon,mlat);
		fprintf(1,'Distance to mooring: %6.0f m \n',etakm)
		if abs(cog-dir2m) > 90 & abs(cog+360-dir2m) > 90 & abs(cog-360-dir2m) > 90;
			dir2m = dir2m+180;
		end
		[drlon,drlat,a21] = m_fdist(mlon,mlat,dir2m,fall_backm);
% Check		[fb_ck,a12,a21] = m_idist(mlon,mlat,drlon,drlat);
		if drlon > 180
			drlon = drlon-360;
		end
		fprintf(1,'Drop %8.2f N  %8.2f W \n',dd2dm(drlat),dd2dm(drlon))
		[etakm,dir2m,dir2s] = m_idist(x1,y1,drlon,drlat);
 		fprintf(1,'Distance to drop: %6.0f m \n',etakm)

%    	check = 1000*m_lldist([mlon,drlon],[mlat,drlat])
		xlon = drlon;
		xlat = drlat;
	else
		xlon = mlon;
		xlat = mlat;
	end
 
	[etakm,dir2m,dir2s] = m_idist(x1,y1,xlon,xlat);
	etakm = etakm/1000;
    etaNM = etakm*60/111.2;
	if abs(mod(cog,360)-mod(dir2m,360)) < 30;
    	etam = 60*etaNM/sog;
	    time_now=now;
	    etaGMT=datestr(time_now+etam/60/24,'HH:MM   dd/mm/yy');
	else 
		etam = NaN;
		etaGMT='???';
		fprintf(1,'*** Not heading towards waypoint/drop no ETA ***')
	end
%   etam = etam*cos(dir2m-pi*cog/180);  
  end

% Scale vectors so in reasonable size range
  vec_scl = (0.2+0.25/sog);
  dx = vec_scl * (1/cos(lats(1)*pi/180))*sog*sin(cog*pi/180)/60;
  dy = vec_scl * sog*cos(cog*pi/180)/60;
 
  x2 = x1+dx;
  y2 = y1+dy;
 
%  wsfa
%  wsps
  wspd = sqrt(wsfa.^2+wsps.^2);
  wdir = atan2(wsps,wsfa)+dagyro*pi/180;

  if  sum(~isnan(zdata.depthm)) >= 1
      corr_struct = mcarter(xdata.lat(end),xdata.long(end),zdata.depthm(end));
      corr_depth = corr_struct.cordep;
      corr_struct = mcarter(xdata.lat(end)*ones(size(zdata.depthm)), ...
                            xdata.long(end)*ones(size(zdata.depthm)),zdata.depthm);
	  zdata.corrd = corr_struct.cordep;
  else
      corr_depth = NaN;
  end

% Scale vectors so in reasonable size range
  vec_scl = (0.2+0.4/wspd);
  dx2 = vec_scl * (1/cos(lats(1)*pi/180))*wspd.*sin(wdir)/60;
  dy2 = vec_scl * wspd.*cos(wdir)/60;
  
  x3 = x1 + dx2;
  y3 = y1 + dy2;
  
% Just gyro heading
  vec_scl = (east-west)/20;
  dx3 = vec_scl * (1/cos(lats(1)*pi/180))*sin(dagyro*pi/180);
  dy3 = vec_scl * cos(dagyro*pi/180);
  
  x4 = x1 + dx3;
  y4 = y1 + dy3;
  

  subplot(1,2,1)
  h1 =  m_plot(lons,lats,'r+-');
  h1b =  m_plot(lons(end),lats(end),'kx','LineWidth',1,'MarkerSize',9);
  h2 =  m_plot([x1 x2], [y1 y2],'g-','LineWidth',1.5);
%%  h3 =  m_plot([x1 x3], [y1 y3],'b-','LineWidth',1.5);
  h4 =  m_plot([x1 x4], [y1 y4],'k-','LineWidth',1.5);

  subplot(4,2,[2 4 6])
  set(gca,'XTickLabel',[])
  set(gca,'YTickLabel',[])
  xx00 = 0.0;
  str1 = sprintf('Speed O.G. %6.2f knots  ',sog);
  ha = text(xx00,0.925,str1);
  set(ha,'BackgroundColor',[.75 .9 .75],'FontSize',14)
  if cog < 0
    cog = cog+360;
  end
  str2 = sprintf('Course O.G.  %5.1f deg  ',cog);
  hb = text(xx00,0.825,str2);
  set(hb,'BackgroundColor',[.75 .9 .75],'FontSize',14);

  str3= sprintf('EM log speed %5.1f knots',wsco);
  hc = text(xx00,0.7,str3);
  set(hc,'BackgroundColor',[.85 .85 .85],'FontSize',14)
%  str2 = sprintf('EM log FA %5.1f knots',wsfa);
%  hc = text(xx00,0.15,str2);
%  set(hc,'BackgroundColor',[.75 .75 .9],'FontSize',14)
%  str3 = sprintf('EM log PS %5.1f knots',wsps);
  str4 = sprintf('Gyro heading %6.1f deg. ',dagyro);
  hd = text(xx00,0.6,str4);
  set(hd,'BackgroundColor',[.85 .85 .85],'FontSize',14)
  
  switch this_moor
  case {'ship'}
	  etam = 60;
  otherwise
  str4a = sprintf('Heading to drop %5.1f deg. ',dir2m);
  hd = text(xx00,0.5,str4a);
  set(hd,'BackgroundColor',[.85 .85 .85],'FontSize',14)
  
   if fall_backm > 0 
       str5 = sprintf('Approx dist. to drop %6.2f Nm ',etaNM);
   else
	   str5 = sprintf('Approx dist. to waypt. %6.2f Nm ',etaNM);
   end
  he = text(xx00,0.375,str5);
  set(he,'BackgroundColor',[.9 .8 .8],'FontSize',14);
  if abs(etam) < 150
	     if fall_backm > 0 
		     str6 = sprintf('Approx time to drop %5.0f mins',etam);
		 else
    		 str6 = sprintf('Approx time to waypt. %5.0f mins',etam);
		 end
  else
	  	if fall_backm > 0 
      	  	str6 = sprintf('Approx time to drop %6.1f hours',etam/60);
  		else
	        str6 = sprintf('Approx time to waypt. %6.1f hours',etam/60);
		end
  end 
  hf = text(xx00,0.275,str6);
  set(hf,'BackgroundColor',[.9 .8 .8],'FontSize',14)
  str7 = sprintf('ETA.     %s GMT ',etaGMT);
  hg = text(xx00,0.175,str7);
  set(hg,'BackgroundColor',[.9 .8 .8],'FontSize',14)
    fprintf(1,' %s \n %s \n', str5,str6);
end

  str8 = sprintf('Corrected depth %5.0f m',corr_depth);
  hh = text(xx00,0.05,str8);
  set(hh,'BackgroundColor',[.75 .75 .9],'FontSize',14)
  
  set(gca,'Visible','off')
	
  if sum(~isnan(zdata.depthm)) > 1
	 subplot(4,2,8)
 	 btime = zdata.time + t_origin;
	 iduse = find(now - btime < 0.5/24);
  	 bathy = zdata.corrd;
 	 idok = ~isnan(zdata.corrd(iduse)) & zdata.corrd(iduse) > 0; 
	 if sum(idok) > 10
  		bathy(iduse(~idok)) = NaN;
		bat_temp = zdata.corrd(iduse(idok));
  		bathyf = medfilt1(bat_temp - median(bat_temp),15) + median(bat_temp);
  		bathy(iduse(idok)) = bathyf;
  		plot(btime(iduse),bathy(iduse))
  		datetick('x');xlim([now-0.5/24 now]);
  		set(gca,'YDir','reverse','YAxisLocation','right');grid on
    	byls = prctile(bathy(iduse),[2.5 97.5]);
    	ylim([min(byls)-25 max(byls)+25])
		title('EA600 Corr depth - lastest 30 min')
		grid on
	end
  end

  
% Print a plot
  set(gcf,'PaperUnits','centimeters','PaperPosition',[-2 0 21 14 ])
  if web_dis == 1
  	 eval(['print -dpng ' work_dir 'latest_plot']);
     close(gcf);
%   Finally move plots to the web directory
    eval(['!\cp  ' work_dir 'latest_plot.png  ' web_dir]);
%   Now wait twait seconds before running again
    if isnan(etam)
    	  twait = 60;
    elseif etaNM < 1
    	  twait = 10;
    elseif etaNM < 2
    	  twait = 30;
    elseif 0.05*etam*60 > 90
    	  twait = 90;
    else
        twait = max(10,0.05*etam*60);
    end
    fprintf(1,'Waiting for %i seconds before running agian \n',floor(twait))
    pause(twait);
	ms_raw_to_sed(57);
	update_allmat;
    fprintf(1,' \n')
% *** End main loop ***
  else
	return
  end
  end

  function 	[sog,cog] = get_sog_cog(xdata);
	  iok = find(~isnan(xdata.lat));
	  i1 = max(1,length(iok)-60);
	  lon1 = xdata.long(iok(i1));
	  lon2 = xdata.long(iok(end));
	  lat1 = xdata.lat(iok(i1));
	  lat2 = xdata.lat(iok(end));
	  [s,a12,a21] = m_idist(lon1,lat1,lon2,lat2);  % dist in meters
	  sog = s/(xdata.time(iok(end))-xdata.time(iok(i1)));  % m per day
	  sog = sog/(1000*24);   % km/hr
	  sog = sog*60/111.2;   % Nm per hour
	  cog = a12;