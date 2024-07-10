function Anchor_seabed_triang(varargin)
%
% Anchor_seabed_triang(varargin)
% Anchor_seabed_triang(moor_loc, release_height, transducer_depth)
%
% e.g.
% >> Anchor_seabed_triang('ebh3', 1, 5)
% or without input arguments,
% >> Anchor_seabed_triang
% will prompt user for inputs
%
% Find the anchor position from 3 or more positions and ranges from text
% file
%   {MOORPROC_G.moordatadir}/raw/{cruise}/moor_positions/{moor}_times.txt
%
% Format of text file should be
%    YYYY MM DD HH MM SS range depth(if applicable)
%
% First line is anchor drop:
%   - enter zero for range
%   - for depth, a POSITIVE value will be used as the estimated
%   (uncorrected) depth; to instead interpolate from the ship's underway
%   depth stream to the anchor drop time, enter a NEGATIVE value on this
%   line
%
% Following lines are triangulation times (enter zero for depth)
%
% e.g.
%
% ebh4l6_times.txt:
% 2015 10 28 10 28 20 0 -1014
% 2015 10 28 10 46 10 1157 0
% 2015 10 28 11 00 30 1549 0
% 2015 10 28 11 11 30 1547 0
%
% Also solves for triangulation position.
%
% Saves results to moor_positions/{moor}_triangle.txt
%
% This version uses calls to mrvdas/mtechsas/mscs to get data from ship
% underway system, and uses m_map toolbox to use m_fdist and m_idist to
% calculate distances on spheroid.
%
% DAS October 2012 Updated for JC103 May 2014
%

global MEXEC_G MOORPROC_G

mcruise = MOORPROC_G.cruise;
opt1 = 'ship'; opt2 = 'datasys_best'; get_cropt %default_navstream
nav_stream = default_navstream;
dirout = fullfile(MOORPROC_G.moordatadir,'raw',mcruise,'moor_positions');
iscor = 0;

% User input
if nargin>0
    loc_name = varargin{1};
    if nargin>1
        rht = varargin{2};
        if nargin>2
            td = varargin{3};
        end
    end
else
    fprintf(1,'\n Enter mooring name (e.g. ebh3)). Times will then be read from <name>_times.txt: \n');
    fprintf(1,' - in the folder %s \n',dirout)
    fprintf(1,' - with one range per line in format  YYYY MM DD HH MM SS range. \n')
    fprintf(1, ' Output will be saved to <name>_triangle.txt in the same directory. \n\n');
    loc_name = input('Enter mooring name (e.g. ebh3): ','s');
end
if ~exist('rht','var')
    rht=input('Enter approximate height of release above seabed in metres: ');
end

if ~isempty(loc_name)
    if ~exist(dirout,'dir')
        warning('making directory %s', dirout)
        mkdir(dirout)
    end
    filein = fullfile(dirout, [loc_name '_times.txt']);
    fileout = fullfile(dirout, [loc_name '_triangle.txt']);
    disp(['Reading times from ' filein]);
    if exist(filein,'file')
        indata = load(filein);
    else
        error('input file %s does not exist',filein)
    end
    disp(['Opening output file ' fileout]);
    disp('If file exists contents will be overwritten')
    iout = fopen(fileout,'w+');
end

% Read in data from text file assume first point is anchor drop
tme = datenum(indata(:,1:6));
range = indata(:,7); range(1) = -999;
if size(indata,2) == 8
    dptvl = indata(1,8);
else
    dptvl = -999;
end

% Now get data from underway data system
switch MEXEC_G.Mshipdatasystem
    case 'techsas'
        dep_stream = 'sim';
        pos = mtload(nav_stream,datevec(tme(1)-0.04), ...
            datevec(tme(end)+0.01),'time lat long');
        pos.time = pos.time + MEXEC_G.uway_torg;
        pos.lon = pos.long;
        if dptvl<0
            dep = mtload(dep_stream,datevec(tme(1)-0.01), ...
                datevec(tme(1)+0.01),'time depthm');
            dep.time = dep.time + MEXEC_G.uway_torg;
        end
    case 'rvdas'
        dep_stream = 'ea640_sddbs';
        dep_vars = 'time waterdepth';
        pos = mrload(nav_stream,...
            datevec(tme(1)-0.04),datevec(tme(end)+0.01),...
            'time latitude longitude');
        pos.lat = pos.latitude;
        pos.lon = pos.longitude;
        pos.time = pos.dnum;
        if dptvl<0
            dep = mrload(dep_stream,datevec(tme(1)-0.01), ...
                datevec(tme(1)+0.01),'time waterdepthmeters');
        end
    case 'scs'
        dep_stream = 'singleb';
        dep_vars = 'time waterdepth';
        pos = msload(nav_stream,...
            datevec(tme(1)-0.04),datevec(tme(end)+0.01),...
            'time GPS-Furuno-GGA-lat GPS-Furuno-GGA-lon');
        pos.lat = pos.GPS_Furuno_GGA_lat;
        pos.lon = pos.GPS_Furuno_GGA_lon;
        [pos.time, IA, ~] = unique(pos.time);
        pos.lat = pos.lat(IA);
        pos.lon = pos.lon(IA);
        if dptvl<0
            dep = msload(dep_stream,datevec(tme(1)-0.01), ...
                datevec(tme(1)+0.01),dep_vars);
            dep.depthm = dep.waterdepth;
            [dep.time, IA, ~] = unique(dep.time);
            dep.depthm = dep.depthm(IA);
        end
end
iscor = 0; %depth is from singlebeam
if ~exist('td','var')
    td =input(sprintf('Enter approximate transducer depth correction to %s, %s in metres: ',dep_stream,dep_vars(end)));
end

% interpolate triangulation positions
lat = interp1(pos.time,pos.lat,tme);
lon = interp1(pos.time,pos.lon,tme);

% information on depth
if dptvl < 0
    % Need to check the depths
    iabs = dep.depthm == 0;
    if sum(~iabs)>2
        wd2 = dinterp1(dep.time(~iabs),dep.depthm(~iabs),tme);
        f1 = figure;
        plot(dep.time(~iabs),dep.depthm(~iabs));
        hold on; grid on;
        y1 = ylim; t1 = tme(1);
        plot([t1 t1],y1,'r');grid on;
        datetick;xlabel('Time'),ylabel('Depth');
    else wd2 = 0;
    end
    fprintf(1,'Uncorrected water depth at anchor drop %5.1f m \n',wd2(1));
    ich = input('Enter 1 to accept or 2 to change values: ');
    if ich == 2
        wd = input('Enter uncorrected water depth in metres: ');
    else
        wd = wd2(1);
    end
    %   close(f1)
    % The corrected water depth
    corr_struct = 	mcarter(lat(1),lon(1),wd);
    wd_corr = corr_struct.cordep;
else
    corr_struct = 	mcarter(lat(1),lon(1),dptvl);
    wd_corr = corr_struct.cordep;
    wd = 2*dptvl-wd_corr;
    if iscor
        wd_corr = dptvl; 
    end
end
%add td here if using waterdepthmetrefromtransducer rather than
%waterdepthmetrefrom surface***


% Separate the anchor drop position
latD = lat(1);
lonD = lon(1);
tmeD = tme(1);
range = range(2:end);
lat = lat(2:end);
lon = lon(2:end);
tme = tme(2:end);
no_fixes=length(range);
% correct the ranges for speed of sound as well (use anchor drop position
% so all are adjusted the same in case of crossing a carter table boundary!)
range_corr = range;
if max(range)<8000
    corr_struct = 	mcarter(lat(1),lon(1),range);
    range_corr = corr_struct.cordep;
end


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

% Work out effective horizontal range (wd_corr is already relative to
% surface, only need to correct for release height)
rangeh = sqrt(range_corr.^2 - (wd_corr-rht)^2);

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
titletext2=sprintf( ...
    'Corrected water depth: %5.0f m. Release Height: %3.0f m. Transducer depth: %3.0f', ...
    wd_corr,rht,td);
title({titletext1;titletext2});
xlabel('Longitude'); ylabel('Latitude')


% Calculate triangulation point by solving least squares problem
% Nb take differences of equations to make linear

%lon0 = mean(lon);
%lat0 = mean(lat);
lat0 = latD;
lon0 = lonD;

[long1,latg1,~,~,reser0] = solve_anchor(lon0,lat0,lon,lat,rangeh,no_fixes);
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
    [dis,az12,~] = m_idist(lon0,lat0,APlon,APlat);
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
fallback1 = sw_dist([APlat lat(1)],[APlon lon(1)],'nm');
fallback2 = sw_dist([APlat lat(1)],[APlon lon(1)],'km');

plot(APlon,APlat,'r+');
plot(APlon,APlat,'ro');
plot(lonD,latD,'g+');
plot(lonD,latD,'go');
titletext3=['Red = anchor seabed position. ',num2str(APlat),'N ',num2str(APlon),'W.'];
titletext4=['Fall back = ' num2str(fallback2) ' km = ' num2str(fallback1) ' nm.'];
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
    plotout = fullfile(dirout, [loc_name '_triangle']);
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
if strcmp(topo_map,'none')
    swtx = sprintf('No swath data available depth used was %7.1f \n',wd_corr);
else
    rdpath = '/local/users/pstar/projects/rpdmoc/bathym_data/from_cruises/';
    load([rdpath topo_map]);
    newdep = interp2(lon,lat,depth,APlon,APlat);
    swtx = sprintf('Depth used was %7.1f but corrected depth from swath map at\n trilaterated position is %7.1f;\n you may want to edit input file and rerun\n',wd_corr,newdep);
end
fprintf(1,'%s',swtx)

% Check what bathymetry was when ship passed over the anchor point
ixt = pos.time > tme(1)-0.5 & pos.time < tme(1) +0.5;
tme_sh = pos.time(ixt);
dis_sh = 1000*111.2*sqrt( (pos.lat(ixt)-APlat).^2+cos(APlat*pi/180)^2*(pos.lon(ixt)-APlon).^2);
im = find(dis_sh == min(dis_sh)); im = im(1);
tp_sh = tme_sh(im);

if dptvl < 0 && sum(~iabs) > 2
    wd_sh = dinterp1(dep.time(~iabs),dep.depthm(~iabs),tp_sh);
    corr_struct = mcarter(lat(1),lon(1),wd_sh);
    wd_corr_sh = corr_struct.cordep;
else
    corr_struct = mcarter(lat(1),lon(1),wd);
    wd_corr_sh = corr_struct.cordep;
end

fprintf(1,'\n Closest ship track to anchor was %6.1f m \n',dis_sh(im))
fprintf(1,' at time %s \n',datestr(tp_sh))
fprintf(1,' Where depth was %6.1f m \n',wd_corr_sh)


% Now data that we save in the output file
fprintf(iout,'%s  %s \n',loc_name,datestr(tmeD,1));
fprintf(iout,'Anchor drop at: %s %8.4f %8.4f Corr. water depth: %6.1f \n', ...
    datestr(tmeD,31),latD,lonD,wd_corr);
fprintf(iout,'%s',swtx)
fprintf(iout,'Date      Time       Lat   Lon  Slant range  Horiz range Residual \n');
for i = 1:no_fixes
    fprintf(iout,'%s  %8.4f %8.4f %7.0f %7.0f %7.0f \n', ...
        datestr(tme(i),31),lat(i),lon(i),range(i),rangeh(i),reser(i));
end
fprintf(iout,'Trilaterated position \n');
fprintf(iout,'%8.4f %8.4f \n',APlat,APlon);
fprintf(iout,'%s \n',title4);
fprintf(iout,'Fallback  %s m \n',fallback);
fprintf(iout,'Comments: %s \n',us_com);


