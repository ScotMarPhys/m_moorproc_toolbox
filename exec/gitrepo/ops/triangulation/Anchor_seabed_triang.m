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
% following lines are triangulation times (enter zero for depth)
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
% Optional a text file called targets_{cruise}.txt can be created to list
% the target positions of the different moorings as decimal degree N and E:
% {moor} {lat} {lon}
%
% This version uses calls to mrvdas/mtechsas/mscs to get data from ship
% underway system, and uses m_map toolbox to use m_fdist and m_idist to
% calculate distances on spheroid.
%
% DAS October 2012 Updated for JC103 May 2014
%

global MEXEC_G MOORPROC_G

pd = moor_inoutpaths('reports');
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
    fprintf(1,' - in the folder %s \n',pd.trilatdir)
    fprintf(1,' - with one range per line in format  YYYY MM DD HH MM SS range. \n')
    fprintf(1, ' Output will be saved to <name>_triangle.txt in %s\n\n',pd.trilatdir);
    loc_name = input('Enter mooring name (e.g. ebh3): ','s');
end
if ~exist('rht','var')
    rht=input('Enter approximate height of release above seabed in metres: ');
end
if ~exist('td','var')
    td =input('Enter approximate transducer depth in metres: ');
end
if isempty(loc_name)
    return
end

%target positions
targ_loc=NaN(2,1);
target_fn = fullfile(pd.targetdir,['targets_',MOORPROC_G.cruise,'.txt']);
if exist(target_fn,'file')
    fileID = fopen(target_fn,'r');
    txt = textscan(fileID,'%s','delimiter','\n');
    for i=1:length(txt{1})
        if strfind(txt{1}{i},[loc_name,' '])
            targ_loc = sscanf(txt{1}{i},[loc_name,' %f %f']);
        end
    end
    fclose(fileID);
end
if numel(find(isnan(targ_loc)))>0      
    fprintf('No target position found for %s \n',loc_name)
    ta_str = input('Would you like to enter a target position? yes=1, no=0 (default=0) :');
    if ta_str
        targ_loc(1)=input('Please enter target latitude as decimal degrees N: ');
        targ_loc(2)=input('Please enter target longitude as decimal degrees E: ');
        fileID=fopen(target_fn,'w');
        fprintf(fileID,'%s %f %f \n',loc_name,targ_loc(1), targ_loc(2));
        fclose(fileID);
    end
end

filein = fullfile(pd.trilatdir, [loc_name '_times.txt']);
if ~exist(filein,'file')
    warning('input file %s not found',filein)
    return
end

%load times and ranges
disp(['Reading times from ' pd.trilatdir]);
indata = load(filein);
tme = datenum(indata(:,1:6));
range = indata(:,7);
if size(indata,2)<8
    %add depth column
    dptvl = -999+zeros(size(range));
else
    dptvl = indata(:,8);
end
dptvl = indata(:,8);

%anchor drop is first line
tmeD = tme(1); tme = tme(2:end);
rangeD = range(1); range = range(2:end);
dptvlD = dptvl(1); dptvl = dptvl(2:end);

% get ship position and depth on ship track between anchor drop and last
% ranging
[pos, dep] = get_uway([tmeD-0.1 max(tme)+0.1]);

% interpolate for positions
latD = interp1(pos.time,pos.lat,tmeD);
lonD = interp1(pos.time,pos.lon,tmeD);
lat = interp1(pos.time,pos.lat,tme);
lon = interp1(pos.time,pos.lon,tme);

%depth at anchor drop
if dptvlD>0
    if iscor %***
        %corrected depth on first line of file
        wd_corr = dptvlD;
    else
        %uncorrected depth on first line of file
        wd = dptvlD;
        corr_struct = mcarter(latD,lonD,wd);
        wd_corr = corr_struct.cordep;
    end

else
   %interpolate data from underway
   wd = dinterp1(dep.time,dep.depthm,tme);
   figure(1); clf
   plot(dep.time,dep.depthm);
   hold on;
   y1 = ylim; t1 = tme(1);
   plot([t1 t1], y1,'r'); grid on;
   datetick; xlabel('Time'); ylabel('Depth');
   fprintf(1,'Uncorrected water depth at anchor drop %5.1f m \n',wd(1));
   ich = input('Enter 1 to accept or 2 to change values: ');
   if ich == 2
       wd=input('Enter uncorrected water depth in metres: ');
   end
   %and correct
   corr_struct = mcarter(lat(1),lon(1),wd);
   wd_corr = corr_struct.cordep;

end

% open output file
fileout = fullfile(pd.trilatdir, [loc_name '_triangle.txt']);
disp(['Opening output file ' fileout]);
disp('If file exists contents will be overwritten'); pause(0.1)
iout = fopen(fileout,'w+');

% If no fixes then just give anchor drop position info
no_fixes=length(range);
if no_fixes == 0
    fprintf(1,'Anchor drop at: %s %8.4f %8.4f Corr. water depth: %6.1f \n', ...
        datestr(tmeD,31),latD,lonD,wd_corr);
    fprintf(1,'%8.2f N %8.2f W \n',dd2dm(latD),dd2dm(lonD))

    fprintf(iout,'%s  %s \n',loc_name,datestr(tmeD,1));
    fprintf(iout,'Anchor drop at: %s %8.4f %8.4f Corr. water depth: %6.1f \n', ...
        datestr(tmeD,31),latD,lonD,wd_corr);
    fprintf(iout,'%8.4f %8.4f \n',dd2dm(latD),dd2dm(lonD));
    isn = lat>0; ise = lon>0;
    latdeg = floor(abs(latD));
    latmin = 60*(abs(latD)-latdeg);
    londeg = floor(abs(lonD));
    lonmin = 60*(abs(lonD)-londeg);
    if isn; latdstr = 'N'; else; latdstr = 'S'; end
    if ise; londstr = 'E'; else; londstr = 'W'; end
    title4 = sprintf('Latitude %i %5.2f %s, Longitude %i %5.2f %s',latdeg,latmin,latdstr,londeg,lonmin,londstr);
    fprintf(iout,'%s',title4);
    return
end

% Work out effective horizontal range 
% Sound speed correction: apply the carter table correction (based on
% vertical profiles of sound speed) as though slant range is depth only if
% the vector is close to vertical; otherwise just use the same % correction
% as for the anchor-drop depth 
mo = (wd_corr-wd)./wd;
rangeC = range*(1+mo);
m = acos(wd./range)<=10/180*pi;
corr_struct = mcarter(lat(m),lon(m),range(m));
rangeC(m) = corr_struct.cordep;
rangeh = sqrt(rangeC.^2 - (wd_corr-rht-td).^2);

% Open a new plot
figure(2); clf
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
% plot ship's track
plot(pos.lon,pos.lat,'k--')

%make circles
xdegr = -180:2:180;
[clon,clat,~] = m_fdist(lon,lat,xdegr,rangeh); %works as vector as long as length(rangeh) is not the same as length(xdegr), which it shouldn't be!
clon = clon-360;
hl = plot([lon lon]',[lat lat]','b+'); if length(hl)>3; set(hl(4:end),'color','r'); end
hlc = plot(clon',clat','b'); if length(hlc)>3; set(hlc(4:end),'color','r'); end
titletext={sprintf('Triangulation Survey for: %s',loc_name);...
    sprintf('Corrected water depth: %5.0f m.',wd_corr);... 
    sprintf('Release Height: %3.1f m. Transducer depth: %3.1f',rht,td)};
title(titletext);
xlabel('Longitude'); ylabel('Latitude')


% Calculate triangulation point by solving least squares problem
% Nb take differences of equations to make linear
[long1,latg1,~,~,reser0] = solve_anchor(lonD,latD,lon(:),lat(:),rangeh(:)',no_fixes);
[APlon,APlat,x,y,reser] = solve_anchor(long1,latg1,lon(:),lat(:),rangeh(:)',no_fixes);
if max(reser0-reser) > 0.1
    fprintf(1,'Something not right')
    keyboard
end

% plot estimated anchor seabed locations and anchor drop point
h1 = plot(APlon,APlat,'mo');
text((east-west)*0.1+west,(north-south)*0.90+south,'Anchor drop','color','g');
plot(lonD,latD,'g+');
plot(lonD,latD,'go');
grid on
%hh = plot(tlon,tlat,'MarkerSize',20,'LineWidth',1.5);

% prompt user to determine anchor position from plot %***allow user to pick
% point using ginput? either way this plot should be in a geographical
% projection or at least be scaled so that angles are about right***
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
latdeg = floor(abs(APlat));
londeg = floor(abs(APlon));
latmin = 60*(abs(APlat)-latdeg);
lonmin = 60*(abs(APlon)-londeg);
if APlat>0
    latstr = 'N';
else 
    latstr = 'S';
end
if APlon>0
    lonstr = 'E';
else
    lonstr = 'W';
end

% Calculate distance from drop to moored positions
fallback=sw_dist([latD APlat],[lonD APlon],'km')*1000;
fallback=sprintf('%5.0f',fallback);
fallbacknm = sw_dist([latD APlat],[lonD APlon],'nm');
fallbacknm = sprintf('%5.2f',fallbacknm);

if numel(~isnan(targ_loc))==2
    dist2targ = sw_dist([targ_loc(1) APlat],[targ_loc(2) APlon],'km')*1000;
else
    dist2targ=NaN;
end

plot(APlon,APlat,'r+');
plot(APlon,APlat,'ro');
plot(lonD,latD,'g+');
plot(lonD,latD,'go');
title4 = sprintf('Latitude %i %5.2f %s, Longitude %i %5.2f %s',latdeg,latmin,latstr,londeg,lonmin,lonstr);
title5 = sprintf('Distance from target: %3.0fm',dist2targ);
titletext = [titletext;...
    title4;...
    title5];
title(titletext);

% Text for plot
text((east-west)*0.1+west,(north-south)*0.95+south,['Fallback = ' fallback 'm/' fallbacknm 'nm'],'color','k');
text((east-west)*0.1+west,(north-south)*0.90+south,'Anchor drop','color','g');
text((east-west)*0.1+west,(north-south)*0.85+south,'Anchor seabed location','color','r');

% Finish plot
xlim([west east]);
ylim([south north]);

% Plot to file is saving
if isempty(loc_name)
    return
else
    plotout = fullfile(pd.trilatposdir, [loc_name '_triangle']);
    print('-depsc', plotout)
    system(['cat ' fileout]);
end

% again check against mapped and ship track bathymetry
[swtx,mdep] = get_mapbathy(loc_name,APlon,APlat,wd_corr);
% bathymetry when ship passed over the calculated anchor point
[pos, dep] = get_uway([tmeD-0.5 max(tme)+0.5]);
[wd_sh, dis_sh, im] = shipdep(pos, dep, APlon, APlat);
fprintf(1,'Depth at closest ship track to anchor position was %6.1f m (%6.1f nm away)\n',wd_sh,dis_sh(im));
if abs(wd_sh-wd)>200 || dis_sh(im)>2
    figure(5); clf; 
    if ~isempty(mdep); contour(mdep.lon,mdep.lat,mdep.depth); end
    hold on
    scatter(pos.lon,pos.lat,10,dep.depthm,'filled'); colorbar
    plot(APlon,APlat,'ok',lonD,latD,'xk'); grid on
    axis([APlon-.1 APlon+.1 APlat-.1 APlat+.1])
    title('raw underway depths (colors), anchor position (o), drop position (x)')
    fprintf(1,'there may be bad depth data;\n');
    a = input('edit data and recalculate (e), select new point (s), or enter to continue  ', 's');
    if strcmp(a,'e')
        clear igood
        fprintf(1,'set igood, index of good points in dep.depthm, then dbcont\n');
        keyboard
        if exist('igood','var')
            dep.time = dep.time(igood); dep.depthm = dep.depthm(igood);
            pos.time = pos.time(igood); pos.lat = pos.lat(igood); pos.lon = pos.lon(igood);
            [wd_sh, dis_sh, im] = shipdep(pos, dep, APlon, APlat);
            fprintf(1,'Depth at closest ship track (with good depth) to anchor position was %6.1f m (%6.1f nm away)\n',wd_sh,dis_sh(im));
        end
    elseif strcmp(a,'s')
        fprintf(1,'zoom if necessary then press enter\n'); pause
        fprintf(1,'choose new point on figure\n');
        [x,y] = ginput(1);
        if ~isempty(x)
            [wd_sh,~,im] = shipdep(pos, dep, x, y);
            dnew = sw_dist([pos.lat(im) APlat],[pos.lon(im) APlon]);
            fprintf(1,'Depth at chosen ship track point was %6.1f m (%6.1f nm away)\n',wd_sh,dnew);
        end
    end
end
corr_struct = mcarter(APlat,APlon,wd_sh);
wd_corr_sh = corr_struct.cordep;


% Now data that we save in the output file
fprintf(iout,'%s  %s \n',loc_name,datestr(tmeD,1));
fprintf(iout,'Anchor drop at: %s %8.4f %8.4f Corr. water depth: %6.1f \n', ...
    datestr(tmeD,31),latD,lonD,wd_corr);
fprintf(iout,'%s',swtx);
fprintf(iout,'Date      Time       Lat   Lon  Slant range  Horiz range Residual \n');
for i = 1:no_fixes
    fprintf(iout,'%s  %8.4f %8.4f %7.0f %7.0f %7.0f \n', ...
        datestr(tme(i),31),lat(i),lon(i),range(i),rangeh(i),reser(i));
end
fprintf(iout,'Release Height: %3.1f m. Transducer depth: %3.1f \n',rht,td);
fprintf(iout,'Trilaterated position \n');
fprintf(iout,'%8.4f %8.4f \n',APlat,APlon);
fprintf(iout,'%s \n',title4);
fprintf(iout,'Fallback  %s m \n',fallback);
fprintf(iout,'%s \n', title5);
fprintf(iout,'Comments: %s \n',us_com);


function [pos, dep] = get_uway(t_range)
% load nav and depth data, interpolate depth to nav times, 
% subsample to 30 s

global MEXEC_G

%where to get position and uncorrected (singlebeam) depth on this ship
opt1 = 'ship'; opt2 = 'datasys_best'; get_cropt %default_navstream

if exist('default_navstream','var')
    nav_stream = default_navstream;
else
    nav_stream = {'techsas' 'cnav'; 'scs' 'posfur'; 'rvdas' 'posmv_gpgga'};
    nav_stream = nav_stream(strcmp(MEXEC_G.Mshipdatasystem,nav_stream(:,1)),2);
end

dv1 = datevec(t_range(1));
dv2 = datevec(t_range(2));

% Now get data from underway data system
switch MEXEC_G.Mshipdatasystem
    case 'techsas'
        dep_stream = 'sim';
        pos = mtload(nav_stream,dv1,dv2,'time lat long');
        pos.time = pos.time + MEXEC_G.uway_torg;
        pos.lon = pos.long;
        dep = mtload(dep_stream,dv1,dv2,'time depthm');
        dep.time = dep.time + MEXEC_G.uway_torg;
        dep.depthm = dep.waterdepth;
    case 'rvdas'
        nav_vars = 'time latitude longitude';
        dep_stream = 'ea640_sddbs';
        dep_vars = 'time waterdepth';
        pos = mrload(nav_stream,dv1,dv2,'time latitude longitude','q');
        pos.lat = pos.latitude;
        pos.lon = pos.longitude;
        pos.time = pos.dnum;
        dep = mrload(dep_stream,dv1,dv2,'time waterdepthmetres','q');
        dep.depthm = dep.waterdepthfromsurface;
        dep.time = dep.dnum;
    case 'scs'
        nav_vars   = 'time GPS-Furuno-GGA-lat GPS-Furuno-GGA-lon';
        dep_stream = 'singleb';
        dep_vars = 'time waterdepth';
        pos = msload(nav_stream,dv1,dv2,nav_vars);
        pos.lat = pos.GPS_Furuno_GGA_lat;
        pos.lon = pos.GPS_Furuno_GGA_lon;
        [pos.time, IA, ~] = unique(pos.time);
        pos.lat = pos.lat(IA);
        pos.lon = pos.lon(IA);
        dep = msload(dep_stream,dv1,dv2,dep_vars);
        dep.depthm = dep.waterdepth;
        [dep.time, IA, ~] = unique(dep.time);
        dep.depthm = dep.depthm(IA);
end
iit = 1:30:length(pos.time);
pos.time = pos.time(iit); pos.lat = pos.lat(iit); pos.lon = pos.lon(iit);
igood = dep.depthm>0;
dep.depthm = interp1(dep.time(igood),dep.depthm(igood),pos.time);
dep.time = pos.time;

function [swtx, mdep] = get_mapbathy(loc_name, APlon, APlat, wd_corr)
% Check bathymetry from swath if available - might want to rerun again after this
if strncmp(loc_name,'wb',2) & ~strcmp(loc_name(3),'4')
    topo_map = 'gr_kn182_plot.mat';
elseif strncmp(loc_name,'ma',2)
    if strncmp(loc_name,'mar3',4)
        topo_map = 'mar34.mat';
    elseif sum(strncmp(loc_name,{'mar2' 'mar1'},4))>0
        topo_map = 'mar12.mat';
    elseif strncmp(loc_name,'mar0',4)
        topo_map = 'mar0_JC064_swath.mat';
    end
else
    topo_map = 'none';
end
if strcmp(topo_map,'none')
    swtx = sprintf('No swath data available depth used was %7.1f \n',wd_corr);
    mdep = [];
else
    rdpath = '/data/pstar/projects/rpdmoc/other_datasets/for_sea/bathym_data/from_cruises/';
    mdep = load(fullfile(rdpath,topo_map),'lon','lat','depth');
    newdep = interp2(mdep.lon,mdep.lat,mdep.depth,APlon,APlat);
    swtx = sprintf('Depth used was %7.1f,\nbut corrected depth from swath map at trilaterated position is %7.1f;\n you may want to edit input file and rerun\n\n',wd_corr,newdep);
end
fprintf(1,'%s',swtx)

function [wd_sh,dis_sh,im] = shipdep(pos, dep, APlon, APlat)

dy = [pos.lat repmat(APlat,size(pos.lat))]'; dx = [pos.lon repmat(APlon,size(pos.lon))]';
dis_sh = sw_dist(dy(:),dx(:)); dis_sh = dis_sh(1:2:end);

m = dis_sh<5; 
dis_sh = dis_sh(m); 
pos.lon = pos.lon(m); pos.lat = pos.lat(m); pos.time = pos.time(m);
dep.time = dep.time(m); dep.depthm = dep.depthm(m);

[~,im] = min(dis_sh);

wd_sh = dinterp1(dep.time,dep.depthm,pos.time(im));
