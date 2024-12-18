function Triangle5
% Plot positions and range circles 
% This version user enters only time information in file  <name>_time.txt
% and position is found from navigation files on Techsas.
% Format of text file shoudl be 
%    YYYY MM DD HH MM SS range  
% 
% Saves plotfile only
% This version uses calls to m_map toolbox to use m_fdist and m_idist to calculate
% distances on spheroid.
% DAS November 2012
ccl = 'rmgybcrmgybc';
global MEXEC_G
global MOORPROC_G
opt1 = 'ship'; opt2 = 'datasys_best'; get_cropt %default_navstream
dep_stream = 'ea600m';
%dirout = '/noc/users/pstar/rpdmoc/rapid/data/moor/raw/d382/moor_positions/';
this_cruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
dirout = fullfile(MOORPROC_G.moordatadir,'raw',this_cruise,'moor_positions');

fprintf(1,'Enter mooring name (e.g. ebh3). Output will be saved to \n%s\n',dirout);
loc_name = input('and times will be read from <name>_time.txt: ','s');

if strmatch(loc_name, '')
    return
else
    if exist(dirout,'dir')
        filein = fullfile(dirout, [loc_name '_times.txt']);
        fileout = fullfile(dirout, [loc_name '_triangle.txt']);
        disp(['Reading times from ' filein]);
        if exist(filein,'file')
            indata = load(filein);
        else
            warning('File of times %s does not exist',filein)
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
wd = input('Enter uncorrected water depth: ');

% Read in data from text file assume first point is anchor drop
for i = 1:size(indata,1)
    dvec = indata(i,1:6);
    tme(i) = datenum(dvec);
    range(i) = indata(i,7);
    if size(indata,2) == 8
        dptvl(i) = indata(i,8);
    end
    if range(i) == 0
        iplt(i) = 0;
    else
        iplt(i) = 1;
    end
end

% No of points in input file
no_fixes=length(range);

% Now get data from Techsas
switch MEXEC_G.Mshipdatasystem
    case 'techsas'
        pos = mtload(default_navstream,datevec(tme(1)-0.04), ...
            datevec(tme(end)+0.01),'time lat long','q');
        pos.lon = pos.long;
        pos.time = pos.time + MEXEC_G.Mtechsas_torg;
    case 'rvdas'
        pos = mrload(default_navstream,datevec(tme(1)-0.04,datevec(tme(end)+0.01)),'time latitude longitude');
        pos.lon = pos.longitude; pos.lat = pos.latitude;
        pos.time = pos.dnum;
end

% interpolate for positions then depths
latp = interp1(pos.time,pos.lat,tme);
lonp = interp1(pos.time,pos.lon,tme);

latM = mean(latp);
lonM = mean(lonp);
%wd_corr = carter(latM,lonM,wd);
corr_struct = 	mcarter(latM,lonM,wd);
wd_corr = corr_struct.cordep;

% Work out effective horizontal range
for i=1:no_fixes
    rangeh(1,i)=sqrt((range(i))^2 - (wd-rht-td)^2);
    rangeh(2,i)=sqrt((range(i))^2 - (wd-25-rht-td)^2);
end
disp(rangeh)
% Open a new plot
figure
% Decide on boudaries for plot
mpd = 111.2*1000;
axylim = 1.05;
mxrangeh = max(rangeh(1,:))';
south = latM -axylim*mxrangeh/mpd;
north = latM +axylim*mxrangeh/mpd;
west =  lonM -axylim*mxrangeh/(mpd*cos(latM*pi/180));
east =  lonM +axylim*mxrangeh/(mpd*cos(latM*pi/180));
m_proj('lambert','lon',[west, east],'lat',[south, north]);
% Following is if use m_map to plot - but this does not allow zooming so have disabled
% m_grid('box','on','color','k','linewidth',[1],'fontsize',[14]);
hold on
grid on
xdegr = -180:2:180;
for i=1:no_fixes
    plot(lonp(i),latp(i),[ccl(i) '+']);
    if range(i) > 0
        [clon,clat,az3] = m_fdist(lonp(i),latp(i),xdegr,rangeh(1,i));
        clon = clon-360;
        plot(clon,clat,ccl(i));
        [clon,clat,az3] = m_fdist(lonp(i),latp(i),xdegr,rangeh(2,i));
        clon = clon-360;
        plot(clon,clat,[ccl(i) '--']);
    end
end

titletext1=['Triangulation Survey for: ',loc_name];
titletext2=sprintf( ...
    'Corrected water depth: %5.0f m. Release Height: %3.0f m. Transducer depth: %3.0f', ...
    wd_corr,rht,td);
title({titletext1;titletext2});
xlabel('Longitude'); ylabel('Latitude')


% determine anchor position from figure(1)
disp('Use figure to determine anchor seabed position:')
APlat = input('Latitude = ');
APlon = input('Longitude = ');

% Degrees and minutes
latdeg = floor(APlat);
londeg = floor(-APlon);
latmin = 60*(APlat-latdeg);
lonmin = 60*(-APlon-londeg);

title4 = sprintf('Latitude %i %5.2f N, Longitude %i %5.2f W',latdeg,latmin,londeg,lonmin);


plot(APlon,APlat,'k+','MarkerSize',20);

title({titletext1;titletext2;title4});

% Finish plot
xlim([west east]);
ylim([south north]);

%iout = 1;
fprintf(iout,'Date      Time       Lat   Lon  Slant range  Horiz range Residual \n');
for i = 1:no_fixes
    fprintf(iout,'%s  %8.4f %8.4f %7.0f %7.0f \n', ...
        datestr(tme(i),31),latp(i),lonp(i),range(1,i),rangeh(i));
end


figure
m_proj('lambert','lon',[west, east],'lat',[south, north]);
m_grid('box','on','color','k','linewidth',1,'fontsize',14);
hold on

if strmatch(loc_name,'rec-mar0')
    load('/noc/users/pstar/rpdmoc/rapid/data/exec/d382/mfiles/rapid_widgit_v2/data/mar0_JC064_swath.mat')
elseif strmatch(loc_name,'rec-mar1l6')
    load('/noc/users/pstar/rpdmoc/rapid/data/exec/d382/mfiles/rapid_widgit_v2/data/mar12.mat')
elseif strmatch(loc_name,'wb1-9')
    load('/noc/users/pstar/rpdmoc/rapid/data/exec/d382/mfiles/rapid_widgit_v2/data/gr_kn182_plot.mat')
else
    dd = 0.1;
    fprintf(1,'using Smith and Sandwell\n')
    [depth,lat,lon]=mygrid_sand([south-dd,north+dd,west-dd,east+dd],1);
    lon = lon - 360;
    depth = -depth;
end

contourd = [500:100:4000];
m_contour(lon,lat,depth,contourd,'k')
for i=1:no_fixes
    m_plot(lonp(i),latp(i),[ccl(i) '+']);
    if range(i) > 0
        [clon,clat,az3] = m_fdist(lonp(i),latp(i),xdegr,rangeh(1,i));
        clon = clon-360;
        m_plot(clon,clat,ccl(i));
        [clon,clat,az3] = m_fdist(lonp(i),latp(i),xdegr,rangeh(2,i));
        clon = clon-360;
        m_plot(clon,clat,[ccl(i) '--']);
    end
end
m_plot(APlon,APlat,'k+','MarkerSize',20);
title({titletext1;titletext2;title4});
%keyboard
