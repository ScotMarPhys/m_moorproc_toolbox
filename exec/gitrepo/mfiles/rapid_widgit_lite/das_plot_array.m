function das_plot_array
% Some examples on use of m_map
mooring.names = {'ebh4','ebh3', 'ebh2', 'ebh1', 'ebhi', 'eb1',  'nog',  'mar3', 'mar2', 'mar1', ...
                  'mar0', 'wb6',  'wb4',  'wbh2', 'wb2',  'wb1',  'wbadcp'};
mooring.lat =  [27.8501, 27.805, 27.6147,27.22,24.933, 23.7577,23.7710,23.807, 24.1775,24.2655, ...
               25.1047,26.4950,26.4782,26.4823,26.5136,26.5012,26.52492];
mooring.lon = -[13.5408, 13.7468,14.2109,15.43,21.2732,24.1582,41.0988,41.0984,49.763, 49.75, ...
               52.0102,70.5223,75.7028,76.6276,76.7379,76.8149,76.86803];
mooring.depth = [1062, 1421,2022, 2997, 4505, 5096, 4247, 5061, 5216, 5215, 5529, 5491, 4687, 4729. 3849, 1387, 608];


load ~/rpdmoc/coastline_data/eez_mat.mat

% Example 2 - a plot of bathymetry along the RAPID line
% Data set currently set to SMith and Sandwell v10.1 (1 minute resolution)
figure
region = [15 35 -80 -10];
region = [26 28.75 -16 -12];

south = region(1);
north = region(2);
west = region(3);
east = region(4);
m_proj('lambert','lon',[west, east],'lat',[south, north]);
% load data covering region
dd = 0.01;
[elevations,elat,elon]=mygrid_sand([south-dd,north+dd,west-dd,east+dd],1);
elon = elon-360;
% Color display of data
%m_pcolor(elon,elat,elevations);
%shading flat;
% use intermediate resolution coastline
m_gshhs_l('patch',[0.43 0.55 0.24]);
m_grid('box','on','color','k','linewidth',[1],'fontsize',[12]);
hold on
keyboard
% Plot soem contours
conts = -[100 200 500 1000 2000 5000];
hmc = m_contour(elon,elat,elevations,conts);
% plot some mooring positions


h = m_plot(mooring.lon,mooring.lat,'k+');
dx = -.1;dy = 0.1;
for i = 1:length(mooring.names)
h2 = m_text(mooring.lon(i)+dx,mooring.lat(i)+dy,char(mooring.names(i)));
end

ieez = eez_lon > west & eez_lon < east &eez_lat > south & eez_lat < north;
m_plot(eez_lon,eez_lat,'r+')
keyboard
