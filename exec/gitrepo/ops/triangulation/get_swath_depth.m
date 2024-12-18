function get_swath_depth(latddmm,londdmm)
% Find depth from swath data
% check_swath(latddmm,londdmm)  
m_setup
this_cruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
rootdir = ['/noc/users/pstar/rpdmoc/rapid/data/exec/' this_cruise];


% To access other mfiles
addpath([rootdir '/mfiles/misc']);

%latddmm = input('Enter latitude in format DDMM.MM: \n');
%londdmm = input('Enter longitude in format DDMM.MM: \n');

latD = dm2dd(latddmm);
lonD = dm2dd(londdmm);

topo_maps = {'gr_kn182_plot.mat','mar34.mat','mar12.mat','mar0_JC064_swath.mat'};
fprintf(1,'Select which map to use: \n');
for i = 1:4 ...
	fprintf(1,' - %s \n ',topo_maps{i})
end
isel = input('enter choise: ');

% Check bathymetry from swath if available - might want to rerun again after this']
if isel == 1;
    topo_map = 'gr_kn182_plot.mat';
elseif isel == 2;
  	topo_map = 'mar34.mat';
elseif isel == 3;
    topo_map = 'mar12.mat';
elseif isel == 4;
	topo_map = 'mar0_JC064_swath.mat';
end

rdpath = [rootdir '/mfiles/rapid_widgit_v2/data/'];
load([rdpath topo_map]);
newdep = interp2(lon,lat,depth,lonD,latD);
 ...
[latdeg,latmin] = dd2dm2(latD);
[londeg,lonmin] = dd2dm2(lonD);
fprintf(1,'Depth at positoin is %7.1f \n',newdep); 
pos_text = sprintf('Latitude %i %5.2f N, Longitude %i %5.2f W \n',latdeg,latmin,londeg,lonmin);
fprintf(1,'%s',pos_text)
