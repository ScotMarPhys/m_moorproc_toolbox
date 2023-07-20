function rapid_plot_moor(varargin)
% This script plots a mooring locations
% *********** PLEASE CHECK POSITIONS ARE CORRECT *************
% The positions are stored in the file:
%  - moor_pos.dat

% There are variable number of inputs:  this_moor,plot_size,cont_int
% Usually only need the first one
%   - this_moor - is a string and can have values:
%   	- 'ship'  Will create plot cetred on ship's position
% 		- 'check' Will list all mooring names, positions and nominal depths
%   - plot_size - is a real number and is size of plot in nautical miles
%   - cont_int - is contour interval usually do not need to set this
%
% Calls following scripts:  rapid_map.m, dirav.m, dm2dd.m, dd2dm2.m
%
% Update DAS April 2014 on JC103

% Default values of input variables
this_moor = 'check';
plot_size = 8;		% Default choice
cont_int = -999;   	% < 0 means automatic choice

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


% Some preliminaries
%m_setup;      % m_setup needed as we will access data from Techsas

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
  if strmatch(this_moor, 'user') == 1
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
dlon = dlat*cos(pi*mlat/180);
north = mlat + dlat/2;
south = mlat - dlat/2;
west = mlon - dlon/2;
east = mlon + dlon/2;

rapid_map(this_moor,mlat,mlon,mdepth,north,south,east,west,cont_int);

