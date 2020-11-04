function rapd_daylight(varargin)
% A script to get hours of daylight at RAPID mooring sites
% calls the suncycle program
% 
% Input: position, time  (both are optional)
% Default time is today otherwise enter as vector [yyyy mm dd]
% Default position is here (using mtposinfo) or enter as either:
%  - a mooring postition name
%  - lat and long in format   ddmm.nn
% e.g.
%   rapid_daylight()
%   rapid_daylight('ebh1',[2014 05 28])
%   rapid_daylight([2713.34 -1525.36],[2014 05 28])

m_setup;
addpath /local/users/pstar/rpdmoc/rapid/data/exec/jc103/mfiles/misc

fprintf(1,'\n\n ****** Check location information is correct and up to date ******** \n')
fprintf(1,' ******  -- Reading data from moor_pos.dat -- ******* \n\n')
fomp = fopen('moor_pos.dat');
aabb = textscan(fomp,'%s%f%f%f');
mooring.names = aabb{1};
mooring.lat = dm2dd(aabb{2});
mooring.lon = dm2dd(aabb{3});
mooring.depth = aabb{4};

if nargin == 0
	[latp lonp] = mtposinfo;
	day = datenum(date);
elseif nargin > 0
	if ~ischar(varargin{1})
		posn = varargin{1};
	    latp = dm2dd(posn(1));
	    lonp = dm2dd(posn(2));
	else
		moor_list = char(mooring.names);
	    this_moor = varargin{1};
  	  	imoor = find(strcmp(this_moor,mooring.names));
		latp = mooring.lat(imoor);
		lonp = mooring.lon(imoor);
	end
	if nargin == 1
	    day = datenum(date)
	elseif nargin == 2
	    day = datenum(varargin{2});
	end
else
	fprintf(1,'Error in input variables');
	return
end

dayvec = datevec(day);
dayvec = dayvec(1:3);

[SunRiseSet] = suncycle(latp,lonp,dayvec);

sunrise = SunRiseSet(1);
sunset  = SunRiseSet(2);

[SRhr,SRmn] = dd2dm2(sunrise);
[SShr,SSmn] = dd2dm2(sunset);
[latd latm] = dd2dm2(latp);
[lond lonm] = dd2dm2(lonp);

if exist('this_moor')
	fprintf('Mooring: %s \n',this_moor)
end
fprintf(1,'For location %i°%4.1f N %i°%4.1f W on %i-%i-%i \n',latd,latm,lond,lonm,dayvec);
fprintf(1,'Sunrise at %2.2i:%2.2i UTC and sunset at %2.2i:%2.2i UTC \n\n',SRhr,floor(SRmn),SShr,floor(SSmn));

