function [lons,lats] = plan_triangle(lon0,lat0,mdist,angls)
% [lons,lats] = plan_triabgle(lon0,lat0,mdist,angls)
% Inputs:
%  - lon0,lat0 = long and lat of expected anchor position in decimal degrees
%  - mdist     = distance in metres from mooring of triangulation waypoints
%  - angls     = vector of bearings, one for each waypoint
% Ouputs
%  - lons,lats = decimal positions of waypoints
%
  angls = angls-90;
  [lons,lats,a21] = m_fdist(lon0,lat0,angls,mdist);
for i = 1:length(lons)
 if lons(i) > 180
   lons(i) = lons(i)-360;
 end
 end
%
  fprintf(1,'%8.4f   %8.4f \n \n',lon0,lat0);
  for i = 1:length(angls)
    fprintf(1,'%8.4f   %8.4f \n',lons(i),lats(i));
    end
%
  fprintf(1,' -- \n')

  fprintf(1,'%8.2f   %8.2f \n \n',dd2dm([lon0 lat0]));

  for i = 1:length(angls)
	apos = dd2dm([lons' lats']);
    fprintf(1,'%8.2f   %8.2f \n',apos(i,:));
    end
   

