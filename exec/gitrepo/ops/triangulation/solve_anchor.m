function [lonA,latA,x,y,resid] = solve_anchor(lon0,lat0,lon,lat,rangeh,no_fixes);
% Solve for positon of anchor from positiona dn hroizontal ranges
% solve equation of form a * sol = b where sol is the correction to hte anchor position

% express as east and north distances in m from first estimate
for j=1:no_fixes
  [dis,az12,az21] = m_idist(lon0,lat0,lon(j),lat(j));
  az = az12*pi/180;
  x(j) = dis*sin(az);
  y(j) = dis*cos(az);
end

% Now calculate a and b
xyr = x.^2+y.^2-rangeh.^2;;
for j = 1:no_fixes-1
  b(j) = xyr(j)-xyr(j+1);
  a(j,1) = 2*(x(j)-x(j+1));
  a(j,2) = 2*(y(j)-y(j+1));
end

j = no_fixes;
b(j) = xyr(j)-xyr(1);
a(j,1) = 2*(x(j)-x(1));
a(j,2) = 2*(y(j)-y(1));

% This should work for 3 or more ranges
sol = b/a';

% Residual errors:
resid  = abs( sqrt((x-sol(1)).^2+(y-sol(2)).^2)-rangeh);

% The answer
disx = sqrt(sol(1).^2+sol(2).^2);
az = atan2(sol(1),sol(2))*180/pi;
 
% now express as longitud and latitude
[lonA,latA,az3] = m_fdist(lon0,lat0,az,disx);
lonA = lonA-360;
