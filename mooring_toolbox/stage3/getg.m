function g=getg(lat,z)
% function g=getg(lat,z)
% Acceleration due to gravity in m/s^2
% function of latitude in degree
% and height in m
% GILL pp 597	(the formula is wrong !!!)
%
% version 1.1.0		last change 17.11.1995

% author unknown
% corrected formula	G.Krahmann, IfM Kiel, Nov 1995

if nargin < 2, z=0; end
if nargin >0
 a=6371e3;
 if (size(lat,1) < size(z) | size(lat,2) < size(z))
  lat = lat(:); 
  lat = ones(length(z),1)*lat';
 end
 lr=lat*pi/180;
% g=(9.78032+0.005172*(sin(lr).^2)-0.00006*(sin(2*lr).^2))./(1+z/a).^2 ;
 g=9.78032*(1+0.005172*(sin(lr).^2)-0.00006*(sin(2*lr).^2))./(1+z/a).^2 ;
else
 g=9.7976;
end
