function [up,down,dm] = thickness(p,lat,sig,int)

% THICKNESS  [up,down,dm] = thickness(p,lat,sig,int)
%
%            calculate the thickness of a density layer
%
%            input  : p        pressure matrix                  [dbar]
%                     lat      latitude                         [deg] 
%                     sig      density (potential or in situ)   [kg/m^3]
%                     int      density layer, e.g. [27.74 27.80]
%
%            output : up       upper density boundary values [kg/m^3]
%                     down     lower density boundary values [kg/m^3]
%                     dm       layer thickness               [m]
%
%            uses   : p2z.m (gravity.m), nearest.m
%
%            version 1.0, d.kieke IfMK, 01.04.1998

% --- CALCULATE DEPTH FROM PRESSURE MATRIX -------------------------------------

[m,n] = size(p);
lat   = lat(:)';
z     = -p2z(p,lat(ones(m,1),:));

% --- FIND DENSITY VALUES AND THEIR INDICES IN DENSITY MATRIX ------------------

up   = ones(n,1);
nxu  = ones(n,1);
down = ones(n,1);
nxd  = ones(n,1);
d1   = ones(n,1);
d2   = ones(n,1);

for i = 1:n
  [up(i),nxu(i)]   = nearest(int(1),sig(:,i));
  [down(i),nxd(i)] = nearest(int(2),sig(:,i));
  d1(i)            = z(nxu(i),i);
  d2(i)            = z(nxd(i),i);
  i = i+1;
end

dm = d1-d2;

