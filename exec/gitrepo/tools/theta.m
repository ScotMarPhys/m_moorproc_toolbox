function ptemp = theta(p,t,s,p0)
%THETA Computes local potential temperature at reference pressure.
%       [PTEMP] = THETA(P,T,S,P0) is the local potential temperature
%       at reference pressure P0 using Bryden 1973 polynomial for
%       adiabatic lapse rate and Runge-Kutta fourth order integration
%       algorithm.
%
%       Units:
%               Pressure        P, P0   dbar
%               Temperature     T       deg C
%               Salinity        S       NSU
%	Defaults:
%		P0		0 dbar
%
%       Checkvalue:
%               THETA(10000,40,40,0) = 36.89072
%

%       18/02/93, C. Mertens, IfM Kiel, changed to Matlab
% added default p0=0dbar	G.Krahmann, IfM Kiel, Mar 1996

[m,n] = size(p) ;
if n == 1 ,
        [m,n] = size(t) ;
        p = p*ones(1,n) ;
end

if nargin<4
  p0=0;
end

p = p/10 ; 
p0 = p0/10 ;
h = p0 - p ;
x = h.*atg(p,t,s) ;
t = t + 0.5*x ;
q = x ;
p = p + 0.5*h ;
x = h.*atg(p,t,s) ;
t = t + 0.29289322*(x - q) ;
q = 0.58578644*x + 0.121320344*q ;
x = h.*atg(p,t,s) ;
t = t + 1.707106781*(x - q) ;
q = 3.414213562*x - 4.121320344*q ;
p = p + 0.5*h ;
x = h.*atg(p,t,s) ;
ptemp = t + (x - 2*q)/6 ;

