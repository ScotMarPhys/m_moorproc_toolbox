function [Gamma] = alr(p,t,s)
% function [Gamma] = alr(p,t,s)
%
% adiabatic lapse rate of sea water
% after Landolt-Boernstein pp 250-252
%
% input  : 	p		: pressure [dbar]
%		t		: in situ temperature [degrees Celsius]
%		s		: salinity [psu]
%
% output :	Gamma		: adiabatic lapse rate [K/dbar]
%
% check value :		alr(10000,30,40) = 2.8774579e-4 K/dbar
%
% version 1.1.0		last change 06.09.1995

p2=p.^2;
t2=t.^2;
t3=t.^3;
s35=s-35;

a=[3.5803e-5,8.52587e-6,-6.83605e-8,6.6228e-10];
b=[1.8932e-6,-4.23935e-8];
c=[1.8741e-8,-6.7795e-10,8.7330e-12,-5.4481e-14];
d=[-1.1351e-10,2.7759e-12];
e=[-4.6206e-13,1.8676e-14,-2.1687e-16];

Gamma=		 a(1) + a(2)*t + a(3)*t2 + a(4)*t3 +...
	s35.*	(b(1) + b(2)*t) +...
	p.*	(c(1) + c(2)*t + c(3)*t2 + c(4)*t3 +...
		 s35 .* (d(1) + d(2)*t)) +...
	p2.*	(e(1) + e(2)*t + e(3)*t2);
