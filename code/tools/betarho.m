function [bet]=betarho(p,t,s)
% function [bet]=betarho(p,t,s)
%
% calculate haline expansion coefficient beta
%
% input  : 	p	- pressure [dbar]
%		t	- in situ temperature [degrees]
%		s	- salinity [psu]
%
% output	bet	- haline expansion coeff [d-rho/d-S]
%
% uses : alpha 
% 
% version 1.0.0		last change 1.8.1995

% G.Krahmann, IfM Kiel, Aug 1995

bet=(1./alpha(p,t,s+0.00001)-1./alpha(p,t,s))/0.00001.*alpha(p,t,s);
