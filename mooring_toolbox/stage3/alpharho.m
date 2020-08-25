function [alp]=alpharho(p,t,s)
% function [alp]=alpharho(p,t,s)
%
% calculate thermal expansion coefficient alpha
%
% input  : 	p	- pressure [dbar]
%		t	- in situ temperature [degrees]
%		s	- salinity [psu]
%
% output	alp	- thermal expansion coeff d-rho/d-T  [K^-1]
%
% check values :	alpharho(10000,40,40) = 4.1810597e-4 K^-1
%
% uses : alpha 
% 
% version 1.0.1		last change 06.09.1995

% G.Krahmann, IfM Kiel, Aug 1995
% added check value	G.Krahmann, Sep 1995

alp=-(1./alpha(p,t+0.00001,s)-1./alpha(p,t,s))/0.00001.*alpha(p,t,s);
