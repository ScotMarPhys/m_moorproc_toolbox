function [rr]=rrho(p,t,s,dt,ds)
% function [rr]=rrho(p,t,s,dt,ds)
%
% calculate density ratio R_rho
%
% input  : 	p	- pressure [dbar]
%		t	- temperature [degrees]
%		s	- salinity [psu]
%		dt	- temperature difference [degrees]
%			  (should be potential temperature difference)
%		ds	- salinity difference [psu]
%
% output	rr	- R_rho = alpha(p,t,s)*dt/beta(p,t,s)/ds
%
% uses : alpha, alpharho, betarho 
% 
% version 1.0.0		last change 1.8.1995

% G.Krahmann, IfM Kiel, Aug 1995

al=alpharho(p,t,s);
be=betarho(p,t,s);
rr=al.*dt./be./ds;
