function sig = sigmaany(p,t,s,pref)
%	function sig=sigmaany(p,t,s,[pref])
%	
% please use 'sigma' !!!!!
%
% input  :	p		: pressure [dbar]
%		t		: in situ temperature [degrees Celsius]
%		s		: salinity [psu]
%		pref	[p]	: optional reference pressure
%
% output :	sig		: density of seawater at pressure P (adiabatic)
%
% uses :	alpha.m
%
% version 1.1.1		last change 03.07.1996

% modified from SIGMATH, Uwe Send, March 1995
% optional without Pref		G.Krahmann, IfM Kiel, Sep 1995
% removed bug			G.Krahmann, IfM Kiel, Jul 1996

disp('WARNING : sigmaany.m is replaced by sigma.m ')
if nargin==4
  sig=sigma(p,t,s,pref);
else
  sig=sigma(p,t,s,p);
end
