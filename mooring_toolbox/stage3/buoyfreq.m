function [n2,n2m,pm]=buoyfreq(p,t,s,lat)
% function [n2,n2m,pm]=buoyfreq(p,t,s,lat)
%
% calculate Buoyancy Frequency N^2
% 
% input  :	p		: pressure [dbar]
%		t		: in situ temperature [degrees Celsius]
%		s		: salinity [psu]
%
% output :	n2		: buoyancy frequency^2 [s^(-2)]
%				  calculated with central differences
%				  onto all grid points (at the end non-central)
%		n2m		: buoyancy frequency^2 [s^(-2)]
%				  calculated between grid points
%		pm		: resulting new grid
%
% uses :	p2z, getg, alpha, theta
%
% version 1.1.0		last change 05.09.1995

% G.Krahmann, IfM Kiel, Sep 1995
% changed formula to adiabatic leveling method (Millard DSR 37, 167-181)
%	G.Krahmann, IfM Kiel, Sep 1995

% transpose if necessary
if size(p,1)==1
  p=p';
  t=t';
  s=s';
  trans=1;
else
  trans=0;
end

% get depth from pressure
z=p2z(p,lat);

% earth acceleration
g=getg(lat,z);

% begin old formula after Pond and Pickard
if 0
rho=1000+sigma(p,t,s);
ss=sspcm(p,t,s);

dummy1=imag(gradient(rho));
dummy2=imag(gradient(-z));

n2=g.*( -1./rho.*dummy1./dummy2 -g./ss.^2 );
return
end
% end old formula

% old grid == original p points
% new grid == p-averages between following points
% prepare averages between old grid
pu=p(1:size(p,1)-1,:);
pl=p(2:size(p,1),:);
su=s(1:size(p,1)-1,:);
sl=s(2:size(p,1),:);
tu=t(1:size(p,1)-1,:);
tl=t(2:size(p,1),:);

% prepare new grid
pm=(pu+pl)/2;
tm=(tu+tl)/2;
sm=(su+sl)/2;
gm=getg(lat,p2z(pm,lat));

% calculate potential steric anomalies
[d1,d2,du3]=alpha1(pm,theta(pu,tu,su,pm),su);
[d1,d2,dl3]=alpha1(pm,theta(pl,tl,sl,pm),sl);

% calculate n2 on new and old grid
n2m=(du3-dl3)./(pu-pl);
if size(p,1)==2
  n2=[n2m(1,:);n2m(1,:)];
else
  n2=[n2m(1,:);(n2m(1:size(n2m,1)-1,:)+...
	n2m(2:size(n2m,1),:))/2;n2m(size(n2m,1),:)];
end
n2m=-gm.^2./alpha1(pm,tm,sm).^2.*n2m/1e4;
n2=-g.^2./alpha1(p,t,s).^2.*n2/1e4;

if trans==1
  n2=n2';
end
