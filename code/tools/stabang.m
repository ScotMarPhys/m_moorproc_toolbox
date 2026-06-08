function [sang,sangm,pm]=stabang(p,t,s);
% function [sangsangm,pm]=stabang(p,t,s);
%
% calculate stability angle profile after
% Washburn & Kaese, JPO 1987, 17, 12-25
%
% input  : 	p		: pressure profile [dbar]
%		t		: in situ temperature profile [degrees Celsius]
%		s		: salinity profile [psu]
%
% output :	sang		: stability angle profile [degrees]
%				  on original grid (central differences)
%		sangm		: stability angle profile [degrees]
%				  on new grid (between old points)
%		pm		: new pressure grid
%
% uses :	alpharho, betarho, alpha
%
% version 1.0.0		last change 04.09.1995

% G.Krahmann, IfM Kiel, Sep 1995

% check input and transpose if necessary
trans=0;
if size(p,1)<2
  if size(p,2)>0
    p=p';
    t=t';
    s=s';
    trans=1;
  else
    disp('A profile is needed for the calculation !')
    return
  end
end

% determine new grid
pm=p(1:size(p,1)-1,:)+diff(p)/2;
tm=t(1:size(t,1)-1,:)+diff(t)/2;
sm=s(1:size(s,1)-1,:)+diff(s)/2;

% calculate expansion coefficients
al=alpharho(p,t,s);
be=betarho(p,t,s);
alm=alpharho(pm,tm,sm);
bem=betarho(pm,tm,sm);

% calculate vertical gradients
tz=imag(gradient(t));
sz=imag(gradient(s));
tzm=diff(t);
szm=diff(s);

% calculate stability angle (negative arguments for z positive upward)
sang=atan2( -(al.*tz-be.*sz) , -(al.*tz+be.*sz) )/pi*180;
sangm=atan2( -(alm.*tzm-bem.*szm) , -(alm.*tzm+bem.*szm) )/pi*180;

% if necessary transpose back
if trans==1
  sang=sang';
  sangm=sangm';
  pm=pm';
end
