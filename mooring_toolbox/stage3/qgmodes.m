function [uvmode,dispmode,rhomode,z ]=qgmodes(z,n2,modes);
% [uvmode,dispmode,rhomode,z] = qgmodes(z,n2,modes);
%
% QGMODES  is a function to calculate quasigeostrophic modes
%          for a given N^2 Profile.
%
% input  : 	z         	: depth vector (must be equidistant !!)
%		n2		: buoyancy frequency profile
%          	modes     	: number of modes to be calculated
%
% output :     	uvmodes   	: velocity modes 
%          	dispmodes 	: displacement modes 
%          	rhomode   	: density pertubation modes
%
% version 1.1.1		last change 20.12.1995

%          translated from FORTRAN to Matlab 3.5
%          by Ulf Garternicht in january 1994 for IfM Kiel
% changed input and header		G.Krahmann, IfM Kiel, Sep 1995
% corrected header, better code		G.Krahmann, IfM Kiel, Dec 1995

% calculate eigenfunctions of tridiagonal matrix

A=n2.^(-1);
dz=z(2)-z(1);
r=1/dz/dz;
len=length(A);
B=r.*(diag(2*A(2:len-1))+diag(-A(3:len-1),-1)+diag(-A(2:len-2),1));
[eigfct,eigval]=eig(B);
eigval=diag(eigval);

[tmp n]=sort(eigval);
eigfct=eigfct(:,n(1:modes));
eigval=eigval(n(1:modes)); 
% normalize displacement modes with relation
% dispmodes(:,i)'*(dispmodes(:,j).*n2) = delta(i,j) ,
% velocity modes with  d dispmode(z) /dz + eigval*uvmode(z) = 0

[m n]=size(eigfct);
ortho=eigfct'*(eigfct.*(n2(2:len-1)*ones(1,n)));
ortho=diag(ortho);
dispmode=-eigfct./(ones(m,1)*sqrt(ortho(:)'))*sqrt(m+1);
dispmode=[zeros(1,n);dispmode;zeros(1,n)];

len=length(z);
zz=z(1:len-1)+dz/2;
%%umode=-diff(dispmode)/dz;
%%uvmode=table1f([zz,umode],z,1);
uvmode=-gradient(dispmode')'/dz;
ortho=diag(uvmode'*uvmode);
uvmode=uvmode./(ones(m+2,1)*sqrt(ortho(:)'))*sqrt(m+1);
rhomode=dispmode.*(n2*ones(1,n));
