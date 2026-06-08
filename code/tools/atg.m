function a = atg(p,t,s)
%ATG Computes adiabatic temperature gradient (required by THETA).
%       A = ATG(P,T,S)

%       VAX 11/750      1983    J.HOLTORFF
%       18/02/93, C. Mertens, IfM Kiel, changed to Matlab

s = s-35.0 ;

a = (((-2.1687E-13*t + 1.8676E-11).*t - 4.6206E-10).*p ...
   + (( 2.7759E-10*t - 1.1351E-08).*s ...
   + ((-5.4481E-12*t + 8.7330E-10).*t - 6.7795E-08).*t + 1.8741E-06)).*p ...
   +  (-4.2393E-07*t + 1.8932E-05).*s ...
   + (( 6.6228E-09*t - 6.8360E-07).*t + 8.5258E-05).*t + 3.5803E-04 ;


