function o2=oxsat(s,t)
%OXSAT Oxygen saturation.
%  O2 = OXSAT(S,T) return the oxygen saturation in ml/l as a function of
%  salinity s (PSS-78) and in situ temperature t (°C).
%
%  References:
%  Weiss, R. F. (1970). The solubility of nitrogen, oxygen and argon in water
%    and seawater. Deep-Sea Res., 17, 721-735.

%  C. Mertens, IfM Kiel
%  $Revision: 1.2 $ $Date: 1996/01/05 13:29:30 $


a = [-173.4292 249.6339 143.3483 -21.8492];
b = [-0.033096 0.014259 -0.0017]; 

t = 0.01*(t + 273.15);
o2 = a(1) + a(2)./t + a(3)*log(t) + a(4)*t + s.*(b(1) + b(2)*t + b(3)*t.*t);
o2 = exp(o2);

