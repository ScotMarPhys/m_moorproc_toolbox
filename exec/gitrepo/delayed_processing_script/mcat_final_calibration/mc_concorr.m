% function ct = mc_concorr(c,pac,[pref])
%
% correct microcat conductivity  (without pressure sensor)
%  
% If microCATs are deployed without build-in pressure sensors
% a reference pressure 'pref' has to be chosen according to the
% depth where the pressure sensor is moored. Conductivity itself
% is pressure dependend for it depends on the geometric layout
% of the conductivity cell. The latter is changed due to changing
% pressures. This function provides a correction that accounts for
% this effect (according to Seabird Electonics Application Note 10).
%
% input:
%         c     ----  microCAT conductivity [any unit]  
%         pac   ----  actual pressure [dbar]
%         pref  ----  setup reference pressure [dbar]
%                     if omitted, pref = 0 dbar  
%
% output:  
%         ct    ----  true conductivity [same unit as c]
%
% Kanzow 22.03.01
 
function ct = mc_concorr(c,p,pref)

if nargin == 2
  pref = 0;
end
   
spref = 1 - 9.57e-8 * pref;

sp    = 1 - 9.57e-8 * p;

ct    = c .* spref ./ sp;

  
