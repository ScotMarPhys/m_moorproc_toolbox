% RO Oceanographic functions library
%
% last change: 09.08.2000
%
% acoef.m               : used by sal78
% alpha.m               : equation of state for seawater, bulk-mod, ster-an
% alpharho.m            : thermal expansion coefficient
% alr.m                 : adiabatic lapse rate
% aou.m                 : Apparent oxygen utilization rate
% astp.m                : sigma-stp from t,s,p,sigma-t
% atg.m	                : adiabatic temperature gradient
% bcoef.m               : used by sal78
% betarho.m             : haline expansion coefficient
% buoyfreq              : buoyancy frequency
% ccoef.m               : used by sal78
% cpsw.m                : specific heat
% cepe.m                : specific heat of sea-water (OBSOLETE, use cpsw.m)
% dsal.m                : used by sal78
% getg.m                : acceleration due to gravity in m/s^2
% gravity               : Gravitational acceleration
% oxsat.m               : oxygen equilibrium saturation in ml/l 
%                         (after Weiss, 1970)
% oxsat_umol.m          : oxygen equilibrium saturation in umol/kg 
%                         (after Benson & Krause, 1984)
% oxsat_percent.m       : oxygen saturation in percent, 
%                         based on Benson & Krause, 1984
% qgmodes.m             : calculate quasigeostrophic modes
% rrho.m                : density gradient
% rt35.m                : used by sal78
% sal.m			: used by sal78
% sal78.m               : convert salinity to conductivity and back
% sigma.m               : density of seawater
% sigmaany.m            : density of seatwater (OBSOLETE, replaced by sigma.m)
% sigmat.m              : density of seatwater (OBSOLETE, replaced by sigma.m
% sigmath.m             : density of seatwater (OBSOLETE, replaced by sigma.m)
% sspcm.m               : sound speed in sea water
% stabang.m             : stability angle
% t68                   : Convert temperature [deg C] fom IPTS-68 
%                         temperature scale to IPTS-90
% t90                   : Convert temperature [deg C] from  IPTS-90 temperature 
%                         scale to ITS-68
% theta.m               : potential temperature
% thickness             : Calculate the thickness of a density layer

