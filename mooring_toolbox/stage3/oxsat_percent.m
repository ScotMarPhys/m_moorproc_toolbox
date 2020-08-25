function ox_pc = oxsat_percent(pt,s,o)

% OXSAT_PERCENT    Calculation of oxygen saturation in percent.
%
%                              C_sea * 100    with C_sat   : saturation
%                  C_sat [%] = -----------         C_sea   : o2-concentration 
%                                C_equi                      [umol/kg]
%                                                  C__equi : oxygen equilibrium
%                                                            concentration
%                                                            [umol/kg]
%
%                  usage  : ox_pc = oxsat_umol(pt,s,o)
%
%                  input  : pt      potential temperature [deg C]
%                           s       salinity              [psu]
%                           o       oxygen concentration  [umol/kg]
%
%                  output : ox_pc   oxygen saturation     [%]   
%
%                  uses   : oxsat_umol.m
%
%                  version 1.0, Wed Jun 21 23:08:28 MET DST 2000, d.kieke
%

% --- caluclate oxygen equilibrium concentration ---

o_equi = oxsat_umol(pt,s);

% --- saturation in % ---

ox_pc = (o*100)./o_equi;



