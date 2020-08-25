function apoxu = aou(pt,s,o)

% AOU          Apparent oxygen utilization, defined as the difference
%              between the observed oxygen content and the saturation oxygen 
%              content of a sample of sea water. This is a method of estimating 
%              the amount of dissolved oxygen utilized by organisms via 
%              respiration, although it is called "apparent" for a reason. 
%              Surface waters may more than likely carry more than the 
%              saturation amount of oxygen due to the nonlinearity in the 
%              solubility of oxygen with temperature. The effects of
%              this nonlinearity are small, though, and the AOU is usually 
%              quite close to TOU, the True Oxygen Utilization. 
%              See Broecker and Peng (1982).
%
%              (from Steve Baum's glossary...)'   
%
%              usage  : apoxu = aou(pt,s,o);
%
%              input  : pt    potential temperature  [deg C]
%                       s     salinity               [psu]
%                       o     oxygen concentration   [umol/kg]                 
%
%              output : aou   apparent oxygen utilization [umol/kg]
%
%              uses   : oxsat_umol.m
%
%              
%              version 1.0, Wed Jun 21 19:15:56 MET DST 2000, d.kieke


o_equi = oxsat_umol(pt,s);
apoxu  = o - o_equi;



