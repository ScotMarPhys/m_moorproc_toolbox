function tc = t90(t)

% T90    tc = t90(t)
%
%        Convert temperature [deg C] fom IPTS-68 temperature scale to ITS-90.
%        This should be used when reporting temperature measurements 
%
%        input  : t   in situ temperature IPTS-68 [deg C]
%
%        output : tc  converted in situ temperature ITS-90 [deg C]
%
%        uses   : ---
%
%
%
%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1996/01/09 11:48:55 $
%
%  04.03.1998, d.kieke, improved header


tc = 0.99976*t;

