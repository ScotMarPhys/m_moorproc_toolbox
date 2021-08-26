function tc = t68(t)

% T68    tc = t68(t)
%
%        Convert temperature [deg C] from ITS-90 temperature scale to IPTS-68.
%        This is needed in computation of properties of seawater.
%
%        input  : t   in situ temperature ITS-90 [deg C]
%
%        output : tc  converted in situ temperature IPTS-68 [deg C]
%
%        uses   : ---
%
%
%
%  C. Mertens, IfM Kiel
%  $Revision: 1.2 $ $Date: 1996/01/09 11:51:40 $
%
%  04.03.1998, d.kieke, improved header

tc = 1.00024*t;

