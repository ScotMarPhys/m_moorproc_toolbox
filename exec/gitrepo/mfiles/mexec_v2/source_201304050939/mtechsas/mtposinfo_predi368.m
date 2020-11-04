function [lat lon] = mtposinfo(dn1,navstream)
% function [lat lon] = mtposinfo(varargin)
% eg
% [lat lon] = mtposinfo;
% [lat lon] = mtposinfo([2009 4 4 12 0 0]);
% [lat lon] = mtposinfo([2009 4 4 12 0 0],'dps116');
% mtposinfo now-1;
% mtposinfo '2009 4 4 12 0 0';
% mtposinfo now dps116;
% mtposinfo;
% 
% Obtain position from a techsas nav file. On JC032 this was hardwired to
% posmvpos. If no argument, the latest data cycle in the latest file is returned.
% The argument can be a Matlab datenum or a Matlab datevec
% 
% 
% first draft by BAK on JC032

m_common

if ~exist('navstream','var'); navstream = MEXEC_G.Mtechsas_default_navstream; end  % use default nav stream name
instream = navstream; % mexec stream short name
tstream = mtresolve_stream(instream);

if ~exist('dn1','var')
    [pdata u] = mtlast(tstream);
    lat = pdata.lat;
    lon = pdata.long;
    dn = pdata.time+MEXEC_G.Mtechsas_torg;
elseif isempty(dn1); 
    [pdata u] = mtlast(tstream);
    lat = pdata.lat;
    lon = pdata.long;
    dn = pdata.time+MEXEC_G.Mtechsas_torg;
else
    if ischar(dn1);
        cmd =['dn1 = [' dn1 '];'];  % if the arg has come in as a string, convert from char to number
        eval(cmd);
    end
    dn = datenum(dn1);
    % load data 5 minutes either side in case the required time falls precisely
    % between two files

    pdata = mtload(tstream,dn-5/1440,dn+5/1440,'time lat long','q');

    tin = pdata.time+MEXEC_G.Mtechsas_torg;
    latin = pdata.lat;
    lonin = pdata.long;
    % jbs on jc054:
    [uni iunique]=unique(tin); 
    lat = interp1(tin(iunique),latin(iunique),dn);
    lon = interp1(tin(iunique),lonin(iunique),dn);
% else
%     m = 'expect precisely zero or one input arg which should be matlab datenum or datvec)';
%     fprintf(MEXEC_A.Mfider,'%s\n',m)
%     lat = nan; lon = nan;
%     return
end

if nargout > 0 return; end
% else print to screen


[latd latm] = m_degmin_from_decdeg(lat);
[lond lonm] = m_degmin_from_decdeg(lon);

dvnow = datevec(now);
yyyy = dvnow(1);
doffset = datenum([yyyy 1 1 0 0 0]);
daynum1 = floor(dn) - doffset + 1;

str1 = datestr(dn,'yy/mm/dd');
str1a = datestr(dn,'HH:MM:SS');
fprintf(MEXEC_A.Mfidterm,'%s\n',tstream);
fprintf(MEXEC_A.Mfidterm,'%s     %8s   %03d %8s\n','time',str1,daynum1,str1a);
fprintf(MEXEC_A.Mfidterm,'%s %10.5f %5.0f %6.2f\n','lat',lat,latd,latm);
fprintf(MEXEC_A.Mfidterm,'%s %10.5f %5.0f %6.2f\n','lon',lon,lond,lonm);

return

