function get_loc_depth(varargin)
	m_setup
% Get position and corrected water depth at time or times
	if nargin== 2
		t1 = datenum(cell2mat(varargin(1)));
		t2 = datenum(cell2mat(varargin(2)));
	elseif nargin == 1
		t1 = datenum(cell2mat(varargin(1)))-30/86400;		
		t2 = t1 + 60/86400;		
	else
		t1 = now -60/86400;
		t2 = now;
	end
	
nav = mtload('posmvpos',t1,t2,'time long lat');
ea600 = mtload('sim',t1,t2,'time depthm');
%em120 = mtload('em120',t1,t2,'time snd');

if length(nav) > 0
    tme = nanmean(nav.time) + MEXEC_G.uway_torg;;
    lat = nanmean(nav.lat);
    lon = nanmean(nav.long);
	if length(ea600) > 0
		if sum(~isnan(ea600.depthm)) >= 1
			dep = nanmean(ea600.depthm);
     	   	corr_struct = mcarter(lat,lon,dep);
     	  	corr_depth = corr_struct.cordep;
  	  	else
    		corr_depth = NaN;
  	  	end
	else
		corr_depth = NaN;
	end
else
	tme = NaN;
	lat = NaN;
	lon = NaN;
end

fprintf(1,'\n %s \n lat = %7.4f lon = %7.4f \n Corr depth = %6.0fm \n',datestr(tme),lat,lon,corr_depth)
[latdd, latmm] = dd2dm2(lat);
[londd, lonmm] = dd2dm2(lon);
fprintf(1,' %4.0f %5.2f N %4.0f %5.2f E  \n',latdd,latmm,londd,lonmm)