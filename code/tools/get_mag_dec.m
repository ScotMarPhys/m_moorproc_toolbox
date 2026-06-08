%function mag_dec=get_mag_dec(date_in,lat,lon)
%
% Written for paths on D382, but can be easily updated to give path to
% geomag70.exe linux executable
%
% function to call the linux routine geomag70.exe that queries the
% IGRF11.COF database (released 2009) to get magnetic information for the
% position(s) entered. This version of geomag will work with dates up to
% 2015. http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
%
% Required inputs: date in [yy mm dd] format - can be multiple rows
%                  lat - latitude in decimal degrees N
%                  lon - longitude in decimal degrees E
% Output: mag_dec array the same length as the inputs with the magnetic
%         declination calcuclated by the geomag70.exe executable.

function mag_dec=get_mag_dec(date_in,lat,lon)
    current_dir=pwd;
    if (length(date_in(:,1))~=length(lat) || length(date_in(:,1))~=length(lon) || length(lat)~=length(lon))
        disp('inputs to get_mag_dec are not the same length. Stopping!')
        return
    end
    cd /noc/users/pstar/rpdmoc/rapid/data/exec/d382/tools/geomag70_linux/
    for i=1:length(date_in(:,1))
        seconds=DMYhhmmss2sec([date_in(i,3),date_in(i,2),date_in(i,1),zeros(size(date_in(i,1:3)))]);
        dec_year=seconds/60/60/24/365.25+date_in(i,1);
        [status,output]=eval(['system(''geomag70.exe IGRF11.COF ' num2str(dec_year) ' D M0 ' num2str(lat(i)) ' ' num2str(lon(i)) ''')']);
        mag_dec_index=strfind(output,'(yr)      (deg min)   (deg min)     (nT)     (nT)     (nT)     (nT)     (nT)');
        mag_dec_string=output(mag_dec_index+90:mag_dec_index+97);
        if strmatch(mag_dec_string(1),'-')
            mag_dec(i)=str2double(mag_dec_string(1:strfind(mag_dec_string,'d')-1))-str2double(mag_dec_string(strfind(mag_dec_string,'d')+1:strfind(mag_dec_string,'m')-1))/60;
        else
            mag_dec(i)=str2double(mag_dec_string(1:strfind(mag_dec_string,'d')-1))+str2double(mag_dec_string(strfind(mag_dec_string,'d')+1:strfind(mag_dec_string,'m')-1))/60;
        end
    end
    cd(current_dir)
end