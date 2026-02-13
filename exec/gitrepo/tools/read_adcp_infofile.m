function info_adcp = read_adcp_infofile(infofile)
% READ_ADCP_INFOFILE Reads ADCP metadata into a structured format
%
% USAGE:
%   info_adcp = read_adcp_infofile(infofile)

    % Define the variables to be loaded from the RODB file
    infovar = 'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:MagDeviation'; 

    if exist(infofile, 'file')
        % Load data using the rodbload utility
        [id, sn, z, s_t, s_d, e_t, e_d, lat, lon, wd,magdev] = rodbload(infofile, infovar);
    else
        warning('No info file found at: %s. Setting fields to NaN.', infofile);
        % Use deal to assign NaN to all outputs if file doesn't exist
        [id, sn, z, s_t, s_d, e_t, e_d, lat, lon, wd,magdev] = deal(NaN);
    end

    % Assign variables to the info_adcp structure
    info_adcp.id  = id;
    info_adcp.sn  = sn;
    info_adcp.z   = z;
    info_adcp.s_t = s_t;
    info_adcp.s_d = s_d;
    info_adcp.e_t = e_t;
    info_adcp.e_d = e_d;
    info_adcp.lat = lat;
    info_adcp.lon = lon;
    info_adcp.wd  = wd;
    info_adcp.magdev = magdev;

end