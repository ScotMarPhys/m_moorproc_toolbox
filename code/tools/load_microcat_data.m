function mc_data = load_microcat_data(mc_infile)
% LOAD_MICROCAT_DATA Loads MicroCAT variables into a structure using rodbload
%
% USAGE:
%   mc_data = load_microcat_data(mc_infile)

    % Define the variable list for rodbload
    mc_vars = 'yy:mm:dd:hh:t:c:p:SerialNumber:Instrdepth';

    if exist(mc_infile, 'file');
        % Load variables from the file
        [yy, mm, dd, hh, t, c, p, sn_mc, depth_mc] = rodbload(mc_infile, mc_vars);

        dum = -9999.0000;
       
        % 1. Apply QC immediately after loading
        t      = dum2nan(t,dum);
        c      = dum2nan(c,dum);
        p      = dum2nan(p,dum);

        % 2. Calculate Salinity 
        SP = gsw_SP_from_C(c, t, p);     % Practical salinity from conductivity
        SP(SP < 0) = NaN; % Additional QC for salinity
       
        % Assign remaining data to structure
        mc_data.time = datenum(yy, mm, dd, hh, 0, 0);
        mc_data.t  = t;
        mc_data.c  = c;
        mc_data.p  = p;
        mc_data.SP   = SP;
        mc_data.sn = sn_mc;
        mc_data.z  = depth_mc;


    else
        warning('MicroCAT file not found: %s. Returning empty structure.', mc_infile);
        mc_data = struct('yy',[],'mm',[],'dd',[],'hh',[],'t',[],'c',[],'p',[],'sn',[],'z',[]);
    end
end