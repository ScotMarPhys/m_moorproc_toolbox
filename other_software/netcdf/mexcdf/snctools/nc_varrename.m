function nc_varrename ( ncfile, old_variable_name, new_variable_name )
% NC_VARRENAME:  renames a NetCDF variable.
%
% NC_VARRENAME(NCFILE,OLD_VARNAME,NEW_VARNAME) renames a netCDF variable from
% OLD_VARNAME to NEW_VARNAME.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_varrename.m 2528 2008-11-03 23:06:25Z johnevans007 $
% $LastChangedDate: 2008-11-03 18:06:25 -0500 (Mon, 03 Nov 2008) $
% $LastChangedRevision: 2528 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nargchk(3,3,nargin);
nargoutchk(0,0,nargout);

if getpref('SNCTOOLS','USE_TMW',false)
    ncid=netcdf.open(ncfile,nc_write_mode);
    try
        netcdf.reDef(ncid);
        varid = netcdf.inqVarID(ncid, old_variable_name);
        netcdf.renameVar(ncid, varid, new_variable_name);
        netcdf.endDef(ncid);
        netcdf.close(ncid);
    catch myException
        netcdf.close(ncid);
        rethrow(myException);
    end

else
    [ncid,status ]=mexnc('OPEN',ncfile,nc_write_mode);
    if status ~= 0
        ncerr = mexnc('strerror', status);
        error ( 'SNCTOOLS:NC_VARGET:MEXNC:OPEN', ncerr );
    end
    
    
    status = mexnc('REDEF', ncid);
    if status ~= 0
        mexnc('close',ncid);
        ncerr = mexnc('strerror', status);
        error ( 'SNCTOOLS:NC_VARGET:MEXNC:REDEF', ncerr );
    end
    
    
    [varid, status] = mexnc('INQ_VARID', ncid, old_variable_name);
    if status ~= 0
        mexnc('close',ncid);
        ncerr = mexnc('strerror', status);
        error ( 'SNCTOOLS:NC_VARGET:MEXNC:INQ_VARID', ncerr );
    end
    
    
    status = mexnc('RENAME_VAR', ncid, varid, new_variable_name);
    if status ~= 0
        mexnc('close',ncid);
        ncerr = mexnc('strerror', status);
        error ( 'SNCTOOLS:NC_VARGET:MEXNC:RENAME_VAR', ncerr );
    end
    
    
    status = mexnc('ENDDEF', ncid);
    if status ~= 0
        mexnc('close',ncid);
        ncerr = mexnc('strerror', status);
        error ( 'SNCTOOLS:NC_VARGET:MEXNC:ENDDEF', ncerr );
    end
    
    
    status = mexnc('close',ncid);
    if status ~= 0
        ncerr = mexnc('strerror', status);
        error ( 'SNCTOOLS:NC_VARGET:MEXNC:CLOSE', ncerr );
    end
end


return;
