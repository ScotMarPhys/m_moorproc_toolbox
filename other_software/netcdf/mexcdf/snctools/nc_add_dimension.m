function nc_add_dimension ( ncfile, dimension_name, dimension_length )
% NC_ADD_DIMENSION:  adds a dimension to an existing netcdf file
%
% USAGE:  nc_add_dimension ( ncfile, dimension_name, dimension_size );
%
% PARAMETERS:
% Input:
%     ncfile:  path to netcdf file
%     dimension_name:  name of dimension to be added
%     dimension_size:  length of new dimension.  If zero, it will be an
%         unlimited dimension.
% Output:
%     none
%
% In case of an error, an exception is thrown.
%
% Because of underlying limitations, this m-file requires mexnc.
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_add_dimension.m 2528 2008-11-03 23:06:25Z johnevans007 $
% $LastChangedDate: 2008-11-03 18:06:25 -0500 (Mon, 03 Nov 2008) $
% $LastChangedRevision: 2528 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nargchk(3,3,nargin);

if getpref('SNCTOOLS','USE_TMW',false)
    nc_add_dimension_tmw ( ncfile, dimension_name, dimension_length )
else
    nc_add_dimension_mex ( ncfile, dimension_name, dimension_length )
end

return



%-----------------------------------------------------------------------
function nc_add_dimension_tmw ( ncfile, dimension_name, dimension_length )

ncid = netcdf.open(ncfile, nc_write_mode );

try
    netcdf.reDef(ncid );
    dimid = netcdf.defDim(ncid, dimension_name, dimension_length);
    netcdf.endDef(ncid );
    netcdf.close(ncid );
catch myException
    netcdf.close(ncid);
    rethrow(myException);
end



return







%-----------------------------------------------------------------------
function nc_add_dimension_mex ( ncfile, dimension_name, dimension_length )
[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if status
    ncerr = mexnc ( 'strerror', status );
    error_id = 'SNCTOOLS:NC_ADD_DIMENSION:openFailed';
    error ( error_id, ncerr );
end

status = mexnc ( 'redef', ncid );
if status
    mexnc ( 'close', ncid );
    ncerr = mexnc ( 'strerror', status );
    error_id = 'SNCTOOLS:NC_ADD_DIMENSION:redefFailed';
    error ( error_id, ncerr );
end

[dimid, status] = mexnc ( 'def_dim', ncid, dimension_name, dimension_length );
if status
    mexnc ( 'close', ncid );
    ncerr = mexnc ( 'strerror', status );
    error_id = 'SNCTOOLS:NC_ADD_DIMENSION:defdimFailed';
    error ( error_id, ncerr );
end

status = mexnc ( 'enddef', ncid );
if status
    mexnc ( 'close', ncid );
    ncerr = mexnc ( 'strerror', status );
    error_id = 'SNCTOOLS:NC_ADD_DIMENSION:enddefFailed';
    error ( error_id, ncerr );
end


status = mexnc ( 'close', ncid );
if status 
    ncerr = mexnc ( 'strerror', status );
    error_id = 'SNCTOOLS:NC_ADD_DIMENSION:closeFailed';
    error ( error_id, ncerr );
end



return












