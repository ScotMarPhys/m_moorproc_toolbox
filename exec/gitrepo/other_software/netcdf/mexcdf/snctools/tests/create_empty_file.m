function create_empty_file ( ncfile )
% CREATE_EMPTY_FILE:  Does just that, makes an empty netcdf file.
%
% USAGE:  create_empty_file ( ncfile );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: create_empty_file.m 2515 2008-07-03 20:36:38Z johnevans007 $
% $LastChangedDate: 2008-07-03 16:36:38 -0400 (Thu, 03 Jul 2008) $
% $LastChangedRevision: 2515 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if getpref('SNCTOOLS','USE_TMW',false)
	ncid_1 = netcdf.create(ncfile, nc_clobber_mode );
	netcdf.close(ncid_1);
else
	[ncid_1, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
	if ( status ~= 0 )
		ncerr_msg = mexnc ( 'strerror', status );
		msg = sprintf ( '%s:  ''create'' failed, error message '' %s ''\n', mfilename, ncerr_msg );
		error ( msg );
	end
	
	%
	% CLOSE
	status = mexnc ( 'close', ncid_1 );
	if ( status ~= 0 )
		error ( 'CLOSE failed' );
	end
end
return
