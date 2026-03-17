function test_varid ( ncfile )
% TEST_VARID
%

[ncid, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end


%
% DIMDEF
[xdimid, status] = mexnc ( 'def_dim', ncid, 'x', 20 );
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end
[ydimid, status] = mexnc ( 'def_dim', ncid, 'y', 24 );
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end
[zdimid, status] = mexnc ( 'def_dim', ncid, 'z', 32 );
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end


%
% VARDEF
[xdvarid, status] = mexnc ( 'def_var', ncid, 'x_double', 'double', 1, xdimid );
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  put_att_double failed on variable attribute, ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end



[varid, status] = mexnc('VARID', ncid, 'x_double');
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end



%
% Are they the same?
if ( varid ~= xdvarid )
	error ( 'VARID did not return the same varid as VARDEF..\n' );
end



%
% Bogus ncid
[status] = mexnc('VARID', -4, 'x_double');
if ( status >= 0 )
	error ( 'Bogus ncid case succeeded for VARID.\n' );
end

%
% Bogus varid
[status] = mexnc('VARID', ncid, 'newname');
if ( status >= 0 )
	error ( 'Bogus name case succeeded for VARID.\n' );
end


%
% test non existance ""
try
	[varid, status] = mexnc('VARID', ncid, '');
	msg = sprintf ( '%s:  INQ_VARID returned non-negative status for empty string variable\n', mfilename );
	error ( msg );
end


%
% ENDEF
[status] = mexnc ( 'enddef', ncid );
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end




status = mexnc ( 'close', ncid );
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  %s\n', mfilename, ncerr );
	error ( msg );
end


fprintf ( 1, 'VARID succeeded.\n' );

return













