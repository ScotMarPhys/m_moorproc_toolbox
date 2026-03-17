function test_inquire ( ncfile )
% TEST_INQUIRE
%
% Tests number of dimensions, variables, global attributes, record dimension for
% foo.nc
%
% Tests bad ncid as well.

[ncid, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status < 0 )
	error ( 'CREATE failed' );
end


%
% DIMDEF
[xdimid, status] = mexnc ( 'def_dim', ncid, 'x', 20 );
if ( status < 0 )
	error ( 'DEF_DIM failed on X' );
end
[ydimid, status] = mexnc ( 'def_dim', ncid, 'y', 24 );
if ( status < 0 )
	error ( 'DEF_DIM failed on y' );
end
[zdimid, status] = mexnc ( 'def_dim', ncid, 'z', 32 );
if ( status < 0 )
	error ( 'DEF_DIM failed on z' );
end


%
% VARDEF
[xdvarid, status] = mexnc ( 'def_var', ncid, 'x_double', 'double', 1, xdimid );
if ( status < 0 )
	error ( 'DEF_VAR failed on x_double' );
end


%
% Define some attributes
attvalue = 'this is a test';
attlen = length(attvalue);
status = mexnc ( 'put_att_text', ncid, xdvarid, 'test_variable_attributes', 'char', attlen, attvalue );
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  put_att_double failed on variable attribute, ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end
attvalue = 'this is a global test';
attlen = length(attvalue);
status = mexnc ( 'put_att_text', ncid, -1, 'test_global_attributes', 'char', attlen, attvalue );
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	err_msg = sprintf ( '%s:  put_att_double failed on variable attribute, ''%s''\n', mfilename, ncerr );
	error ( err_msg );
end


%
% ENDEF
[status] = mexnc ( 'enddef', ncid );
if ( status < 0 )
	error ( 'ENDEF failed with write' );
end



[ndims, nvars, natts, recdim, status] = mexnc('INQUIRE', ncid);
if ( status < 0 )
	error ( 'INQUIRE failed on nowrite' );
end

if ndims ~= 3
	msg = sprintf ( 'INQUIRE returned %d dimensions when there should only have been 1\n', ndims );
	error ( msg );
end

if nvars ~= 1 
	msg = sprintf ( 'INQUIRE returned %d variables when there should have been 18\n', nvars );
	error ( msg );
end

if natts ~= 1
	msg = sprintf ( 'INQUIRE returned %d attributes when there should have been 1\n', natts );
	error ( msg );
end

if recdim ~= -1
	msg = sprintf ( 'INQUIRE returned an unlimited dimension when there should not have been one\n' );
	error ( msg );
end



%
% Try a bogus case
[ndims, nvars, natts, recdim, status] = mexnc('INQUIRE', -1);
if ( status >= 0 )
	error ( 'INQUIRE return status did not signal an error on bogus ncid case' );
end
fprintf ( 1, 'INQUIRE succeeded\n' );





status = mexnc ( 'close', ncid );
if ( status < 0 )
	ncerr = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  %s\n', mfilename, ncerr );
	error ( msg );
end

return











