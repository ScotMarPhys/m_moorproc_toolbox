function [d h] = mloadoriginal(varargin)

% load data and header contents of mstar NetCDF file into structure
% arrays
% eg
% [d h] = mload('filename.nc');
% or type mload to be prompted for answers

m_common
MEXEC_A.Mprog = 'mload';
m_proghd;


varg = varargin;

MEXEC_A.MARGS_IN = [varg MEXEC_A.MARGS_IN];
varargin = {};

fn = m_getfilename(varargin{:}); % this uses the optional input argument if there is one
ncfile.name = fn;
ncfile = m_ismstar(ncfile); %check it is an mstar file and that it is not open

h = m_read_header(ncfile);
m_print_header(h);

if nargout == 0
    fprintf(fider,'\n%s\n',' data won''t be saved in calling program unless called with at least one argument');
elseif nargout == 1
    fprintf(fider,'\n%s\n',' header won''t be saved in calling program unless called with at least two arguments');
end

endflag = 0;
while endflag == 0
    m = sprintf('%s\n','Type variable names or numbers to load (0 or return to finish, ''/'' for all):');
    var = m_getinput(m,'s');
    if strcmp(' ',var) == 1; endflag = 1; continue; end
    if strcmp('0',var) == 1; endflag = 1; continue; end

    vlist = m_getvlist(var,h);
    m = ['list is ' sprintf('%d ',vlist) ];
    disp(m);


    for k = 1:length(vlist)
        vname = h.fldnam{vlist(k)};
        vname2 = m_check_nc_varname(vname);
        cmd = ['d.' vname2 ' = nc_varget(ncfile.name,vname);'];
        eval(cmd);
        m = [sprintf('%15s',vname) ' loaded as d.' vname2];
        disp(m)
    end

end


