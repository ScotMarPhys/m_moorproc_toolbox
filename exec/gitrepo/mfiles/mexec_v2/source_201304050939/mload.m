function [d h] = mload(varargin)

%sourceroot = '/local/users/pstar/cruise/'; % parent of /local/users/pstar/cruise/sw/mexec_v2/source

setpref('MEXNC','USE_TMW',true); % use some matlab builtin netcdf functions
setpref('SNCTOOLS','USE_TMW',true);

%addpath([sourceroot '/sw/mexec_v2/source/msubs']);
%addpath([sourceroot '/sw/mexec_v2/source/unfinished']);
%addpath([sourceroot '/sw/mexec_v2/source/msubs_old']);
%addpath([sourceroot '/sw/netcdf/mexcdf/snctools']);
%addpath([sourceroot '/sw/netcdf/mexcdf/mexnc']);

% load data and header contents of mstar NetCDF file into structure
% arrays
% eg
% [d h] = mload('filename.nc');
% or type mload to be prompted for answers

m_common
m_margslocal
m_varargs

MEXEC_A.Mprog = 'mload';
m_proghd;


% varg = varargin;nargi
% 
% if ~isempty(varg) > 0
%     v1 = varg{1};
%     varg(1) = [];
% else
%     v1 = {};
% end
% 

fn = m_getfilename; % this uses the optional input argument if there is one
ncfile.name = fn;
ncfile = m_ismstar(ncfile); %check it is an mstar file and that it is not open

h = m_read_header(ncfile);
m_print_header(h);
% keyboard
if nargout == 0
    fprintf(MEXEC_A.Mfider,'\n%s\n',' warning: data will be saved as variable ''ans'' in calling program unless called with at least one argument');
elseif nargout == 1
    fprintf(MEXEC_A.Mfider,'\n%s\n',' warning: header won''t be saved in calling program unless called with at least two arguments');
end

endflag = 0;
while endflag == 0
    if exist('d','var') == 1 & length(fieldnames(d)) == h.noflds; 
        m = 'All variables now loaded';
        fprintf(MEXEC_A.Mfidterm,'\n%s\n',m);
        endflag = 1; 
        continue; 
    end % h.noflds vars have been loaded so we assume that's all of them. No point asking for more names
    %if (length(varargin)<2)
    m = sprintf('%s\n','Type variable names or numbers to load (0 or return to finish, ''/'' for all):');
    %     if ~isempty(varg) > 0
    %         var = varg{1}; varg(1) = [];
    %     else
    var = m_getinput(m,'s');
    %     end
    %else 
    %var=varargin{2};
    %endflag=1;
    %end
    if strcmp(' ',var) == 1; endflag = 1; continue; end
    if strcmp('0',var) == 1; endflag = 1; continue; end

    vlist = m_getvlist(var,h);
    m = ['list is ' sprintf('%d ',vlist) ];
    disp(m);


    for kk = 1:length(vlist)
        vname = h.fldnam{vlist(kk)};
        vname2 = m_check_nc_varname(vname);
        cmd = ['d.' vname2 ' = nc_varget(ncfile.name,vname);'];
        eval(cmd);
        m = [sprintf('%15s',vname) ' loaded as d.' vname2];
        disp(m)
    end

end


