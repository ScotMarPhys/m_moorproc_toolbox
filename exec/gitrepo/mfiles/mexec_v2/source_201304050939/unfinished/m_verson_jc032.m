function m_verson(ncfile,increment)


% allow user to specify version file increment. Can specify as zero when
% loading pstar files to mstar.

if nargin == 1; increment = 1; end

m_common

MEXEC.versfile = MEXEC_G.VERSION_FILE;  %at the moment, the version file is a mat file

load(MEXEC.versfile); %contains 'datanames' which is a cell array  and 'versions' which is a double array
n = length(datanames);



file_dataname = nc_attget(ncfile.name,nc_global,'dataname');
file_version = nc_attget(ncfile.name,nc_global,'version');


kmatch = strmatch(file_dataname,datanames,'exact');

if length(kmatch) > 1
    error(['problem with multiple occurrences of dataname ' file_dataname ' in version file'])
end

if isempty(kmatch) %new dataname at this site
    index = n+1;
    datanames{index} = file_dataname;
%     versions(index) = file_version;
    new_version = 1;
else
    index = kmatch;
%     new_version = increment+max(versions(index),file_version);
    new_version = increment+versions(index);
end
% keyboard
versions(index) = new_version;
save(MEXEC.versfile,'datanames','versions');
    
nc_attput(ncfile.name,nc_global,'version',new_version);
nc_attput(ncfile.name,nc_global,'mstar_site',MEXEC_G.SITE);

return