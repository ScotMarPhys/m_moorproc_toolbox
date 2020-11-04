function dim_names = m_unpack_dimnames(ncfile)

% Unpack dimension names from netcdf file

%If the metadata isn't passed in, then read it from the file
toc
if ~isfield(ncfile,'metadata')
    metadata = nc_info(ncfile.name);
    ncfile.metadata = metadata;
end
toc
metadata = ncfile.metadata;
toc
if isfield (metadata, 'Dimension' )
    num_dims = length(metadata.Dimension);
else
    num_dims = 0;
end
toc
dim_names = cell(num_dims,1);

toc
for k = 1:num_dims
    dim_names{k} =  metadata.Dimension(k).Name;
end
toc
return