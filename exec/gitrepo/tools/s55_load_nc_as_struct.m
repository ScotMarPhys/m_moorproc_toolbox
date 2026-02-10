function DS = s55_load_nc_as_struct(filename)
% s55_load_nc_as_struct Loads NetCDF data into a MATLAB structure.
%   DS = s55_load_nc_as_struct(filename) reads all variables and global 
%   attributes from the specified NetCDF file (.nc) and returns them 
%   as fields in a dynamic structure (DS).
%
%   INPUT:
%       filename - String or char array containing the full path to 
%                  the .nc file.
%
%   OUTPUT:
%       DS - A structure containing:
%            - [VariableNames]: Each NetCDF variable is saved as a 
%              field (e.g., DS.velocity, DS.time).
%            - Attributes: A sub-structure containing all global 
%              metadata (e.g., DS.Attributes.mooring_id).
%
%   EXAMPLE:
%       data = s55_load_nc_as_struct('rhadcp_01_2020_stage1.nc');
%       plot(data.time, data.pressure);
%
%   See also NCREAD, NCINFO.

    % Load all variables and attributes into a structure
    info = ncinfo(filename);
    DS = struct();
    
    % Load Variables
    for i = 1:length(info.Variables)
        varName = info.Variables(i).Name;
        DS.(varName) = ncread(filename, varName);
    end
    
    % Load Global Attributes (Metadata)
    DS.Attributes = struct();
    for i = 1:length(info.Attributes)
        attName = info.Attributes(i).Name;
        DS.Attributes.(attName) = info.Attributes(i).Value;
    end
end