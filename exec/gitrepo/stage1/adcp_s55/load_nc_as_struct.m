function DS = load_nc_as_struct(filename)
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