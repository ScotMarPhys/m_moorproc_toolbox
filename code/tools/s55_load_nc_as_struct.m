function DS = s55_load_nc_as_struct(filename,varargin)
% s55_load_nc_as_struct Loads NetCDF data into a MATLAB structure.
    if nargin==1
        var_attrs = 1;
    else
        var_attrs = varargin{1};
    end

    % Load file schema details
    info = ncinfo(filename);
    DS = struct();
    
    % Load Variables and their Variable-Specific Attributes
    for i = 1:length(info.Variables)
        varName = info.Variables(i).Name;
        validVarName = matlab.lang.makeValidName(varName); % Safe field formatting
        
        % Read the actual data matrix/vector
        DS.(validVarName) = ncread(filename, varName);
        
        % --- TIME CONVERSION LOGIC START ---
        % Look for attributes to see if this is a time variable
        var_info = info.Variables(i);
        is_time_var = false;
        base_date_str = '';
        
        for k = 1:length(var_info.Attributes)
            % Check if the attribute name is 'units'
            if strcmpi(var_info.Attributes(k).Name, 'units')
                units_val = var_info.Attributes(k).Value;
                % Match pattern like "days since 1950-01-01" or "hours since..."
                if contains(lower(units_val), 'since')
                    is_time_var = true;
                    % Extract everything after the word "since "
                    split_str = regexp(units_val, 'since\s+(.*)', 'tokens');
                    if ~isempty(split_str)
                        base_date_str = split_str{1}{1};
                    end
                end
            end
        end
        
        % If it is a time variable, convert the loaded data into MATLAB datenum
        if is_time_var && ~isempty(base_date_str)
            try
                % Strip out timezone text like "UTC" if present so datenum can read it
                clean_date_str = regexprep(base_date_str, '\s*UTC\s*', '');
                basetime = datenum(clean_date_str);
                
                % Convert the raw NetCDF data to MATLAB datenum
                DS.(validVarName) = DS.(validVarName) + basetime;
            catch
                warning('Could not parse time units string: %s', units_val);
            end
        end
        % --- TIME CONVERSION LOGIC END ---
        
        if var_attrs
            % Extract individual attributes for this specific variable
            if ~isempty(var_info.Attributes)
                attrStructName = [validVarName, '_attrs'];
                DS.(attrStructName) = struct();
                
                % Loop through all attributes belonging to this variable
                for j = 1:length(var_info.Attributes)
                    attName = var_info.Attributes(j).Name;
                    validAttName = matlab.lang.makeValidName(attName);
                    
                    DS.(attrStructName).(validAttName) = var_info.Attributes(j).Value;
                end
            end
        end
    end
    
    % Load Global Attributes (Metadata)
    DS.Attributes = struct();
    for i = 1:length(info.Attributes)
        attName = info.Attributes(i).Name;
        validGlobalAttName = matlab.lang.makeValidName(attName);
        DS.Attributes.(validGlobalAttName) = info.Attributes(i).Value;
    end
end
