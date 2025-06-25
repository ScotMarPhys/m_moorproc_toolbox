function d = read_cnv_data(filename,start_date_cast)

fid = fopen(filename, 'r');

if fid == -1
    error('Cannot open file: %s', filename);
end

% Step 1: Read header and parse variable names
varNames = {};
while ~feof(fid)
    line = fgetl(fid);
    if startsWith(line, '# name')
        % Parse line like: "# name 0 = prDM: Pressure, Digiquartz [db]"
        tokens = regexp(line, '# name \d+ = ([^:]+):', 'tokens');
        if ~isempty(tokens)
            varName = strtrim(tokens{1}{1});
            % Clean variable name to make it a valid struct field
            varName = matlab.lang.makeValidName(varName);
            varNames{end+1} = varName;
        end
    elseif contains(line, '*END*')
        break;
    end
end

% Step 2: Read numeric data
data = [];
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || all(isnan(str2double(strsplit(line))))
        continue;
    end
    values = str2double(strsplit(line));
    data = [data; values];
end

fclose(fid);

% Step 3: Create struct d with field names from varNames
d = struct();
numVars = min(length(varNames), size(data, 2));  % In case of mismatch
for i = 1:numVars
    d.(varNames{i}) = data(:, i);
end

d.cond=d.c0S_m*10;
d.cond1=d.c0S_m*10;
d.cond2=d.c1S_m*10;
d.temp=d.t090C;
d.temp1=d.t090C;
d.temp2=d.t190C;
d.press=d.prDM;
d.oxygen=d.sbeox0Mm_L;

d.mtime=repmat(start_date_cast,[length(d.timeS),1]);

% Convert to date vector
dv = datevec(d.mtime);

% Extract year, month, day
year = dv(:,1);
month = dv(:,2);
day = dv(:,3);

% Compute decimal hour
decimalHour = dv(:,4) + dv(:,5)/60 + dv(:,6)/3600;


% convert seconds to decimal hour
dechour=d.timeS/3600;

dechour=decimalHour+dechour;

% Display results
d.yd=julian(year, month, day, dechour);

fields = {'timeS', 'prDM', 'depSM', 't090C', 't190C', 'c0S_m',...
    'c1S_m', 'sal00', 'sal11', 'sbeox0V', 'sbeox0Mm_L', 'v0', 'v1',...
    'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'pumps', 'flag'};
d = rmfield(d,fields);

end