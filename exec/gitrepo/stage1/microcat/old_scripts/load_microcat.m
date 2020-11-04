% load_microcat.m
% load a converted SBE 37 microcat data file and return
% Matlab serial date and data structure with whatever variables are there
% Adam Houk
% 4/17/2011
function [mc,cfg] = load_microcat(filename)
fid = fopen(filename,'r');

% first line should contain the type of SeaBird used
headerline = fgetl(fid);

if any(findstr(headerline,'SBE37'))
     cfg.instr_type = 'SBE37';
else
     disp('Wrong file-type!')
end

cfg.full_path =  fgetl(fid);
cfg.sw_version = cell2mat(regexp(fgetl(fid),'\d+\.?\d*','match'));
cfg.temp_sn = cell2mat(regexp(fgetl(fid),'\d+\.?\d*','match'));
cfg.cond_sn = cell2mat(regexp(fgetl(fid),'\d+\.?\d*','match'));
cfg.upload_time = char(regexp(fgetl(fid),'(?<=[\=]).*','match'));

i=1;
str = fgetl(fid);
while (~strncmp(str,'*END*',5));
	if strncmp(str,'* ds',4)
	% system header/status
		str = fgetl(fid);
        if length(str) > 4
			array1 = regexp(str,'\d+.?\S*','match');
            cfg.firmware_ver = array1{2};
            cfg.instrument_sn = array1{3};
        end
	  
		while ~strncmp(str,'* S>',4)
	       disp(str)
	       if ~isempty(strfind(str,'sample interval = '))
		    cfg.sample_int = str2num(char(regexp(str,'(?<=[\=]\>).\d+','match')));
           end
           if ~isempty(strfind(str,'samplenumber = '))
		    cfg.sample_num = str2num(char(regexp(str,'(?<=[\=]\>).\d+','match')));
	       end
	       str = fgetl(fid);
		end
    elseif strncmp(str,'**',2)
	  comments{i} = str;
    elseif isempty(str)
	  % empty line
    else
    end
    str = fgetl(fid);
    i=i+1;
end

cfg.start_time = char(regexp(fgetl(fid),'(?<=[\=]).*','match'));
junk = fgetl(fid); junk = fgetl(fid);

% read in the full data record
[mc] = textscan(fid,'%f %f %f %s %s','delimiter',',');

fclose(fid);




