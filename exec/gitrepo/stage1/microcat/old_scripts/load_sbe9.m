% load_sbe9.m
% load a converted SBE 9 data file and return
% Matlab serial date and data structure with whatever variables are there
% Adam Houk
% 4/17/2011
function [sb,cfg] = load_sbe9(filename);
fid = fopen(filename,'r');

% first line should contain the type of SeaBird used
headerline = fgetl(fid);

if any(findstr(headerline,'SBE 9'))
     cfg.instr_type = 'SBE 9';
elseif any(findstr(headerline,'SBE 19'))
     cfg.instr_type = 'SBE 19';
else
     disp('Unknown file-type!')
end

% system header
j=1;
str = fgetl(fid);
while strcmp('*',str(1))
	cfg.sys(j,:) = regexp(str,'\=','split');
	j=j+1;
	str = fgetl(fid);
end

% system parameters
i=1;j=1;
while (~strncmp(str,'*END*',5));
	if strncmp(str,'# nquan',7)
		cfg.ncol = sscanf(str,'%*9c %f');
	elseif strncmp(str,'# nvalues',9)
		cfg.nrow = sscanf(str,'%*11c %f');
	elseif strncmp(str,'# name',6)
		cfg.vars=sscanf(str(7:10),'%d',1);
		cfg.vars=cfg.vars+1;
		cfg.varlabel{cfg.vars}=str(findstr(str,'=')+2:min([findstr(str,':')-1 ...
			findstr(str,'/')-1 findstr(str,'-')-1]));
	elseif (strncmp(str,'# sensor',8))  
		sens=sscanf(str(10:11),'%d',1);
		cfg.sensors{sens+1}=str;
	elseif (strncmp(str,'# span',6))  
		sens=sscanf(str(8:9),'%d',1);
		cfg.spans{sens+1}=sscanf(str(findstr(str,'=')+1:end),'%f,%f');
	elseif (strncmp(str,'# bad_flag',10))  
		isub=13:length(str);
		bad_flag=sscanf(str(isub),'%g',1);
	elseif isempty(str)
		% empty line
     else
		% ignore line
	end
     str = fgetl(fid);
     i=i+1;
end

% read in the full data record (two ways to do it)
data = textscan(fid,'%f');
% or alternatively --> data = fscanf(fid,'%f',str2double([cfg.ncol,cfg.nrow]));
data2 = reshape(data{1},cfg.ncol,cfg.nrow)';
fclose(fid);

% clean up stuff
fclose all;
clear array1 sens isub str var fid i j

% assign header labels to data columns
for i=1:cfg.ncol
     eval(['sb.' cfg.varlabel{i} ' = data2(:,i);'])
end
