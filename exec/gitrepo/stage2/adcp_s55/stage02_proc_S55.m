% stage02_proc_S55(moor,varargin)
% basic preprocessing for stage 1 Signature 55 ADCP data (nc-format)
% Please make sure raw ADCP data in .mat format as exported by Signature
% Viewer were QC by calling stage01_read_qc_S55 before calling this
% function
%
% required inputs: moor - mooring name e.g. 'rhadcp_01_2020'
%                  
% optional inputs:     dataindir = varargin{1};
%                      infofile = varargin{2};
%                      logfile = varargin{3};
%                      outdir = varargin{4}; #output dir (fig, data, logs)
% functions called:    rodbload, gsw_z_from_p
% features
% 

function stage02_proc_S55(moor, varargin)

if nargin==0
    help stage02_proc_S55
    return
end

if nargin==1
    global MOORPROC_G
    operator = MOORPROC_G.operator;
    pd = moor_inoutpaths('adcp_S55',moor);
    dataindir = pd.stage1path;
    infofile = pd.infofile;
    logfile = pd.stage1log;
    outdir = pd.stage1path;
    infile_form = pd.stage1form;
    ouput_form = pd.stage2form;
else
    operator = getenv('COMPUTERNAME');    
    dataindir = varargin{1};
    infofile = varargin{2};
    logfile = varargin{3};
    outdir = varargin{4};
    infile_form =  [moor '_%d_stage1.nc'];
    ouput_form = [moor '_%d_stage2.nc'];
end

if ~exist(outdir,'dir')
    mkdir(outdir)
end

% read in meta data
info_adcp = read_adcp_infofile(infofile);

% start log entry
fidlog   = fopen(logfile,'w');
fprintf(fidlog, '\n==== START ENTRY  =====\n');
fprintf(fidlog,['Processing steps of Signature 55 ADCP data taken by ' ...
    mfilename ':\n']);
% fprintf(fidlog,'  1. eliminate lauching and recovery period\n'); %edit
% fprintf(fidlog,'  2. save data to rodb file\n');
fprintf(fidlog,'\n Operated by:%s on %s\n',operator,datestr(clock)); 
fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',info_adcp.lat);
fprintf(fidlog,'Longitude %6.3f \n\n\n',info_adcp.lon);

info_adcp.start_timestamp = datenum(datetime([info_adcp.s_d.' , info_adcp.s_t.' , 0])); %start date and time gb
info_adcp.end_timestamp = datenum(datetime([info_adcp.e_d.' , info_adcp.e_t.' , 0])); %end date and time ed
test = mfilename
if isnan(info_adcp.id)
   fprintf('No serial number given, use default 200044');
   serial_nums=200044;
else
    vec=find((info_adcp.id>=319) & (info_adcp.id <=328)); % Possible ADCP codes - taken from IMP moorings package
    serial_nums=info_adcp.sn(vec);
end

%% Load data: Config, Data, Description
for i = 1:length(serial_nums)
    fprintf('Processing sn %d',serial_nums(i));

outfile = fullfile(outdir,sprintf(ouput_form,serial_nums(i)));
infile = fullfile(dataindir,sprintf(infile_form,serial_nums(i)));
DS = s55_load_nc_as_struct(infile);

fprintf(fidlog,'infile : %s\n',infile);
fprintf(fidlog,'ADCP serial number  : %d\n',serial_nums(i));

end