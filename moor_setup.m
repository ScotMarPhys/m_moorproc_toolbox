function moor_setup(varargin)
% moor_setup(MOORPROC_G, parameter, value)
% moor_setup(parameter, value)
% set paths and basic information for moored data processing,
% stored in global MOORPROC_G
%
% input parameter-value input pairs can include structure MOORPROC_G itself,
% and/or can include: 
%   basedir (top-level prefix like '/local/users/pstar/projects/rpdmoc'; if not set, query)
%   cruise (e.g. 'jc238'; otherwise, if MEXEC_G set, use that, else query)
%   YEAR (e.g. 2022; otherwise, if MEXEC_G set, use that, else query)
%   cruise_ctd (defaults to same as cruise)
%
% will try to set paths as usual for RAPID and OSNAP projects, otherwise
% will query the user where to find ctd data, moored data, and where to put
% summaries/plots
% to avoid the default paths for RAPID and OSNAP, do:
% MOORPROC_G.basedir = '';
% moor_setup(MOORPROC_G)
%
%% set paths and information like deployment year for moored data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global MEXEC_G MOORPROC_G
if isempty(MOORPROC_G)
    clear MOORPROC_G
end

%parameters input to function
if nargin>0
    if isstruct(varargin{1})
        MOORPROC_G = varargin{1};
        varargin(1) = [];
    end
    for vno = 1:2:length(varargin)
        MOORPROC_G.(varargin{vno}) = varargin{vno+1};
    end
end

% Add to the path mooring functions to go with this file:
pathgit = fileparts(which(mfilename));
addpath(genpath(fullfile(pathgit, 'code')))

%cruise, cruise_ctd, YEAR, project, operator
if isstruct(MEXEC_G)
    if ~isfield(MOORPROC_G,'cruise')
        MOORPROC_G.cruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
    end
    if ~isfield(MOORPROC_G,'cruise_ctd')
        MOORPROC_G.cruise_ctd = MEXEC_G.MSCRIPT_CRUISE_STRING;
    end
else
    if ~isfield(MOORPROC_G,'cruise')
        MOORPROC_G.cruise = input('mooring cruise name?   ','s');
    end
    if ~isfield(MOORPROC_G,'cruise_ctd')
        MOORPROC_G.cruise_ctd = input('CTD cruise name? ','s');
    end
end

if ~isfield(MOORPROC_G,'YEAR')
    if isstruct(MEXEC_G)
        MOORPROC_G.YEAR = MEXEC_G.MDEFAULT_DATA_TIME_ORIGIN(1);
    else
        %MOORPROC_G.YEAR = input('start year?   ');
        MOORPROC_G.YEAR = input('start year of this cruise?   ');
    end
end

if ~isfield(MOORPROC_G,'project')
    MOORPROC_G.project = input('project?   ','s');
end
if ~isfield(MOORPROC_G,'operator')
    [s,u] = system('whoami');
    if s==0
        MOORPROC_G.operator = replace(u,newline,'');
    else
        MOORPROC_G.operator = input('operator?    ','s');
    end
end

% Define where to find the mooring and other data (formerly pathdata)
dirs_req = {'ctddatadir', 'reportdir', 'moordatadir'};
fn = fieldnames(MOORPROC_G);
if isfield(MOORPROC_G,'basedir') && isempty(MOORPROC_G.basedir)
    MOORPROC_G = rmfield(MOORPROC_G,'basedir');
elseif sum(ismember(fn,dirs_req))<length(dirs_req)
    MOORPROC_G = project_default_dirs(MOORPROC_G);
end
fn = fieldnames(MOORPROC_G);
if sum(ismember(fn,dirs_req))<length(dirs_req)
    if ~isfield(MOORPROC_G,'ctddatadir')
        MOORPROC_G.ctddatadir = input('directory containing 1 Hz CTD data   ','s');
    end
    if ~isfield(MOORPROC_G,'moordatadir')
        MOORPROC_G.moordatadir = input('directory containing moored data (should have subdirectories proc, proc_calib, raw)   ','s');
    end
    if ~isfield(MOORPROC_G,'reportdir')
        MOORPROC_G.reportdir = input('directory for data reporting (statistics and figures)   ','s');
    end
end
if ~exist(MOORPROC_G.reportdir,'dir')
    c = input(sprintf('data report directory %s does not exist; create it? (y/n)  ',strrep(MOORPROC_G.reportdir, filesep, '/')),'s');
    if strcmp(c,'y')
        mkdir(MOORPROC_G.reportdir)
    end
end

%check directories at this stage
n = 0;
while n<4 && ~exist(MOORPROC_G.moordatadir,'dir')
    c = input(sprintf('directory for mooring data, %s, not found\n; create (1), change setting (2), or skip(3)?  ',strrep(MOORPROC_G.moordatadir, filesep, '/')),'s');
    if c==1
        mkdir(MOORPROC_G.moordatadir);
    elseif c==2
        MOORPROC_G.moordatadir = input('directory for mooring data:    ','s');
    elseif c==3
        n = 6;
    end
    n = n+1;
end
n = 0;
while n<4 && (~isfield(MOORPROC_G,'ctddatadir') || ~exist(MOORPROC_G.ctddatadir,'dir'))
    c = input(sprintf('directory for ctd data, %s, not found\n; create (1), change setting (2), or skip(3)?  ',strrep(MOORPROC_G.ctddatadir, filesep, '/')));
    if c==1
        mkdir(MOORPROC_G.ctddatadir);
    elseif c==2
        MOORPROC_G.ctddatadir = input('directory for mooring data:    ','s');
    elseif c==3
        n = 6;
    end
    n = n+1;
end
if ~exist(MOORPROC_G.moordatadir,'dir')
    warning('timed out on setting valid MOORPROC_G.moordatadir; expect mooring processing to fall over')
end
if ~exist(MOORPROC_G.ctddatadir,'dir')
    warning('timed out on setting valid MOORPROC_G.ctddatadir; expect mooring processing to fall over')
end

disp ([upper(MOORPROC_G.cruise) ', ' num2str(MOORPROC_G.YEAR)])
disp ('-----------------------------------------------------------------------------------------------')
fprintf(1,'ran the %s setup file to open the paths to\n%s\n%s\n%s\n',MOORPROC_G.cruise,pathgit,MOORPROC_G.moordatadir,MOORPROC_G.ctddatadir);

% display version of m_moorproc_toolbox
d = pwd; cd(pathgit)
[s,c] = system('git branch | grep "*"');
if s==0 && ~contains(c, 'fatal:')
    c = replace(c,newline,'');
    fprintf(1,'using m_moorproc_toolbox branch %s\n',c); 
end
cd(d)

disp('MOORPROC_G contains:')
disp(MOORPROC_G)


function MOORPROC_G = project_default_dirs(MOORPROC_G)

if ~isfield(MOORPROC_G,'project') || ~isfield(MOORPROC_G,'location')
    return
end

pl = [MOORPROC_G.project '_' MOORPROC_G.location];
switch pl
    case 'RAPID_atnoc'
        if ~isfield(MOORPROC_G,'basedir')
            MOORPROC_G.basedir = input('base data directory (e.g. /data/pstar/projects/rpdmoc, containing subdirectories cruise_data and rapid) ','s');
        end
        MOORPROC_G.cruisedir = fullfile(MOORPROC_G.basedir,'cruise_data',MOORPROC_G.cruise);
        if ~exist(MOORPROC_G.cruisedir,'dir')
            MOORPROC_G.cruisedir = fullfile(MOORPROC_G.basedir,MOORPROC_G.cruise);
        end
        MOORPROC_G.ctddatadir = fullfile(MOORPROC_G.cruisedir,'mcruise','data','ctd');
        MOORPROC_G.reportdir = fullfile(MOORPROC_G.cruisedir,'report_tables');
        MOORPROC_G.moordatadir = fullfile(MOORPROC_G.basedir,'rapid','data','moor');
    case 'RAPID_atsea'
        if ~isfield(MOORPROC_G,'basedir')
            MOORPROC_G.basedir = input('base data directory (e.g. /data/pstar/projects/rpdmoc, containing subdirectories cruises and moorings) ','s');
        end
        MOORPROC_G.cruisedir = fullfile(MOORPROC_G.basedir,'cruises',MOORPROC_G.cruise);
        MOORPROC_G.ctddatadir = fullfile(MOORPROC_G.cruisedir,'data','ctd');
        MOORPROC_G.reportdir = fullfile(MOORPROC_G.cruisedir,'report_tables');
        MOORPROC_G.moordatadir = fullfile(MOORPROC_G.basedir,'moorings');
    case 'OSNAP_atsea'
        if ~isfield(MOORPROC_G,'basedir')
            MOORPROC_G.basedir = input('base data directory (e.g. /data/pstar/projects/osnap, containing subdirectories cruises and moorings) ','s');
        end
        MOORPROC_G.cruisedir = fullfile(MOORPROC_G.basedir,'cruises',MOORPROC_G.cruise);
        MOORPROC_G.ctddatadir = fullfile(MOORPROC_G.cruisedir,'data','ctd');
        MOORPROC_G.reportdir = fullfile(MOORPROC_G.cruisedir,'report_tables');
        MOORPROC_G.moordatadir = fullfile(MOORPROC_G.basedir,'moorings');
end
