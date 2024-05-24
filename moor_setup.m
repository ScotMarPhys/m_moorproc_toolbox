function moor_setup(varargin)
% moor_setup(MOORPROC_G, parameter, value)
% moor_setup(parameter, value)
% set paths and basic information for moored data processing,
% stored in global MOORPROC_G
%
% input parameter-value input pairs can include structure MOORPROC_G itself,
% and/or can include: 
%   datadir (top-level prefix like '/local/users/pstar/projects/rpdmoc/rapid'; if not set, query)
%   cruise (e.g. 'jc238'; otherwise, if MEXEC_G set, use that, else query)
%   YEAR (e.g. 2022; otherwise, if MEXEC_G set, use that, else query)
%   cruise_ctd (defaults to same as cruise)
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
addpath(genpath(fullfile(pathgit,'exec','gitrepo')))

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
        MOORPROC_G.YEAR = input('start year?   ');
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
if ~isfield(MOORPROC_G,'moordatadir') || ~isfield(MOORPROC_G,'reportdir') || ~isfield(MOORPROC_G,'ctddir')
    %define or request these directories
    if ~isfield(MOORPROC_G,'datadir') || isempty(MOORPROC_G.datadir)
        MOORPROC_G.datadir = input('base data directory (e.g. /data/pstar/projects/osnap, or D:\osnap) containing subdirectories data and documents ','s');
    end
    MOORPROC_G.moordatadir = fullfile(MOORPROC_G.datadir,'data','moor');
    MOORPROC_G.reportdir = fullfile(MOORPROC_G.datadir,'documents','datareports'); %for RAPID
    if ~exist(MOORPROC_G.reportdir,'dir')
        MOORPROC_G.reportdir = fullfile(MOORPROC_G.reportdir,'Documents','datareports'); %for OSNAP
    end
    if ~exist(MOORPROC_G.reportdir,'dir')
        MOORPROC_G.reportdir = input('directory to contain datareports (e.g. /data/pstar/projects/osnap/Documents/datareports) ','s');
    end
    if strcmp(MOORPROC_G.reportdir(end),'/') || strcmp(MOORPROC_G.reportdir(end),'\')
        MOORPROC_G.reportdir = MOORPROC_G.reportdir(1:end-1);
    end
    [~,d1] = fileparts(MOORPROC_G.reportdir);
    if ~strcmp(d1,MOORPROC_G.cruise)
        MOORPROC_G.reportdir = fullfile(MOORPROC_G.reportdir,MOORPROC_G.cruise);
    end
    MOORPROC_G.ctddir = fullfile(fileparts(MOORPROC_G.datadir),MOORPROC_G.cruise,'mcruise','data','ctd'); %RAPID
    if ~exist(MOORPROC_G.ctddir,'dir')
        MOORPROC_G.ctddir = fullfile(MOORPROC_G.datadir,'cruise_data',MOORPROC_G.cruise,'data','ctd'); %OSNAP
    end
    if ~exist(MOORPROC_G.ctddir,'dir')
        MOORPROC_G.ctddir = input('directory containing CTD data (e.g. /data/pstar/projects/osnap/cruise_data/jc238/mcruise/data/ctd ','s');
    end
    MOORPROC_G = rmfield(MOORPROC_G,'datadir');
end

if ~exist(MOORPROC_G.reportdir,'dir')
    mkdir(MOORPROC_G.reportdir);
end
MOORPROC_G.logdir  = fullfile(MOORPROC_G.moordatadir, 'logs');
if ~exist(MOORPROC_G.logdir,'dir')
    mkdir(MOORPROC_G.logdir);
end
if ~exist(MOORPROC_G.reportdir,'dir')
    mkdir(MOORPROC_G.reportdir);
end

%check directories at this stage
n = 0;
while n<4 && ~exist(MOORPROC_G.moordatadir,'dir')
    c = input(sprintf('directory for mooring data, %s, not found\n; create (1), change setting (2), or skip(3)?  ',MOORPROC_G.moordatadir),'n');
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
while n<4 && (~isfield(MOORPROC_G,'ctddir') || ~exist(MOORPROC_G.ctddir,'dir'))
    c = input(sprintf('directory for ctd data, %s, not found\n; create (1), change setting (2), or skip(3)?  ',MOORPROC_G.ctddatadir),'n');
    if c==1
        mkdir(MOORPROC_G.ctddir);
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
if ~exist(MOORPROC_G.ctddir,'dir')
    warning('timed out on setting valid MOORPROC_G.ctddir; expect mooring processing to fall over')
end

disp ([upper(MOORPROC_G.cruise) ', ' num2str(MOORPROC_G.YEAR)])
disp ('-----------------------------------------------------------------------------------------------')
fprintf(1,'ran the %s setup file to open the paths to\n%s\n%s\n%s\n',MOORPROC_G.cruise,pathgit,MOORPROC_G.moordatadir,MOORPROC_G.ctddir);

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




