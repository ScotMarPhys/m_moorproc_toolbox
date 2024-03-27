function moor_setup(varargin)
% moor_setup(parameter, value)
% set paths and basic information for moored data processing,
% stored in global MOORPROC_G
%
% input parameter-value input pairs can include structure MOORPROC_G itself,
% and/or can include: 
%   datadir (top-level prefix like '/local/users/pstar/projects/'; if not set, query)
%   progdir (top-level prefix like '/local/users/pstar/programs/'; if not set, query)
% and may include
%   mcruise (e.g. 'jc238' if MEXEC_G set, use that, else query)
%   YEAR (e.g. 2022 if MEXEC_G set, use that, else query)
%   ctdcruise (defaults to same as mcruise)
%
%% set paths and information like deployment year for moored data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear MOORPROC_G
global MEXEC_G MOORPROC_G
for vno = 1:2:length(varargin)
    eval([varargin{vno} ' = varargin{vno+1};']);
end

if ~isstruct(MOORPROC_G) || ~isfield(MOORPROC_G,'cruise')
    if exist('cruise','var')
        MOORPROC_G.cruise = cruise;
    elseif isstruct(MEXEC_G)
        MOORPROC_G.cruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
    else
        MOORPROC_G.cruise = input('mooring cruise name?   ','s');
    end
end
if ~exist('YEAR','var')
    if exist('MEXEC_G','var') && isstruct(MEXEC_G)
        YEAR = MEXEC_G.MDEFAULT_DATA_TIME_ORIGIN(1);
    else
        YEAR = input('start year?   ','n');
    end
end
YEAR
if ~isfield(MOORPROC_G,'cruise_ctd')

    if ~exist('mcruise','var')
        if isstruct(MEXEC_G)
            mcruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
        else
            mcruise = input('cruise name?   ','s');
        end
    end
    MOORPROC_G.cruise = mcruise;
    if exist('ctdcruise','var')
        MOORPROC_G.cruise_ctd = ctdcruise;
    else
        %default: same as mcruise
        MOORPROC_G.cruise_ctd = mcruise;
    end
    MOORPROC_G.YEAR = YEAR;
end

%if exist('operator', %MOORPROC_G.operator = input('operator?   ','s');
MOORPROC_G.project = input('project?   ','s');
%MOORPROC_G.project = 'RAPID';
MOORPROC_G.operator = 'ylf';

% Add to the path osnap mooring functions to go with this file:
pathgit = fileparts(which(mfilename));
addpath(genpath(fullfile(pathgit,'exec','gitrepo')))

% Define where to find the mooring and other data
% pathdata   = 'C:\Users\sa01ld\OneDrive - SAMS\OSNAP_mooring_processing\osnap';
% pathdata   = 'D:\osnap';
if ~exist('datadir','var') || ~exist(datadir,'dir') %should be set in startup for a given user/computer
    datadir = '~';
end
switch lower(MOORPROC_G.project)
    case 'rapid'
        basedir = fullfile(datadir, 'rpdmoc');
        MOORPROC_G.moordatadir = fullfile(basedir,'rapid','data','moor');
        MOORPROC_G.reportdir = fullfile(basedir,'rapid','documents','datareports',MOORPROC_G.cruise);
        MOORPROC_G.ctddir = fullfile(basedir,MOORPROC_G.cruise,'data','ctd');
    case 'osnap'
        basedir = fullfile(datadir, 'osnap');
        MOORPROC_G.moordatadir = fullfile(basedir,'data','moor');
        MOORPROC_G.reportdir = fullfile(basedir,'Documents','datareports',MOORPROC_G.cruise);
        MOORPROC_G.ctddir = fullfile(basedir,'cruise_data',MOORPROC_G.cruise,'data','ctd');
    otherwise
        basedir = fullfile(datadir, lower(MOORPROC_G.project));
        MOORPROC_G.moordatadir = fullfile(basedir,'data','moor');
        MOORPROC_G.reportdir = fullfile(basedir,'documents','datareports',MOORPROC_G.cruise);
        MOORPROC_G.ctddir = fullfile(basedir,MOORPROC_G.cruise,'data','ctd');
end
MOORPROC_G.logdir  = fullfile(MOORPROC_G.moordatadir, 'logs');
if ~exist(MOORPROC_G.logdir,'dir')
    mkdir(MOORPROC_G.logdir)
end
if ~exist('MOORPROC_G.reportdir','dir')
    mkdir(MOORPROC_G.reportdir)
end

disp ([upper(mcruise) ', ' num2str(YEAR)])
disp ('-----------------------------------------------------------------------------------------------')
fprintf(1,'%s\n',['ran the ' mcruise ' setup file to open the paths to'], ...
    pathgit,MOORPROC_G.moordatadir,MOORPROC_G.ctddir);

% display version of m_moorproc_toolbox
d = pwd; cd(pathgit)
[s,c] = system('git branch | grep "*"');
if s==0 && ~contains(c, 'fatal:'); fprintf(1,'using m_moorproc_toolbox branch %s',c); end
cd(d)

disp('MOORPROC_G contains:')
disp(MOORPROC_G)
%check directories at this stage
n = 0; while n<4 && ~exist(MOORPROC_G.moordatadir,'dir')
    c = input('directory for ctd data, %s, not found\n; create (1), change setting (2), or skip(3)?  ','n');
    if c==1
        mkdir(MOORPROC_G.moordatadir);
    elseif c==2
        MOORPROC_G.moordatadir = input('directory for mooring data:    ','s');
    elseif c==3
        n = 6;
    end
    n = n+1;
end
n = 0; while n<4 && ~exist(MOORPROC_G.ctddir,'dir')
    if c==1
        mkdir(MOORPROC_G.ctddir);
    elseif c==2
        MOORPROC_G.moordatadir = input('directory for mooring data:    ','s');
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



