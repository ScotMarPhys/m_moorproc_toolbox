%% startup file for cruise data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define GLOBAL variables
global MEXEC_G
global MOORPROC_G

MOORPROC_G.operator = 'ylf';

if isstruct(MEXEC_G) && 1 %can change this flag if you have m_setup for a different cruise
    mcruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
    YEAR = MEXEC_G.MDEFAULT_DATA_TIME_ORIGIN(1);
else
    if ~exist('mcruise','var')
        mcruise = input('cruise name?','s');
        YEAR = input('start year?','n');
    end
end
MOORPROC_G.cruise = mcruise;
MOORPROC_G.cruise_ctd = mcruise;
MOORPROC_G.YEAR = YEAR;


% Add to the path osnap mooring functions:
%pathgit     = 'C:\Users\sa01ld\Desktop\OSNAP\m_moorproc_toolbox';   
% pathgit     = 'C:\Users\SA01LD\m_moorproc_toolbox';   
pathgit     = '/local/users/pstar/programs/m_moorproc_toolbox'; 
d = pwd; cd(pathgit)
[s,c] = system('git branch | grep "*"');
if s==0 && ~contains(c, 'fatal:'); fprintf(1,'using m_moorproc_toolbox branch %s',c); end
cd(d)

% Define where to find the mooring and other data
% pathdata   = 'C:\Users\sa01ld\OneDrive - SAMS\OSNAP_mooring_processing\osnap';
% pathdata   = 'D:\osnap';
basedir = '/local/users/pstar/projects/';
pathdata = fullfile(basedir,'rpdmoc','rapid','data','moor');
%pathdata = '/local/users/pstar/projects/osnap/data/moor';
if ~exist(pathdata,'dir'); pathdata = pwd; end
pathctd = fullfile(basedir,'rpdmoc',mcruise,'mcruise','data','ctd');

disp ([upper(mcruise) ', ' num2str(YEAR)]) 
disp ('-----------------------------------------------------------------------------------------------')
fprintf(1,'%s\n',['this is the ' mcruise ' startup file to open the paths to'], ...
    pathgit,pathdata,pathctd);

MOORPROC_G.moordatadir = pathdata;
MOORPROC_G.ctddir = pathctd;
MOORPROC_G.execdir = fullfile(pathgit,'exec');
MOORPROC_G.logdir = basedir;

baselogdir  = fullfile(pathdata,'logs');
if ~exist(baselogdir,'dir')
    mkdir(baselogdir)
end
addpath(genpath(fullfile(pathgit,'exec','gitrepo')))

MOORPROC_G
disp('if necessary, MOORPROC_G fields (e.g. operator) can be changed in workspace for this session')
