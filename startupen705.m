%% startup file for cruise data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define GLOBAL variables
global basedir execdir pathgit pathdata this_cruise cruise_ctd YEAR

% Add to the path osnap mooring functions:
% pathdata   = 'C:\Users\sa01ld\OneDrive - SAMS\OSNAP_mooring_processing\osnap';
% pathdata   = 'D:\osnap';
pathdata = '/local/users/pstar/projects/rpdmoc/'; if ~exist(pathdata,'dir'); pathdata = pwd; end
%pathgit     = 'C:\Users\sa01ld\Desktop\OSNAP\m_moorproc_toolbox';   
% pathgit     = 'C:\Users\SA01LD\m_moorproc_toolbox';   
pathgit     = '/local/users/pstar/programs/m_moorproc_toolbox'; 

disp ('EN705, 2023') 
disp ('-----------------------------------------------------------------------------------------------')
disp (['this is the en705 startup file to open the RAPID paths to '...
    pathdata '  ' pathgit])

if exist('MEXEC_G','var') && 1 % change this flag if you have mexec but it's set for a different cruise
    cruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
    this_cruise = cruise;
    cruise_ctd = cruise;
    YEAR = MEXEC_G.MDEFAULT_DATA_TIME_ORIGIN(1);
else %set manually here
    cruise = 'jc238';
    this_cruise = cruise;
    cruise_ctd = 'JC238';
    YEAR = 2022;
end

moordatadir     = [ pathdata filesep 'rapid' filesep 'data' filesep 'moor'];
ctddir = [pathdata filesep cruise_ctd filesep 'mcruise' filesep 'data' filesep 'ctd']
execdir     = [ pathgit filesep 'exec/']; 

baselogdir  = [basedir];
addpath(genpath([execdir]))

% CD to working directory
cd([pathgit '/exec/gitrepo']) 
                            
p=path;

