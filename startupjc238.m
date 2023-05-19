%% startup file for cruise data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define GLOBAL variables
global basedir datadir execdir pathgit pathosnap

% Add to the path osnap mooring functions:
pathosnap   = 'C:\Users\sa01ld\m_moorproc_toolbox';
pathgit     = 'C:\Users\sa01ld\m_moorproc_toolbox';   

disp ('JC238, 2022')
disp ('-----------------------------------------------------------------------------------------------')
disp (['this is the JC238 startup file to open the OSNAP paths to '...
    pathosnap '  ' pathgit])

basedir     = [ pathosnap filesep]; 
datadir     = [ pathgit filesep 'data\' ]; 
execdir     = [ pathgit filesep 'exec\']; 

this_cruise = 'jc238';
baselogdir  = [basedir];
addpath(genpath([execdir]))

% CD to working directory
cd([pathgit '\exec\gitrepo']) 
                            
p=path;
