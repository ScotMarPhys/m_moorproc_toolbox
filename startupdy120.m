%% startup file for cruise data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define GLOBAL variables
global basedir datadir execdir pathgit pathosnap

% Add to the path osnap mooring functions:
pathosnap   = 'C:\Users\sa01ld\Desktop\m_moorproc_toolbox';
pathgit     = 'C:\Users\sa01ld\Desktop\m_moorproc_toolbox';   

disp ('DY120, 2020')
disp ('-----------------------------------------------------------------------------------------------')
disp (['this is the dy120 startup file to open the OSNAP paths to '...
    pathosnap '  ' pathgit])

basedir     = [ pathosnap filesep]; 
datadir     = [ pathgit filesep 'data\' ]; 
execdir     = [ pathgit filesep 'exec\']; 

this_cruise = 'dy120';
baselogdir  = [basedir];
addpath(genpath([execdir]))

% CD to working directory
cd([pathgit '\exec\gitrepo']) 
                            
p=path;
