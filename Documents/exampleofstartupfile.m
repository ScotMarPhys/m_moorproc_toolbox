% Add to the path osnap mooring functions:
pathosnap   = '/Users/locupe/Dropbox/Work/osnap';
pathgit = '/Users/locupe/Dropbox/Work/Python/Repos_perso/m_moorproc_toolbox';   

%% startup file for cruise data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp ('Welcome back to Matlab!')
disp ('-----------------------------------------------------------------------------------------------')
disp('    ')
disp ('DY120, 2020')
disp ('-----------------------------------------------------------------------------------------------')
disp (['this is the dy120 startup file to open the OSNAP paths to '...
    'osnap/data/moor/ and osnap/exec/dy120, '])

basedir = [ pathosnap filesep]; 
datadir = [ pathgit filesep 'data' filesep ]; 
execdir = [ pathgit filesep 'exec' filesep]; 

this_cruise = 'dy120';
baselogdir = [basedir];
addpath(genpath([execdir]))
cd([pathgit filesep 'exec' filesep 'gitrepo']) 

p=path;
