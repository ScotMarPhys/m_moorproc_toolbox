% Add to the path osnap mooring functions:
pathosnap   = '/Users/locupe/Dropbox/Work/osnap';
pathgit = '/Users/locupe/Dropbox/Work/Python/Repos_perso/m_moorproc_toolbox';   

%% startup file for cruise data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp ('Welcome back to Matlab!')
disp ('------------------------------------------------------------------------')
disp('    ')
disp ('DY120, 2020')
disp('    ')
disp ('------------------------------------------------------------------------')
disp (['  OSNAP data can be accessed by runing the command: '])
disp ([               '"cd(pathosnap)" in matlab;'])
disp (['  The mooring processing scripts can be accessed by '])
disp ([               ' runing the command : "cd(pathgit)"'])
disp('    ')
disp ('------------------------------------------------------------------------')
basedir = [ pathosnap filesep]; 
execdir = [ pathgit filesep 'exec' filesep]; 

this_cruise = 'dy120';
baselogdir = [basedir];

% Add git repository to matlab path 
disp('    ')
disp([pathgit filesep ' and'])
disp(['subfolders added to matlab path'])
disp ('------------------------------------------------------------------------')
disp('    ')
addpath(genpath([execdir]))

p=path;


