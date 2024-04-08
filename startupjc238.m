%% startup file for cruise data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define GLOBAL variables
global basedir datadir execdir pathgit pathosnap

% Add to the path osnap mooring functions:
[s,r] = system('whoami');
if s==0 && strncmp(r,'yvonng',6)
%    pathosnap = '/Users/yvonng/projects/osnap';
    pathosnap = '/Users/yvonng/programs/m_moorproc_toolbox';
    pathgit = '/Users/yvonng/programs/m_moorproc_toolbox';
else
    pathosnap   = 'C:\Users\sa01ld\Desktop\m_moorproc_toolbox';
    pathgit     = 'C:\Users\sa01ld\Desktop\m_moorproc_toolbox';   
end


disp ('JC238, 2022')
disp ('-----------------------------------------------------------------------------------------------')
disp (['this is the JC238 startup file to open the OSNAP paths to '...
    pathosnap '  ' pathgit])

basedir     = pathosnap; 
datadir     = fullfile(pathgit, 'data');
execdir     = fullfile(pathgit, 'exec');

this_cruise = 'jc238';
baselogdir  = basedir;
addpath(genpath(execdir))

% CD to working directory
cd(fullfile(execdir,'gitrepo'))
                            
p=path;
