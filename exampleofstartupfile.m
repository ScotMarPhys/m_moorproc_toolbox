%This is the startup file for processing hydro and/or RAPID and/or OSNAP (etc.) mooring data on NOC-OCP workstations
%Append it to, or call it from, your own MATLAB startup.m file, changing progdir and/or datadir (below) as necessary to reflect top-level paths on your computer
%note that after en705 much of the setting of processing paths has been broken out into moor_setup.m (although addpath is still done in this file), so that it is not necessary to have a startup{cruise}.m (e.g. startupjc238.m); computer/user-specific top-level paths are set in this file, and cases can be added to moor_setup.m for new project configurations. 

progdir = '/local/users/pstar/programs/';
datadir = '/local/users/pstar/projects/';
%example other paths for NOCS Linux machines, or a PC
%progdir = '/noc/mpoc/software/';
%datadir = '/noc/mpoc/';
%progdir = 'C:\Users/sa01ld\Desktop\OSNAP\';
%datadir = 'D:\osnap\';

% setup for hydro data processing
cdir = pwd;
[~,u] = system('whoami');
if strcmp(u,'pstar_devel') && contains(cdir,'test') && contains(cdir,'pstar_devel')
    addpath('/local/users/pstar_devel/programs/ocp/ocp_hydro_matlab')
    warning('you are running from test directory %s\n and using the pstar_devel copy/commit of ocp_hydro_matlab',cdir)
else
    addpath(fullfile(progdir,'ocp','ocp_hydro_matlab'))
end
m_setup; global MEXEC_G
mcruise = MEXEC_G.MSCRIPT_CRUISE_STRING;

% setup for mooring procssing
moorsw = 'moorproc';
switch moorsw

    case 'moorproc0'
        % setup for using m_moorproc_toolbox, SAMS version with separate
        % startup{cruise} file
        pathgit = fullfile(progdir,'m_moorproc_toolbox');
        addpath(pathgit)
        moorstart = ['startup' mcruise];
        if isempty(which(moorstart))
            warning('add m_moorproc_toolbox containing startup%d to path, enter to continue',MEXEC_G.MSCRIPT_CRUISE_STRING)
            pause
        end
        eval(moorstart)

    case 'moorproc'
        % setup for using m_moorproc_toolbox since en705
        pathgit = fullfile(progdir,'m_moorproc_toolbox');
        addpath(pathgit)
        if isempty(which('moor_setup'))
            warning('add m_moorproc_toolbox containing moor_setup to path, enter to continue',MEXEC_G.MSCRIPT_CRUISE_STRING)
            pause
        end
        moor_setup(datadir, progdir)

end

