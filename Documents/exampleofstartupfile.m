%This is the Matlab startup file for processing hydro and/or RAPID and/or OSNAP (etc.) mooring data at sea

%%% change the next lines to reflect your directory structure %%%
project = 'RAPID';
basedir = '/data/pstar/projects/'; %contains osnap, or rpdmoc
progdir = '/data/pstar/programs/gitvcd/'; %contains m_moorproc_toolbox
use_mexec = 0;
if use_mexec
    MEXEC_G_user.other_programs_root = '/data/pstar/programs/others/'; %gsw, etc.
    MEXEC_G_user.mexec_data_root = '/data/pstar/cruise/data'; %mexec hydro data
end
%example other paths for NOCS Linux machines
%basedir = '/noc/mpoc/';
%progdir = '/noc/mpoc/drake/programs/';
%or for a PC
%basedir = 'D:\';
%progdir = 'D:\programs\';

% setup for hydro data processing including running m_setup
addpath(fullfile(progdir,'ocp_hydro_matlab'))
if use_mexec
    path_choose = m_setup(MEXEC_G_user); %m_setup returns 1 if cruise options/user selects to process LADCP rather than moored data
    m_common %global variables to workspace
else
    path_choose = 2;
end

if path_choose==0 || path_choose==2
    % setup for mooring procssing
    % setup for using m_moorproc_toolbox since en705
    pathgit = fullfile(progdir,'m_moorproc_toolbox');
    addpath(pathgit)
    if isempty(which('moor_setup'))
        warning('add m_moorproc_toolbox containing moor_setup to path, enter to continue',MEXEC_G.MSCRIPT_CRUISE_STRING)
        pause
    end
    switch project
        case 'RAPID'    
	    moor_setup('datadir', fullfile(basedir,'rpdmoc/rapid'), 'project', project);
        case 'OSNAP'
            moor_setup('datadir', fullfile(basedir,'osnap'), 'project', project);
        otherwise
            moor_setup('project',project) %prompt for datadir
    end
    expa = which('mc_call_caldip');
    if contains(expa,mcruise) && ~contains(expa,'gitrepo')
        warning('is this where you intend your mooring processing tools to come from? (%s)',fileparts(fileparts(expa)))
    end

end
