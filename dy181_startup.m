%This is the Matlab startup file for processing hydro and/or RAPID and/or OSNAP (etc.) mooring data at sea

%%% change the next lines to reflect your directory structure %%%
project = 'OSNAP';
basedir = '/Users/tiadot/Library/CloudStorage/OneDrive-NOC/WORK/DATA'; %contains osnap, or rpdmoc
progdir = '/Users/tiadot/Library/CloudStorage/OneDrive-NOC/WORK/PROGRAMS'; %contains m_moorproc_toolbox
use_mexec = 0;
cruise = 'dy181';
if use_mexec
    MEXEC_G_user.other_programs_root = '/data/pstar/programs/others/'; %gsw, etc.
    MEXEC_G_user.mexec_data_root = '/Volumes/mpoc/osnap/dy181/backup_20240729112645/data'; %mexec hydro data
end
%example other paths for NOCS Linux machines
%basedir = '/noc/mpoc/';
%progdir = '/noc/mpoc/drake/programs/';
%or for a PC
%basedir = 'D:\';
%progdir = 'D:\programs\';

% setup for hydro data processing including running m_setup
addpath(genpath(fullfile(progdir,'ocp_hydro_matlab')))
global MOORPROC_G
if use_mexec
    path_choose = m_setup(MEXEC_G_user); %m_setup returns 1 if cruise options/user selects to process LADCP rather than moored data
    m_common %global variables to workspace
else
    path_choose = 2;
    MOORPROC_G.cruise = 'dy181';
    MOORPROC_G.cruise_ctd = 'dy181';
end
MOORPROC_G.YEAR = 2024;
MOORPROC_G.ctddir = '/Volumes/mpoc/osnap/dy181/backup_20240729112645/data/ctd';

if path_choose==0 || path_choose==2
    % setup for mooring procssing
    % setup for using m_moorproc_toolbox since en705
    pathgit = fullfile(progdir,'m_moorproc_toolbox');
    addpath(genpath(pathgit))
    if isempty(which('moor_setup'))
        warning('add m_moorproc_toolbox containing moor_setup to path, enter to continue',MEXEC_G.MSCRIPT_CRUISE_STRING)
        pause
    end
    MOORPROC_G.project = project;
    switch project
        case 'RAPID'
            MOORPROC_G.datadir = fullfile(basedir,'rpdmoc/rapid');
        case 'OSNAP'
            MOORPROC_G.datadir = fullfile(basedir,'osnap');
    end
    moor_setup(MOORPROC_G)
    expa = which('mc_call_caldip');
    global MOORPROC_G
    if contains(expa,MOORPROC_G.cruise) && ~contains(expa,'gitrepo')
        warning('is this where you intend your mooring processing tools to come from? (%s)',fileparts(fileparts(expa)))
    end

end
