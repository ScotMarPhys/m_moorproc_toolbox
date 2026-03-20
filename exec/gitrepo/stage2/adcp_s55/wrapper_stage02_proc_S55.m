% WRAPPER_STAGE02_PROC_S55: Processes S55 ADCP mooring data.
%
% This script manages environment setup and calls stage02_proc_S55.
%
% CONFIGURATION MODES:
%   1. AUTOMATIC: Uses global MOORPROC_G (via m_moorprog_startupfiles).
%   2. MANUAL: Prompts for Cruise ID and uses PC-specific local paths.
%
% MAINTENANCE:
%   - NEW CRUISES: Add a 'case' in Section 1 to map Cruise ID to Mooring Name.
%   - NEW USERS: Add a 'case' in Section 2 to define paths for your 'COMPUTERNAME'.
%
% DEPENDENCIES:
%   - m_moorproc_toolbox (for rodbload.m)
%   - GSW Oceanographic Toolbox (v3.06+)
%
% Authors: K Burmeister, S Jones | Date: 02/2026

close all; clear;
global MOORPROC_G

%% 1. Identify Cruise and Map to Mooring Name
if ~isempty(MOORPROC_G) && isfield(MOORPROC_G, 'cruise')
    cruise = lower(MOORPROC_G.cruise);
else
    % Prompt user for input if global setup is missing
    cruise = lower(input('Enter Cruise ID (e.g., jc238, dy181, dy214): ', 's'));
end

% Automatic Mooring Name Mapping: [moorname_deployment#_year]
switch cruise
    case 'jc238', moor = 'rhadcp_01_2020';
    case 'dy181', moor = 'rhadcp_02_2022';
    case 'dy214', moor = 'rhadcp_03_2024';
    otherwise,    error('Cruise ID "%s" not recognized in mapping.', cruise);
end

%% 2. Path Configuration by Computer
pc_name = getenv('COMPUTERNAME');

switch pc_name
    case 'SA07KB-3JN9YY2'
        basedir      = 'C:\Users\sa07kb\Projects\Moor_Data_Proc\';
        pathgit      = fullfile(basedir, 'm_moorproc_toolbox');
        dataindir    = fullfile(basedir, 'moor_examples\osnap\data\moor\proc\', moor, 'adcp_S55');
        infofile     = fullfile(basedir, 'moor_examples\osnap\data\moor\proc\', moor, [moor 'info.dat']);
        outdir       = dataindir;
        toolbox_path = 'C:\Users\sa07kb\Matlab\toolboxes\gsw_matlab_v3_06_16';
        filename     = '200044_data';
        scotiafile = fullfile(basedir, 'moor_examples\osnap\scotia/SCOTIA_monthly_clim_V8.nc');
    case 'SA01SJ-G9WC2J3'
        basedir      = 'D:\Work_computer_sync\OSNAP_postdoc\Python\';
        pathgit      = fullfile(basedir, 'm_moorproc_toolbox');   
        dataindir    = 'D:\Work_computer_sync\OSNAP_postdoc\Mooring\RHADCP\plots\'; 
        infofile     = fullfile('E:\OSNAP\RHADCP\DY181\S200044A012_RHAD2_JC238\conversion2\', [moor 'info.dat']);
        outdir       = 'D:\Work_computer_sync\OSNAP_postdoc\Mooring\RHADCP\plots\';
        toolbox_path = 'D:\Work_computer_sync\MATLAB_functions';
        filename     = 'S200044A012_RHAD2_JC238';
        scotiafile = '';
    otherwise
        error('PC name "%s" not recognized. Add your paths to Section 2.', pc_name);
end

%% 3. Environment Setup and Execution
addpath(genpath(fullfile(pathgit, 'exec', 'gitrepo')));
addpath(genpath(toolbox_path));
logfile = fullfile(outdir, [moor '_ADCP_stage2.log']);

if ~isempty(MOORPROC_G)
    % Standard toolbox execution
    stage02_proc_S55(moor);
else
    % Manual execution with local variables
    stage02_proc_S55(moor, dataindir, infofile, logfile, outdir,scotiafile);
end
