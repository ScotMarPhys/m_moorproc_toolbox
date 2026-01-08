% Call function to read and quality controle S55 ADCP data
%
% K Burmeister, S Jones 12/2025

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PARAMETER PRAEMBLE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% Need to be updated for:
% - directories and file names
% - add additional deployment periods (add moorX = 'filename')
% - turn off/on checkplots
% - edit start/end of total time period
% - data version
% - set depth of shallowest instrument (idepth)
% - preamble for despiking
% - if you add a deployment period, you need to edit steps 1-3
%---------

close all
clear

global MOORPROC_G

moor = 'rhadcp_01_2020';

% TO UPDATE <------------------------
% in- and output directories
pc_name = getenv('COMPUTERNAME');
if strcmp(pc_name,'SA07KB-3JN9YY2');
    basedir = 'C:\Users\sa07kb\Projects\Moor_Data_Proc\';
    dataindir = [basedir,'moor_examples\osnap\data\moor\raw\jc238\adcp_s55\'];
    pathgit = [basedir 'm_moorproc_toolbox\'];
    figureoutdir = fullfile(basedir,'moor_examples\osnap\data\moor\proc', ...
        moor,'adcp_S55');
    logfile = fullfile(figureoutdir,[moor '_ADCP_stage1.log']);
    infofile = fullfile(basedir,'moor_examples\osnap\data\moor\proc', ...
        moor, [moor 'info.dat']);
    filename = '200044_data';
    addpath(genpath(fullfile(pathgit,'exec','gitrepo')))
elseif strcmp(pc_name,'SA01SJ-G9WC2J3')
    dataindir = 'E:\OSNAP\RHADCP\DY181\S200044A012_RHAD2_JC238\conversion2\';
    pathgit = 'D:\Work_computer_sync\OSNAP_postdoc\Python\m_moorproc_toolbox\';    
    figureoutdir = ['D:\Work_computer_sync\OSNAP_postdoc\Mooring\RHADCP\plots\'];
    addpath(genpath('D:\Work_computer_sync\MATLAB_functions')); % General functions
    filename = 'S200044A012_RHAD2_JC238';
    infofile = '';
    logfile = fullfile(figureoutdir,[moor '_ADCP_stage1.log']);
    addpath(genpath(fullfile(pathgit,'exec','gitrepo')))
else
    error('Please add your path above')
end



if ~isempty(MOORPROC_G)
    addpath(genpath('C:\Users\sa07kb\Matlab\toolboxes\gsw_matlab_v3_06_16'))
    stage01_read_qc_S55(moor)
else
    stage01_read_qc_S55(moor,dataindir,filename,...
                 infofile,logfile,figureoutdir)
end