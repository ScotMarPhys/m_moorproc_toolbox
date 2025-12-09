%%  Code for the quality control of Signature 55 ADCP

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
% - praembel for despiking
% - if you add a deployment period, you need to edit steps 1-3
%---------
close all
clear all

addpath(genpath('../functions'));

% TO UPDATE <------------------------
% in- and output directories
pc_name = getenv('COMPUTERNAME');
if strcmp(pc_name,'SA07KB-3JN9YY2');
    basedir = 'C:\Users\sa07kb\OneDrive - SAMS\data\data_OSNAP\German_53N_Array\';
    pathgit = 'C:\Users\sa07kb\Projects\Oi\';
elseif strcmp(pc_name,'SA01SJ-G9WC2J3')
    basedir = 'E:\Oi\general\data_OSNAP\From_Fehmi\';
    pathgit = 'D:\Work_computer_sync\OSNAP_postdoc\Python\Oi\OSNAPi\';
    addpath(genpath('D:\Work_computer_sync\MATLAB_functions'));
    datapath = addpath(genpath('D:\Work_computer_sync\OSNAP_postdoc\Python\Oi\Additional_OSNAPi_datasets'));
else
    error('Please add your path above')
end
grdatdir = [pathgit 'data\interim\mooring\53n_array\hydro\'];