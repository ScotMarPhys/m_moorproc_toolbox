% Script that load the stage3 processing routine cm_edit_NOCS_v3 to correct  (on all the different 
% currentmetes of a specific mooring) : the magnetic deviation, remove spikes, produce a low pass filter dataset, 
% Also possibility of manual QC 
% 
%       by loh, SAMS, 04/01/2016


close all, 
clearvars  -except pathosnap datadir execdir

%moor = 'rtwb1_02_2015';
% moor = 'rteb1_05_2018';
% moor = 'rtwb1_05_2018';
moor = 'rtwb2_05_2018';
% 
% if exist('/Volumes/rpdmoc/rapid/data/exec/jc103/stage1/microcat/mc_call_caldip_jc103_v3.m','file')
%     % using DR Mac with mount to banba on JC103
%     basedir = '/Volumes/rpdmoc/rapid/data/';
% elseif exist('/home/mstar/osnap/','dir')
%     basedir = '/home/mstar/osnap/';
%  elseif exist('m:/Mar_Phys/OSNAP_mooring_data_processing/osnap/','dir')
%     basedir = 'm:/Mar_Phys/OSNAP_mooring_data_processing/osnap/';   
% end
% 
% procpath = [basedir 'data/moor/proc/' ];

procpath = [pathosnap '/data/moor/proc/' ];

cm_edit_NOCS_v5(moor,'procpath',procpath)

