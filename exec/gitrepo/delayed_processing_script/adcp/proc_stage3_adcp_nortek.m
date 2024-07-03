% Script that load the stage3 processing routine cm_edit_NOCS_v3 to correct  (on all the different 
% currentmetes of a specific mooring) : the magnetic deviation, remove spikes, produce a low pass filter dataset, 
% Also possibility of manual QC 
% 
%       by loh, SAMS, 04/01/2016


close all
global MOORPROC_G
clearvars -except MEXEC_G MOORPROC_G
cruise   = MOORPROC_G.cruise;
operator = MOORPROC_G.operator;
moor = input('mooring deployment (e.g. ebh2_15_2022) to process:   ','s');
plot_interval=[]; %automatic based on available times

cm_edit_NOCS(moor)

