% Script that load the stage3 processing routine cm_edit_NOCS_v3 to correct  (on all the different 
% currentmetes of a specific mooring) : the magnetic deviation, remove spikes, produce a low pass filter dataset, 
% Also possibility of manual QC 
% 
%       by loh, SAMS, 04/01/2016


close all
clearvars -except MEXEC MEXEC_A MEXEC_G pathosnap;

%moor = 'rtwb1_01_2014';
%moor = 'rtwb2_04_2017';
%moor = 'rtwb1_04_2017';
%moor  = 'ib5_01_2018';
moor  = 'ib3_01_2018';
moor  = 'rteb1_05_2018';

if exist('pathosnap','var')
    basedir = [pathosnap filesep 'data' filesep];
else
    basedir = '/local/users/pstar/osnap/data/';
end

procpath = [basedir 'moor/proc/' ];


cm_edit_NOCS(moor,'procpath',procpath)

