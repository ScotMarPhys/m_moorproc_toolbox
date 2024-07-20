% process_seaphox.m is a script to process the seaphox data.

% adapt from process_nor.m

close all, 
%clearvars -except MEXEC MEXEC_A MEXEC_G pathosnap;
%cruise   = 'dy181';
%operator = 'tsd';
%moor = 'rteb1_07_2022';

global MOORPROC_G
clearvars -except MEXEC_G MOORPROC_G
cruise   = MOORPROC_G.cruise;
operator = MOORPROC_G.operator;
moor = input('mooring deployment (e.g. ebh2_15_2022) to process:   ','s');


seaphox2rodb_01(moor,'cruise',cruise)

seaphox_raw2use_01(moor,'cruise',cruise)

