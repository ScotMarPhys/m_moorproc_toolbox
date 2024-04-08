% process_seaphox.m is a script to process the seaphox data.

% adapt from process_nor.m

close all, 
clearvars -except MEXEC MEXEC_A MEXEC_G pathosnap;
cruise   = 'ar30';
operator = 'lad';
moor = 'rteb1_04_2017';

seaphox2rodb_01(moor,'cruise',cruise)