% process_seaphox.m is a script to process the seaphox data.

% adapt from process_nor.m

close all, 
clearvars -except MEXEC MEXEC_A MEXEC_G pathosnap;
cruise   = 'dy120_test';
operator = 'lad';
moor = 'test_00_0000';

seaphox2rodb_01(moor,'cruise',cruise)