% process_seaphox.m is a script to process the seaphox data.

% adapt from process_nor.m

close all, 
clearvars -except MEXEC MEXEC_A MEXEC_G pathosnap;
cruise   = 'dy120';
operator = 'lad';
moor = 'rteb1_05_2018';

seaphox2rodb_01_indev(moor,'cruise',cruise)