clearvars -except MEXEC MEXEC_A MEXEC_G;
close all;

% pathosnap = '/home/mstar/osnap';

cruise = 'dy120';
castnber = '3'; 

seaphox2rodb_01(['cast' castnber],'cruise',cruise)