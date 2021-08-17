clearvars -except MEXEC MEXEC_A MEXEC_G;
close all;

% pathosnap = '/home/mstar/osnap';

cruise = 'ar304';
castnber = '6'; 

seaphox2rodb_01(['cast' castnber],'cruise',cruise)