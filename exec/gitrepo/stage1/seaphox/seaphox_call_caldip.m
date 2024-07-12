%clearvars -except MEXEC MEXEC_A MEXEC_G;
close all;

% pathosnap = '/home/mstar/osnap';

% === TSD Jul 2024: modified here to have global MOORPROC_G
global MOORPROC_G

basedir = MOORPROC_G.moordatadir;
cruise= MOORPROC_G.cruise;

% cruise = 'ar304';
castnber = '6'; 

seaphox2rodb_01(['cast' castnber],'cruise',cruise)