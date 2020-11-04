% PROCESS_NORS_rb1201.m is a script to process the nortek data.

% 03.Apr.2010 ZB Szuts
% Updated 2012 Feb 15
% To run change basedir ; moor ; inpath
% rteb1_01_2014_filenames.txt must also exist in EB1_recovery directory

close all, 
clearvars -except MEXEC MEXEC_A MEXEC_G pathosnap;
cruise   = 'ar30';
operator = 'sj';
%moor = 'rtwb2_04_2017';
%moor = 'rtwb1_04_2017';
moor = 'rteb1_04_2017';

plot_interval=[2017 05 10 00, 2018 07 10 00];

if exist('pathosnap','var')
    basedir = [pathosnap filesep 'data' filesep];
else
    basedir = '/home/mstar/osnap/data/';
end


inpath   = [basedir 'moor/raw/' cruise '/nortek/'];
procpath = [basedir 'moor/proc/' ];
outpath  = [procpath moor '/nor/'];


nortek2rodb_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)
nortek_raw2use_02(moor,'outpath',outpath,'procpath',procpath,'plot_interval',plot_interval)
