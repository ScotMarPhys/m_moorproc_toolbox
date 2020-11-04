% PROCESS_NORS_rb1201.m is a script to process the nortek data.

% 03.Apr.2010 ZB Szuts
% Updated 2012 Feb 15
% To run change basedir ; moor ; inpath
% rteb1_01_2014_filenames.txt must also exist in EB1_recovery directory

close all, clear all

cruise   = 'dy053';
operator = 'loh';
moor = 'rtwb2_02_2015';

plot_interval=[2015 06 01 00, 2016 07 05 00];

if exist('/Volumes/rpdmoc/rapid/data/exec/jc103/stage1/microcat/mc_call_caldip_jc103_v3.m','file')
    % using DR Mac with mount to banba on JC103
    basedir = '/Volumes/rpdmoc/rapid/data/';
else
    basedir = '/home/mstar/osnap/';
end
inpath   = [basedir 'data/moor/raw/' cruise '/nortek/'];
procpath = [basedir 'data/moor/proc/' ];
outpath  = [procpath moor '/nor/'];


nortek2rodb_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)
nortek_raw2use_02(moor,'outpath',outpath,'procpath',procpath,'plot_interval',plot_interval)
