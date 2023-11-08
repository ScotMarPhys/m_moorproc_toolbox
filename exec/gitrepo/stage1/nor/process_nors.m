ed% PROCESS_NORS_rb1201.m is a script to process the nortek data.

% 03.Apr.2010 ZB Szuts
% Updated 2012 Feb 15
% To run change basedir ; moor ; inpath
% rteb1_01_2014_filenames.txt must also exist in EB1_recovery directory

close all
global MOORPROC_G
clearvars -except MEXEC_G MOORPROC_G
cruise   = MOORPROC_G.cruise;
operator = MOORPROC_G.operator;
moor='wb1_16_2023a';
plot_interval=[2023 02 07 0, 2023 07 21 0];


inpath   = fullfile(MOORPROC_G.moordatadir,'raw',cruise,'nor');%'nortek'
procpath = fullfile(MOORPROC_G.moordatadir,'proc');
outpath  = fullfile(procpath,moor,'nor');
if ~exist(outpath,'dir')
    mkdir(outpath)
end

nortek2rodb_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)

nortek_raw2use_02(moor,'outpath',outpath,'procpath',procpath,'plot_interval',plot_interval)
