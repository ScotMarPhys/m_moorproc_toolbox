% PROCESS_NORS.m is a script to process the nortek data.

% 03.Apr.2010 ZB Szuts
% Updated 2012 Feb 15
% To run change basedir ; moor ; inpath
% rteb1_01_2014_filenames.txt must also exist in EB1_recovery directory

close all
global MOORPROC_G
clearvars -except MEXEC_G MOORPROC_G
cruise   = MOORPROC_G.cruise;
operator = MOORPROC_G.operator;
moor = input('mooring deployment (e.g. ebh2_15_2022) to process:   ','s');
%moor='wb1_16_2023a';
%plot_interval=[2023 02 07 0, 2023 07 21 0];
plot_interval=[]; %automatic based on available times

pd = moor_inoutpaths('nor',moor);

nortek2rodb_01(moor,pd)

nortek_raw2use_02(moor,pd,'plot_interval',plot_interval)
