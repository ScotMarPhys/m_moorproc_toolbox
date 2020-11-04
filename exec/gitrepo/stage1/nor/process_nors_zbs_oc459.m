% PROCESS_NORS_ZBS_OC459 is a script to process the nortek data.

% 03.Apr.2010 ZB Szuts

close all, clear all

cruise   = 'oc459';
operator = 'zszuts';
mooring  = 'wbh2_3_200912';


basedir  = '/Volumes/surman/rpdmoc/rapid/data/moor/';
inpath   = [basedir 'raw/oc459/nor/'];
procpath = [basedir 'proc/'];
outpath  = [procpath mooring '/nor/'];


nortek2rodb_01(mooring,'inpath',inpath,'outpath',outpath,'procpath',procpath)

nortek_raw2use_01(mooring,'outpath',outpath,'procpath',procpath)
