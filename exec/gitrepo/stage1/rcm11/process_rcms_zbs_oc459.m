% PROCESS_RCMS_ZBS_OC459 is a script to process the RCM11 data

% 03.Apr.2010 ZB Szuts

close all, clear all

cruise   = 'oc459';
operator = 'zszuts';
mooring  = 'wb2_7_200907';
%mooring  = 'wb1_6_200906';


basedir  = '/Volumes/surman/rpdmoc/rapid/data/moor/';
inpath   = [basedir 'raw/oc459/rcm/'];
procpath = [basedir 'proc/'];
outpath  = [procpath mooring '/rcm/'];


rcm2rodb_05_oc459(mooring,'procpath',procpath,'inpath',inpath,...
                  'outpath',outpath);

rcm11raw2use(mooring,'procpath',procpath,'outpath',outpath);

