% PROCESS_RCMS_ZBS_OC459 is a script to process the RCM11 data

% 03.Apr.2010 ZB Szuts
% 13 Apr 2011 EFW
close all, clear all

cruise   = 'kn200_4';
operator = 'efw';
moor = 'wb1_7_201008';
moor = 'wb2_8_201003';
%moor = 'wbh2_4_201004';
%moor = 'wb4_7_201026';
%moor = 'wb6_4_201001';


basedir  = '/Users/hydrosea5/kn200-4/rpdmoc/rapid/data/moor/';
inpath   = [basedir 'raw/kn200_4/rcm/'];
procpath = [basedir 'proc_kn200_4/'];
outpath  = [procpath moor '/rcm/'];


rcm2rodb_05_kn200(moor,'procpath',procpath,'inpath',inpath,...
                  'outpath',outpath);

rcm11raw2use(moor,'procpath',procpath,'outpath',outpath);

