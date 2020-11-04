% PROCESS_NORS_kn200.m is a script to process the nortek data.

% 03.Apr.2010 ZB Szuts
% Updated 2011 Apr 13

close all, clear all

cruise   = 'kn200_4';
operator = 'efw';
moor = 'wb1_7_201008';
moor = 'wb2_8_201003';
moor = 'wbh2_4_201004';
moor = 'wb4_7_201026';
moor = 'wb6_4_201001';

basedir  = '/Users/hydrosea5/kn200-4/rpdmoc/rapid/data/moor/';
inpath   = [basedir 'raw/kn200_4/nor/'];
procpath = [basedir 'proc/'];
outpath  = [procpath mooring '/nor/'];


nortek2rodb_01(mooring,'inpath',inpath,'outpath',outpath,'procpath',procpath)

nortek_raw2use_01(mooring,'outpath',outpath,'procpath',procpath)
