% PROCESS_NORS_rb1201.m is a script to process the nortek data.

% 03.Apr.2010 ZB Szuts
% Updated 2012 Feb 15

close all, clear all

cruise   = 'rb1201';
operator = 'ad';
moor = 'wb1_8_201113';
%moor = 'wb2_9_201114';
%moor = 'wbh2_5_201116';
%moor = 'wb4_8_201115';
%moor = 'wb2_8_201003';
%moor = 'wbh2_4_201004';
%moor = 'wb4_7_201026';
%moor = 'wb6_4_201001';

basedir  = '/Users/hydrosea5/Desktop/RB1201/rapid/data/moor/';
inpath   = [basedir 'raw/rb1201/nor/'];
procpath = [basedir 'proc/'];
outpath  = [procpath moor '/nor/'];

%nortek2rodb_01(mooring,'inpath',inpath,'outpath',outpath,'procpath',procpath)
%nortek_raw2use_01(mooring,'outpath',outpath,'procpath',procpath)

nortek2rodb_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)
nortek_raw2use_01(moor,'outpath',outpath,'procpath',procpath)