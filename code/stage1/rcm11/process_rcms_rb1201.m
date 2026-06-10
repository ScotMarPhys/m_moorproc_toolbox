% PROCESS_RCMS_ZBS_RB1201 is a script to process the RCM11 data

% 03.Apr.2010 ZB Szuts
% 13 Apr 2011 EFW
% 15 Feb 2012 AD
close all, clear all

cruise   = 'RB1201';
operator = 'ad';
%moor = 'wb6_5_201117';
moor = 'wb4_8_201115';
moor = 'wb2_9_201114';
%moor = 'wbh2_4_201004';
%moor = 'wb4_7_201026';
%moor = 'wb6_4_201001';

basedir  = '/Users/hydrosea5/Desktop/RB1201/rapid/data/moor/';
%basedir  = '/Volumes/RB1201/rapid/data/moor/';
inpath   = [basedir 'raw/rb1201/rcm/'];
%procpath = [basedir 'proc_kn200_4/'];  % modifications !!
procpath = [basedir 'proc/'];
outpath  = [procpath moor '/rcm/'];


rcm2rodb_05(moor,'procpath',procpath,'inpath',inpath,...
                  'outpath',outpath);

rcm11raw2use(moor,'procpath',procpath,'outpath',outpath);

