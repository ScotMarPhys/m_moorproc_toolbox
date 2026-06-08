%seaguard_raw2use_call.m
%batch file for calling seaguard_raw2use for more than one file.
%seaguard_raw2use(Serial_number,infile,infofile,varargin)
%seaguard_raw2use('114','/local/users/pstar/Data/rpdmoc/rapid/data/moor/proc_calib/d344/cal_dip/seaguard/cast4/calib4_114.raw','/local/users/pstar/Data/rpdmoc/rapid/data/moor/proc_calib/d344/cal_dip/cast4info.dat')

% serial_number='114';
% infile = '/noc/users/pstar/rpdmoc/rapid/data/moor/proc/ebhi_6_200934/seaguard/ebhi_6_200934_114.raw';
% infofile = '/noc/users/pstar/rpdmoc/rapid/data/moor/proc/ebhi_6_200934/ebhi_6_200934info.dat';
% outpath ='/noc/users/pstar/rpdmoc/rapid/data/moor/proc/ebhi_6_200934/seaguard/';

serial_number='114';
infile = '/Users/hydrosea5/kn200-4/rpdmoc/rapid/data/moor/proc/ebhi_6_200934/seaguard/ebhi_6_200934_114.raw';
infofile = '/Users/hydrosea5/kn200-4/rpdmoc/rapid/data/moor/proc/ebhi_6_200934/ebhi_6_200934info.dat';
outpath ='/Users/hydrosea5/kn200-4/rpdmoc/rapid/data/moor/proc/ebhi_6_200934/seaguard/';


seaguard_raw2use(serial_number,infile,infofile,'outpath',outpath)

