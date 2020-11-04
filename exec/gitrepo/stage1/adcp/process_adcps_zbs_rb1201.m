% PROCESS_ADCPS_ZBS_RB1201 is a script to process the adcp data.

% 03.Apr.2010 ZB Szuts

close all, clear all

cruise   = 'rb1201';
operator = 'zszuts';
mooring  = 'wbadcp_6_200909';


basedir  = '/Volumes/RB1201/rapid/data/moor/';
inpath   = [basedir 'raw/rb1201/adp/'];
procpath = [basedir 'proc/'];
% outpath  = [procpath mooring '/adcp/'];
outpath  = [procpath mooring '/adp/'];


adcp2rodb_01(mooring,'inpath',inpath,'outpath',outpath,'procpath',procpath)




