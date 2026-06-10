% PROCESS_ADCPS_ZBS_OC459 is a script to process the adcp data.

% 03.Apr.2010 ZB Szuts

close all, clear all

cruise   = 'oc459';
operator = 'zszuts';
mooring  = 'wbadcp_6_200909';


basedir  = '/Volumes/surman/rpdmoc/rapid/data/moor/';
inpath   = [basedir 'raw/oc459/adcp/'];
procpath = [basedir 'proc/'];
outpath  = [procpath mooring '/adcp/'];


adcp2rodb_01(mooring,'inpath',inpath,'outpath',outpath,'procpath',procpath)




