% PROCESS_ADCPS_ZBS_RB1201 is a script to process the adcp data.

% 03.Apr.2010 ZB Szuts

close all, clear all

cruise   = 'd382';
operator = 'das';
moor  = 'wbadcp_9_201206';


basedir  = '/noc/users/pstar/rpdmoc/rapid/data/moor/';
inpath   = [basedir 'raw/d382/adcp/'];
procpath = [basedir 'proc/'];
% outpath  = [procpath mooring '/adcp/'];
outpath  = [procpath moor '/adp/'];


dislpay('Startign stage 1 adcp2rodb_01')
adcp2rodb_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)

display('Stage 1 complete')
display('Stage 2 adcp_raw2use')
adcp_raw2use_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)
