% PROCESS_ADCPS_ZBS_RB1201 is a script to process the adcp data.

% 03.Apr.2010 ZB Szuts

close all, clear all

cruise   = 'dy120';
operator = 'lad';
moor  = 'ib5_01_2018';


basedir  = '/local/users/pstar/osnap/data/moor/';
inpath   = [basedir 'raw/' cruise '/adcp/'];
procpath = [basedir 'proc/'];
% outpath  = [procpath mooring '/adcp/'];
outpath  = [procpath moor '/adcp/'];


display('Starting stage 1 adcp2rodb_01')
adcp2rodb_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)

display('Stage 1 complete')
display('Starting stage 2 adcp_raw2use')
adcp_raw2use_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)
