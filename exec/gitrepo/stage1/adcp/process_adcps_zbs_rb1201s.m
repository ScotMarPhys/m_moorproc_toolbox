% PROCESS_ADCPS_ZBS_RB1201 is a script to process the adcp data.

% 03.Apr.2010 ZB Szuts

close all, clear all

cruise   = 'rb1201';
operator = 'ad';
moor  = 'wbadcp_8_201118';


basedir  = '/Users/hydrosea5/Desktop/RB1201/rapid/data/moor/';
inpath   = [basedir 'raw/rb1201/adp/'];
procpath = [basedir 'proc/'];
% outpath  = [procpath mooring '/adcp/'];
outpath  = [procpath moor '/adp/'];


adcp2rodb_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)

adcp_raw2use_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)