% PROCESS_ADCPS_ZBS_RB1201 is a script to process the adcp data.

% 03.Apr.2010 ZB Szuts

close all, clear all

cruise   = 'dy120';
operator = 'lad_loh';
% moor  = 'nocm1_03_2016'; 
% moor  = 'nocm2_03_2016'; 
% moor  = 'nocm3_03_2016'; 
% moor  = 'nocm4_03_2016'; 
% moor  = 'ib5_01_2018'; 
% moor  = 'ib4_01_2018'; 
moor  = 'ib3_01_2018'; 

basedir  = '/local/users/pstar/osnap/';


inpath   = [basedir 'data/moor/raw/' cruise '/adcp/'];
procpath = [basedir 'data/moor/proc/'];
% outpath  = [procpath mooring '/adcp/'];
outpath  = [procpath moor '/adp/'];

infofile = [basedir 'data/moor/proc/' moor '/' moor 'info.dat'];

% Flag bad data 
read_flag_raw_adcp(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)

% Stage 1
display('Starting stage 1 adcp2rodb_02')
adcp2rodb_02(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)

display('Stage 1 complete')

% Stage 2
display('Starting stage 2 adcp_raw2use')
adcp_raw2use_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)
