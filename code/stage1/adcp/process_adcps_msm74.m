% PROCESS_ADCPS_ZBS_RB1201 is a script to process the adcp data.

% 03.Apr.2010 ZB Szuts

close all, clear all

cruise   = 'msm74';
operator = 'nph';
% moor  = 'nocm1_03_2016'; 
% moor  = 'nocm2_03_2016'; 
% moor  = 'nocm3_03_2016'; 
% moor  = 'nocm4_03_2016'; 
moor  = 'nocm5_03_2016'; 


% basedir  = '/home/mstar/osnap/data/moor/';
% MSM74
basedir = '/Users/ukosnap/Documents/aaaaMSM74/new_mooring_processing/moorings/';


inpath   = [basedir 'data/moor/raw/' cruise '/adcp/'];
procpath = [basedir 'data/moor/proc/'];
% outpath  = [procpath mooring '/adcp/'];
outpath  = [procpath moor '/adp/'];

infofile = [basedir 'data/moor/proc/' moor '/' moor 'info.dat'];

display('Starting stage 1 adcp2rodb_01')
adcp2rodb_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)

display('Stage 1 complete')
display('Starting stage 2 adcp_raw2use')
adcp_raw2use_01(moor,'inpath',inpath,'outpath',outpath,'procpath',procpath)
