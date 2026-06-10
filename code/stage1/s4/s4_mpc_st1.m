%D324 to run stage 1 s4 data processing. mpc
clear all;close all;
moor   = 'mar1_3_200621'; % change mooring
sn_s4  = '35612573'; % serial number
stage  = '1'; % stage 1 or 2
operator = 'mpc';
channels=5;
SRBchannels=5;

% mpc define paths for d324: 
inpath   = ['/data32/d324/rapid/data/moor/raw/d324/s4/'];
infile   = [inpath,'',num2str(sn_s4),'_data.cap'];
outpath  = ['/data32/d324/rapid/data/moor/proc/',moor,'/s4/'];
outfile  = [outpath,'',moor,'_',num2str(sn_s4),'.raw'];
infofile = ['/data32/d324/rapid/data/moor/proc/',moor,'/',moor,'info.dat'];
log = [outpath,'s4_stage',num2str(stage),'.log'];

s42rodb_v5(infile,outfile,infofile,log,channels,SRBchannels);

%function s42rodb_v5('/data32/d324/rapid/data/moor/raw/d324/s4/35612567_data.cap','/data32/d324/rapid/data/moor/proc/mar3_3_200623/s4/mar3_3_200623_35612567.raw','/data32/d324/rapid/data/moor/proc/mar3_3_200623/mar3_3_200623info.dat','/data32/d324/rapid/data/moor/proc/mar3_3_200623/s4/s4_stage1.log',5,5)

