% D324. to run stage 2 s4 data processing. mpc

moor   = 'mar1_3_200621'; % change mooring
sn_s4  = '35612573'; % serial number
stage  = '2'; % stage 1 or 2
operator = 'mpc';
channels=5;
SRBchannels=5;

% mpc define paths for d324: 
%inpath   = ['/data32/d324/rapid/data/moor/raw/d324/s4/'];
%inpath   = [inpath1,'',num2str(sn_s4),'_data.cap']
inpath   = ['/data32/d324/rapid/data/moor/proc/',moor,'/s4/'];
outpath  = ['/data32/d324/rapid/data/moor/proc/',moor,'/s4/'];
%outpath  = [outpath1,'',moor,'_',num2str(sn_s4),'.raw']
infofile = ['/data32/d324/rapid/data/moor/proc/',moor,'/',moor,'info.dat'];
log = [outpath,'s4_stage',num2str(stage),'.log']

plot_interval = [2006 05 25 18;   % start time of time axis on plot
		 2007 10 25 8];  % end time of time axis on plot

s4raw2use_v2(moor,inpath,outpath,infofile,plot_interval);
