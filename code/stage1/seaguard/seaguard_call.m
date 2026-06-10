% seaguard_call.m
% batch file to read Seaguard curret meter data into rodb format using
% seaguard2rodb function.

inpath = '/Users/hydrosea5/kn200-4/rpdmoc/rapid/data/moor/proc/ebhi_6_200934/ebhi_6_200934info.dat';
%seaguard2rodb('/local/users/pstar/Data/rpdmoc/rapid/data/moor/proc_calib/d344/cal_dip/cast4info.dat','d344',1);
seaguard2rodb(inpath,'d359',0)