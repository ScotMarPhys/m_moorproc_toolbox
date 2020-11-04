% Batch file for processing all the S4 files over the weekend.

addpath /data/jrd/hydro10/rapid/data/exec/moor/stage1/s4/

s42rodb_v2('/data/jrd/hydro10/rapid/data/moor/raw/WB1_2005/s4/2565data.S4A','/data/jrd/hydro10/dr400/WB1/35612565.raw','/data/jrd/hydro10/rapid/data/moor/proc/wb1_2_200527/wb1_2_200527info.dat','/data/jrd/hydro10/dr400/WB1/35612565.log',5,5);
s42rodb_v2('/data/jrd/hydro10/rapid/data/moor/raw/WB1_2005/s4/2568data.S4A','/data/jrd/hydro10/dr400/WB1/35612568.raw','/data/jrd/hydro10/rapid/data/moor/proc/wb1_2_200527/wb1_2_200527info.dat','/data/jrd/hydro10/dr400/WB1/35612568.log',5,5);
s42rodb_v2('/data/jrd/hydro10/rapid/data/moor/raw/WB1_2005/s4/2569data.S4A','/data/jrd/hydro10/dr400/WB1/35612569.raw','/data/jrd/hydro10/rapid/data/moor/proc/wb1_2_200527/wb1_2_200527info.dat','/data/jrd/hydro10/dr400/WB1/35612569.log',5,5);
s42rodb_v2('/data/jrd/hydro10/rapid/data/moor/raw/WB1_2005/s4/2570data.S4A','/data/jrd/hydro10/dr400/WB1/35612570.raw','/data/jrd/hydro10/rapid/data/moor/proc/wb1_2_200527/wb1_2_200527info.dat','/data/jrd/hydro10/dr400/WB1/35612570.log',5,5);
