% PROCESS_RCMS_ZBS_RB1201 is a script to process the RCM11 data

% 03.Apr.2010 ZB Szuts
% 13 Apr 2011 EFW
% 15 Feb 2012 AD
close all, clear all
fclose all % DR added as need to close open info.dat file if reads when not completely filled in as subsequent runs of the routine continue to read the same text as before the file was changed

cruise   = 'jc103';
operator = 'dr400';
moor = 'ebh3_9_201233';
%moor = 'wb6_5_201117';
%moor = 'wb4_8_201115';
%moor = 'wb2_9_201114';
%moor = 'wbh2_4_201004';
%moor = 'wb4_7_201026';
%moor = 'wb6_4_201001';


if exist('/Volumes/rpdmoc/rapid/data/exec/jc103/stage1/microcat/mc_call_caldip_jc103_v3.m','file')
    % using DR Mac with mount to banba on JC103
    basedir = '/Volumes/rpdmoc/rapid/data/';
else
    basedir = '/local/users/pstar/rpdmoc/rapid/data/';
end
%basedir  = '/Users/hydrosea5/Desktop/RB1201/rapid/data/moor/';
%basedir  = '/Volumes/RB1201/rapid/data/moor/';
%basedir
inpath   = [basedir 'moor/raw/' cruise '/rcm/'];
%procpath = [basedir 'proc_kn200_4/'];  % modifications !!
procpath = [basedir 'moor/proc/'];
outpath  = [procpath moor '/rcm/'];


rcm2rodb_05(moor,'procpath',procpath,'inpath',inpath,...
                  'outpath',outpath);

rcm11raw2use(moor,'procpath',procpath,'outpath',outpath);

