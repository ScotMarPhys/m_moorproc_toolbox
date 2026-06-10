function [P,T,C,S,jd,meas,sampling_rate,qflag] = load_seacat(file,log);
% function [P,T,C,S,jd,meas,sampling_rate,qflag] = load_seacat(file,log);
%
% uses julian hms2h

% Kanzow
 %%slsh = '/';
 %%cln = ':';
 begin_mark = '*END*';
 ret = sprintf('\n');
 
% ------ load data ------------

fid = fopen(file,'r');
disp('load_seacat.m: loading data')
zeile = fscanf(fid,'%c');  %read data into string

disp('complete')
fclose(fid);


% ------- convert data ---------

%nonum = findstr(zeile,slsh); 
%zeile(nonum) = ' ' ;
%nonum = findstr(zeile,cln);
%zeile(nonum) = ' ' ;

markI = findstr(zeile,begin_mark);
retI  = findstr(zeile,ret);

begin = find(retI>markI);
zeile = zeile(retI(begin(1))+1:end);

data = str2num(zeile);

[m,n] = size(data);

data(1,:)

jc = input('Julian day column ');
reference_year = input('Reference Year ');
pc = input('Pressure column ');
tc = input('Temperature column ');
cc = input('Conductivity column ');

%%jd   = julian([data(:,[4 2 3]) hms2h(data(:,[5 6 7]))]);
jd   = data(:,jc) + julian(reference_year,1,1,0);
P    = data(:,pc); % psi to dbar conversion;
T    = data(:,tc);
C    = data(:,cc);
meas = data(:,1);
S    = NaN;

%if n == 9   % no conductivity
% C = NaN;
% S = NaN;
%elseif n>9
% C = data(:,10);
% S = data(:,11);
%end


sampling_rate = round(1/median(diff(jd))); 
qflag = [];

 