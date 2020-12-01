function  [ti,dti] = ctd_impact(file,var,threshold,cruise)

% function  [ti,dti] =ctd_impact(file,[var],[threshold])
%
% find time when seabird ctd first hits the surface
%
% input:  file       --- ctd data file (raw,'.CNV')
%         var        --- 'c' /  'p' for variable used to 
%                        identify impact [default 'c']
%
%         threshold  --- threshold value of variable 'var'       
%                        to identify impact [default 4 S/m]
%
% 
%
% output: ti         --- ctd in water time [year month day hour min sec]
%         dti        --- time that elapsed from begin of data acquisition
%                        till 'ti'   
%
% kanzow 28.3.01
%     uses microcat_month.m (kanzow)
%          hms2h.m (kanzow) 
% Loic Houpert 15.10.2015
%     add condition to fix file with a timestamp problem for the header variable #start_time

% --- setting --  
ctime = 8; % number of time column
timeJdata=0;

if nargin == 1
   
  var        = 'c';
  threshold = 30;  % sea surface conductivity (30(T=0) - 50(T=25)) 
 
elseif nargin == 2
  
   threshold = 30;  % sea surface conductivity (30(T=0) - 50(T=25)) 
end 


% --- get start time from header  

fid = fopen(file,'r');

if fid < 0
 ti = NaN; dti = NaN;
 disp(['could not open ',file]) 
 return 
end

%%vars = ['ti';'pr';'t0';'sa';'c0';'t1';'c1';'ox';'fl';'sb'];
vars = ['ti';'pr';'t0';'sa';'c0';'t1';'c1';'ox';'fl';'sb';'de';'po';'si';'sv';'T2';'C2';'pt';'pu';'sc';'ba';'up';'la';'lo';'al'];

i = 1;

%%% fix to add path correctly for microcat_month
% addpath('/noc/mpoc/rpdmoc/rapid/data/exec/d382/stage1/microcat/preOC459/historical/')

while 1 
  zeile = fgetl(fid);
  if ~isempty(findstr(zeile,'# start_time ='))  
    st = zeile(16:length(zeile));
    ii = find(st ~=' ');
    ti(2)   = microcat_month(st(ii(1:3)));
    ti(3)   = str2num(st(ii(4:5)));
    ti(1)   = str2num(st(ii(6:9)));
    ti(4)   = hms2h(st(ii(10:17)));  
 

    ii =0;
    start_time = zeile(ii+13:ii+32);
    start_time   = datevec(datenum(start_time));
    start_time   = [start_time(1:3) hms2h(start_time(4:6))];
    start_time   = julian(start_time);
  end
  
   if ~isempty(findstr(zeile,'NMEA UTC (Time)'));  
    jj=0;
    nmea_time  = zeile(jj+18:jj+37);
    nmea_time   = datevec(datenum(nmea_time));
    nmea_time   = [nmea_time(1:3) hms2h(nmea_time(4:6))];
    nmea_time   = julian(nmea_time);

   end
     
  
 name =  ['name ',num2str(i-1)];
  ln = length(name);
  a = findstr(zeile,name);
  
  if ~isempty(a)
    for  j = 1 : size(vars,1)
      %if  strcmp(zeile(12:13),vars(j,1:2)) == 1  
      varstr = zeile(a+ln+3:a+ln+4); 
      if  strcmp(varstr ,vars(j,1:2)) == 1  
        varn(i) = j;
        if strcmp(zeile(a+ln+3:a+ln+7),'timeJ')
            varn(i) = 200;
            timeJdata = 1;
        end
        i = i + 1  ;
        break
      end
    %%%i=i+1;
    end     % end for
  end


  if ~isempty(findstr(zeile,'*END'))
    break
  end
end  % end while

if timeJdata == 0
ctime =find(varn==1); % time column 
else
ctime= find(varn==200); % timeJ data format
end


if abs(start_time-nmea_time) >10 % When pb in start NMEA time, (e.g. cruise PE399)
    ti = gregorian(nmea_time);
end

%ctime =find(varn==1); % time column 

ctime = ctime(1);
if strcmp(var,'p')
  ccond = find(varn == 2); % pressure column
elseif strcmp(var,'c') 
  ccond = find(varn == 5); % pressure column
elseif strcmp(var,'s') 
  ccond = find(varn == 4); % pressure column 
end


% ---- get time of ctd impact
sttime='';
while 1
  
  zeile = fgetl(fid);
  if zeile == -1 
    disp('no ctd impact found -- please check if column in .cnv file is correct')       
    break
  end
  data  = str2num(zeile);
  if isempty(sttime)
      sttime=data(ctime);
  end
  if data(ccond) > threshold  % ctd impact!!!
    dti    = data(ctime);
    break
  end 
end

if timeJdata == 1
    dti=dti-sttime;
end


ti(4) = ti(4) + dti/3600; 

ti    = gregorian(julian([ti]));

if length(ti) <5      % if gregorian output is [yy mm dd hh] 
                      % convert to [yy mm dd hh min sec] 
  [ti(4) ti(5) ti(6)] = h2hms(ti(4));
end


  


