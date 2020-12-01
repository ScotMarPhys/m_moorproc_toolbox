% function [bottle_data] = read_rosfile(ros_file) 
%
% read data from seabird ctd rosette file
%  output:  
%        bottle_data (struct array) containing:
%        jd     --- time of bottle marker [julian days]
%        p      --- pressure [dbar]
%        s      --- salinity [psu] 
%        t0,t1  --- temperature from primary/secondary sensor [�C]   
%        c0,c1  --- conductivity from primary/secondary sensor [mS/cm]
%  
%        variable extensions exist for p,t0,t1,c0,c1,s,ox,tm:
%           extensions: av(average),sn(standard deviation),mn(minimum),mx(maximum) 
%           refer to duration  of each bottle stop                 
% kanzow 01.08.05   
%
% uses hms2h.m (kanzow)
%   
% modfied 04.07.2006 Kanzow: distinction between timeS and timeJ
% modfied 19.10.2015 L. Houpert: distinction between timeQ (NMEA) and timeJ, and add 
%  a correction if the start_time in the .ros file is bad due to a problem in the NMEA frame (e.g. PE399)
% modified 19.05.2016 L. Houpert: fix bug when timeJ and timeS are present
% at the same time
% 

function bottle = read_rosfile(ros_file) 
 vartime=[];
  dummy = 9999;
  timeS = 'timeS';
  timeJ = 'timeJ';
  timedef = 0 ; 
fid1 = fopen(ros_file,'r');


disp(['loading ',ros_file,'  ...'])
zeile = fscanf(fid1,'%c');  %read data into string

fclose(fid1);
disp(['loading complete'])

ret = sprintf('\n');
retx = findstr(zeile,ret);  % car. return indices

if length(retx)< 10
  disp('input file does not contain data ')
  fprintf(fidlog,'input file does not contain data \n')
  return
end

%---------------------------------------------------------
%           scan header
%--------------------------------------------------------

vars = ['ti';'pr';'t0';'c0';'sa';'t1';'c1';'ox';'fl';'sb';'de';'po';'si';'sv';'T2';'C2';'pt';'la';'lo'];

i = 1;
while 1
  name =  ['name ',num2str(i-1)];
  ln   = length(name);
  a = findstr(zeile,name); 
  if isempty(a)
    break
  end

  for  j = 1 : size(vars,1)
    varstr =   zeile(a+ln+3:a+ln+4);
    if  strcmp(varstr,vars(j,1:2)) == 1
     if strcmp(varstr,'ti')  
       if strcmp(varstr,'ti') & strcmp(zeile(a+ln+3:a+ln+7),timeS)
          varn(i) = j;
          vartime{i} = timeS;            
          timedef = 1;
       elseif strcmp(varstr,'ti') & strcmp(zeile(a+ln+3:a+ln+7),timeJ)
          varn(i) = j;
          vartime{i} = timeJ;   
          timedef = 1;   
       else 
        varn(i) = 0;           
       end
     else 
        varn(i) = j; 
     end  
    end
  end     % end for
  i = i + 1;
end     % end while

if timedef == 0 
   disp('Check time variable used in read_rosfile.m')
end
ii         = find(varn==0);
varn(ii)   = dummy;

%datainfo.variables = vars(varn,:);

ii         = findstr(zeile,'start_time');

start_time = zeile(ii+13:ii+32);
start_time   = datevec(datenum(start_time));
start_time   = [start_time(1:3) hms2h(start_time(4:6))];
start_time   = julian(start_time);

jj         = findstr(zeile,'NMEA UTC (Time)');
nmea_time  = zeile(jj+18:jj+37);
nmea_time   = datevec(datenum(nmea_time));
nmea_time   = [nmea_time(1:3) hms2h(nmea_time(4:6))];
nmea_time   = julian(nmea_time);


% ---- extract data   ---------------- 


header_end = findstr(zeile,'*END*');
ii         = find(retx>header_end);
data_begin = retx(ii(1));

data       = str2num(zeile(data_begin:end));


for i = 1: 9
  ii               = find(varn == i);
  if ~isempty(ii)
        eval([vars(i,:),'= data(:,ii(1));']) 
        if i == 1
            vartimeselec = vartime{ii(1)};
        end
  else
    eval([vars(i,:),'= NaN;'])
  end    
end

if abs(start_time-nmea_time) >10 % When pb in NMEA time, (e.g. cruise PE399)
    start_time = nmea_time; 
    gtime = gregorian(start_time);
    t00 = julian(gtime(1),gtime(2),gtime(3));
    jd         = t00 + (ti-1);
else
  if strcmp(vartimeselec,timeJ) %time format : timeJ
    gtime= gregorian(start_time);
    jd = julian(gtime(1),1,1) + ti - 1;
  else  
    jd         = start_time + ti/86400;
  end
end
pos        = 1:length(jd);



bottle     = struct('pos',pos,'jd',jd,'p',pr,'s',sa,'t0',t0,'c0',c0,'t1',t1,'c1',c1,'start_time',start_time);







