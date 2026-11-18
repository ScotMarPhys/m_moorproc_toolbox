% function bottle = read_rosfile(ros_file) 
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
% Bug fix: 6th Dec 20220, L. Drysdale Julian function was producing 12 hour offset
% on start time of .ros file, CTD cast 19, JC238. Replaced with matlabs own function,
% juliandate and datetime lines approx 100-120
%
%
%

function bottle = read_rosfile(ros_file) 

 vartime=[];
  dummy = 9999;
  timeS = 'timeS';
  timeJ = 'timeJ';
  timedef = 0 ; 
fid1 = fopen(ros_file,'r');

%scan the header for variable names by column
foundtimes = 0;
while 1
    tline = fgetl(fid1);
    if isempty(tline)
        break
    end
    if strcmp(tline(1),'*') || strcmp(tline(1),'#')
        %header line
        ii1 = strfind(tline,'name ');
        if isempty(ii1) & foundtimes<2
            iist = strfind(tline,'start_time');
            iint = strfind(tline,'NMEA UTC (Time)');
            if ~isempty(iist)
                start_time = tline(iist+13:iist+32);
                start_time = juliandate(datetime(start_time,'InputFormat','MMM dd yyyy HH:mm:ss'));
                foundtimes = foundtimes + 1;
            elseif ~isempty(iint)
                nmea_time  = tline(iint+18:iint+37);
                nmea_time   = juliandate(datetime(nmea_time,'InputFormat','MMM dd yyyy HH:mm:ss'));
                foundtimes = foundtimes + 1;
            end
        else
            ii2 = strfind(tline,'=');
            if ~isempty(ii2)
                vno = str2double(tline(ii1+5:ii2-1));
                if isfinite(vno)
                    vno = vno+1;
                    varnames{vno} = tline(ii2+(2:3));
                    if strcmp(varnames{vno},'ti') && ~exist('tdef','var')
                        tdef.var = tline(ii2+(2:6));
                        tdef.col = vno;
                    end
                end
            end
        end
    else
        %first data line; store then break in order to use textscan
        dat1 = cellfun(@(x) str2double(x), strsplit(strip(tline,"left")));       
        break
    end
end
%read the rest of the lines
a = textscan(fid1,'%f');
a = reshape(a{1},length(varnames),length(a{1})/length(varnames))';
%prepend the first one
a = [dat1; a];

if ~exist('tdef','var') || ~(strcmp(tdef.var,timeS) | strcmp(tdef.var,timeJ))
    disp('Check time variable used in read_rosfile.m')
end

if abs(start_time-nmea_time) >10 % When pb in NMEA time, (e.g. cruise PE399)
    start_time = nmea_time;
    gtime = gregorian(start_time);
    t00 = julian(gtime(1),gtime(2),gtime(3));
    jd         = t00 + (a(:,tdef.col)-1);
else
    if strcmp(tdef.var,timeJ) %time format : timeJ
        gtime= gregorian(start_time);
        jd = julian(gtime(1),1,1) + a(:,tdef.col) - 1;
    else
        jd         = start_time + a(:,tdef.col)/86400;
    end
end

varnames = ['pos','jd',varnames];
a = [(1:length(jd))' jd(:) a];
vars_keep = {'pos','jd','pr','sa','t0','c0','t1','c1'};
vars_new = {'pos','jd','p','s','t0','c0','t1','c1'};
[~,~,ib] = intersect(vars_keep,varnames,'stable');
bottle = table2struct(array2table(a(:,ib),'VariableNames',vars_new),'ToScalar',true);
bottle.start_time = start_time;
