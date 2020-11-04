% Function: rcm2rodb_z('infile','outfile','infofile',fidlog,vec,[toffset])
% By: Hao ZUO (adapted from rcm2rodb.m by Aazani) 18/03/2006

% Reads ACSII output from Aanderaa RCM 11 (310) and converts to RODB format

% Input format contain (YY, MM, DD, HH, MM, SS, reference, current speed, current direction, temperature, conductivity, pressure, tilt, signal strength)
% Output format contain (YY, MM, DD, time decimel hr, u speeds, v speeds, temperature, conductivity, pressure, tilt, signal strength)  
%  
% input: 
%   outfile  : path and name of RODB outputfile
%   infile   : path and name of RCM11 ASCII input file
%   infofile : path and name of mooring info.dat file  
%   fidlog   : file identifier for log file
%   vec      : serial number of the instrument
%   toffset  : offset of recorded instrument time rel to GMT (decimal days),
%              if omitted, default toffset = 0 is set 
%
% opened with: rcm2roodb_scuzuo.m (by Hao Zuo)
% uses: hms2h.m and rodbsave.m (by T.Kanzow), sw_c3515 and sw_salt (CSIRO Sea routines)
% 

%--------------------------------------------------------------------------
% By: HaoZuo
% 18/03/2006 work on the main program
% 19/03/2006 data format transform
% 20/03/2006 run the program with a test input file
% 21/03/2006 basic functions of the program work now with the test file
% 22/03/2006 check for the bad records in the infile
% 23/03/2006 fill the missing time record of 00:00:00
% 24/03/2006 data format transform done and toffset check done

% 27/4/06 DR modified code so that conductivity plot is in mS/cm not S/m
%            modified rodb head variables to use TLT instead of T, MSS
%            instead of S, and REF instead of R.
%            Fixed bug whereby conductivity value saved in .raw file was
%            not the corrected value for the user selected range.
%            Fixed error in converting conductivity using selected ranges -
%            was using 0-100 as base range but should have been 0-75.


function rcm2rodb_z(infile,outfile,infofile,fidlog,vec,c_low,c_upp,toffset)
 
 if nargin < 5
   disp('not enough input arguments')
 end
 
 if nargin < 6
   toffset = 0;
 end

%----------------------------------------------------------
% specifying instr/column and data format info for headers
%----------------------------------------------------------

Instrument   = 'Aanderaa RCM11 (310)'; % Instrument info for rodb header
colheader    = 'YY:MM:DD:HH:REF:U:V:T:C:P:TLT:MSS'; % column info for rodb header 
                % Year, Month, Day, Decimal Hour, Reference number, U & V speeds, Temperature, Conductivity, Pressure, Tilt and Signal Strength)
formatheader = '%4.4d  %2.2d  %2.2d  %7.5f  %3.3d %7.4f  %7.4f  %6.4f  %6.4f  %5.1f %7.6f %6.4f'; %data output format


%----------------------------------------------------------
% check if infile and infofile exist
%----------------------------------------------------------

disp(['Checking if your infile and infofile exist'])

if exist(infofile) ~= 2
   disp(['infofile:  ', infofile,' does not exist'])
   return
end

if exist(infile) ~= 2
   disp(['infile:  ', infile,' does not exist'])
   return
end

if exist(outfile) == 2
   disp(['outfile:  ',outfile,' Already exists!!'])
     overwrite =  input('Overwrite y or n = ','s');
    if overwrite ~='y'
     disp('data conversion stop')
     return
   end
end



%--------------------------------------------------------
% open input file, read data into string
%-------------------------------------------------------- 

SerialNumber=[]; Start_Time=[]; Start_Date=[];

fprintf(fidlog,'infile: %s\n',infile);
fprintf(fidlog,'outfile: %s\n',outfile);

fid = fopen(infile,'r');

if fid == -1
   disp(['unable to open infile:  ',infile])
   return
else
  disp(['loading ',infile]);
end 

zeile = fscanf(fid,'%c');  %read data into string
fclose(fid); %close file 

ret = sprintf('\n');
retx = findstr(zeile,ret);  % car. return indices

if length(retx)< 1
  disp('input file does not contain data ')
  fprintf(fidlog,'input file does not contain data \n')  
  return
end

%---------------------------------------------------------
% find serial number inheader
%-------------------------------------------------------- 

SerialNumber=vec;
fprintf(fidlog,'serial number: %d\n',SerialNumber);


%---------------------------------------------------------------
% get missing header variables from info.dat file
%---------------------------------------------------------------

infovar ='Mooring:Latitude:Longitude:Waterdepth:id:sn:z:StartDate:StartTime:EndDate:EndTime';
[mo,la,lo,wd,id,sn,z,sdate,stime,edate,etime]=rodbload(infofile,infovar);

if isempty(id) | isnan(id)
infovar ='Mooring:Latitude:Longitude:Waterdepth:instrument:SerialNumber:z:StartDate:StartTime:EndDate:EndTime';
[mo,la,lo,wd,id,sn,z,sdate,stime,edate,etime]=rodbload(infofile,infovar);
end

if iscell(mo)
mo = deal(mo{:}); % convert cell array
end

ii = find(SerialNumber == sn);
z  = z(ii);         % instrument depth


%-------------------------------------------
% get data 
%-------------------------------------------

% ---- fill in the time record of 00:00:00 into the gap of the data

data_length=length(retx); % data column length

fprintf(1,'Records number = %d \n', data_length);

for i=1:data_length-1
    row_index(i)=retx(i)+1;  % from the second row
end

temp_outfile=[outfile,'_T'];     % Temp outfile for a recovered '00:00:00' time records input file
if exist(temp_outfile) == 2
   delete(temp_outfile);
end
fid = fopen(temp_outfile,'w+');
disp(['Create new ',temp_outfile]);

string=(zeile(1:retx(1)));

if isempty(findstr(string,':'))   % find the time record missing in the first raw
    slash=findstr(string,'/');
    string=[string(1:slash(2)+5) '00:00:00' string(slash(2)+6:end)];
end

fprintf(fid,'%s',string);

for i=1:data_length-1
    string=zeile(row_index(i):retx(i+1));
    if isempty(findstr(string,':'))   % find the time record missing in the rest raws
        slash=findstr(string,'/');
        string=[string(1:slash(2)+5) ' 00:00:00 ' string(slash(2)+6:end)];
    end
    fprintf(fid,'%s',string);
end

frewind(fid);       % move the file indicator to the beginning

zeile_t = fscanf(fid,'%c');  % reread data into string
fclose(fid);                 % close file 

% -----------define data stream-------------------

dt = zeile_t;

ii = findstr(dt,'\t'); % replace tab by space
dt(ii) = ' '; 
ii = findstr(dt,'/'); % replace / by space
dt(ii) = ' ';  
ii= findstr(dt,':');  % replace : by space
dt(ii)= ' ';

ii = findstr(dt,' .'); % temp. check for missing number

if ~isempty(ii)
     disp('insert dummies for missing elements')
     disp(['missing element element ',num2str(ii+1)])
       fprintf(fidlog,'warning: missing element %s \n',num2str(ii+1));
     for i = 1 : length(ii),
	dt(ii(i):ii(i)+1) = '99';   % All missing data number will be marked as 99
     end
     
end

% --- convert data string-----------

disp('Convert datastr into numbers...')
dt = str2num(dt);

sz = size(dt);
if isempty(dt)
  msgbox(['MC',sprintf('%3.3d',SerialNumber),': deviation from element number'],'conversion stopped')
  fprintf(fidlog,'conversion stopped - deviation from element number \n');
  return
end

dt1=sz(1);              % number of raws
dt2=sz(2);              % number of columns

fprintf(1,'number of rows =%d \n number of columns =%d  \n',dt1,dt2);

% ---- determine time -------

if dt2 == 15
  HH         = hms2h(dt(:,5),dt(:,6),dt(:,7)); % decimal hour
  Start_Date = dt(1,[4 3 2]); % Y, M, D
  End_Date   = dt(dt1,[4 3 2]);
  jd         = julian([dt(:,[4 3 2]) HH]) - toffset;
else
  fprintf(1,'\n WARNING: input file format wrong, data transform stoped');
  return
end
 

Start_Time   = HH(1);
End_Time     = HH(dt1);

bottomstart     = julian([sdate(:)' hms2h([ stime(:)' 0])]);
bottomstop      = julian([edate(:)' hms2h([ etime(:)' 0])]);

if jd(end) < bottomstop
  fprintf(1,'\n\n W A R N I N G: Record already ends before end of mooring deployment period: \n Recorded time could be wrong!! \n')
end 
if jd(1) > bottomstart
  fprintf(1,'\n\n W A R N I N G: Record only starts after beginning of mooring deployment period: \n Recorded time could be wrong!! \n')
end 


% ------converting Speed and Directions to U and V components

% Formula for U and V, from speed, S  and direction, D
% Dr = (90-D) * (pi/180) (% converted to radian)
% U =  S *cos (Dr)
% V =  S *sin (Dr)

disp(['Calculating U and V speeds....'])

for i=1:dt1;
    U(i,1) = (dt(i,9)) *cos((90-dt(i,10))* (pi/180));
    V(i,1)= (dt(i,9)) *sin((90-dt(i,10))* (pi/180));
end

% ---------converting conductivity to real value --------
if c_low == 0 & c_upp == 74 
   % If range is 0-74 then it is likely to be an old 3619 conductivity
   % sensor that has been used so no need to convert scale. This eliminates
   % the problem of having old sensors calibrated on a 0-74 scale and the
   % new 3919 sensors calibrated on a 0-75 scale. DR 27/4/06. 
   for i=1:dt1;
       cond(i,1) = dt(i,12);
   end
else
    for i=1:dt1;
        cond(i,1) = c_low+(c_upp-c_low)*(dt(i,12))/75;
        %cond(i,1) = c_low+(c_upp-c_low)*(dt(i,12))/100;  
    end
end

% ---- Graphics -------------------------


figure(33);

 dl = length(jd);
 
  subplot(6,1,1) % temperature
  hold off
  plot(jd-jd(1),dt(:,11))
     set(gca,'Ygrid','on')
     set(gca,'Xgrid','on')
  ylabel('Temp. (deg C)')
  xlabel('Time [Days]')
  mt = median(dt(:,11));
  fprintf(fidlog,'Median Temperature %5.2f\n',mt);
%   val = find(abs(dt(:,11))<40);
%  sdt = std(dt(val(100:end-100),11)); 
%  set(gca,'ylim',[mt-2.2*sdt mt+2.2*sdt ])
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  yl = get(gca,'ylim');
  xl = get(gca,'xlim');
  tx =text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['mean: ',sprintf('%3.2f',mt)]) ;
  set(tx,'FontWeight','bold','FontSize',11)
  
  subplot(6,1,2)  % conductivity
  hold off
  %plot(jd-jd(1),cond/10)    % input data units is mS/cm
  plot(jd-jd(1),cond)
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  %ylabel('Cond. [S/m]')
  ylabel('Cond. (mS/cm)')
  xlabel('Time [Days]')
  %mt = median(cond(:,1)/10); 
  mt = median(cond(:,1)); 
  fprintf(fidlog,'Median Conductivity %5.2f\n',mt);
%  val = find(abs(dt(:,12))<7);
%  sdt = std(dt(val(100:end-100),12)); 
%  set(gca,'ylim',[mt-2.2*sdt mt+2.2*sdt ])
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  yl = get(gca,'ylim');
  xl = get(gca,'xlim');
  tx=text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['mean: ',sprintf('%3.2f',mt)]); 
  set(tx,'FontWeight','bold','FontSize',11)

  subplot(6,1,3)   % Pressure
  hold off
  plot(jd-jd(1),dt(:,13)*100)       % units transform from Mpa to db
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  ylabel('Pressure (db)')
  xlabel('Time [Days]')
  mt = median(dt(:,13)*100);
  fprintf(fidlog,'Median Pressure %5.1f\n',mt);
%  val = find(abs(dt(:,13))<9000);
%  sdt= std(dt(val(100:end-100),13)); 
%  try
%    set(gca,'ylim',[mt-2.5*sdt mt+2.5*sdt])
%  catch
%   1;
%  end
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  yl = get(gca,'ylim');
  xl = get(gca,'xlim');
  tx=text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['mean: ',sprintf('%3.2f',mt)]); 
  set(tx,'FontWeight','bold','FontSize',11)

  subplot(6,1,4) % compute salinity from measured pressures
  hold off
 % sal=salin78(dt(:,13),dt(:,11),dt(:,12)*10,42.9140,0);
 % SCU JYM 20/4/05 Changed to CSIRO seawater routines
  c3515 = sw_c3515;
  c_ratio = cond(:,1)/c3515;
%  sal = sw_salt(c_ratio,dt(:,11),z*ones(length(dt),1));
  sal = sw_salt(c_ratio,dt(:,11),dt(:,13)*100);
 % end SCU/JYM
 
  plot(jd-jd(1),sal)
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  ylabel('Sal.')
  xlabel('Time [Days]')
  mt = median(sal);
  fprintf(fidlog,'Median Salinity %5.2f\n',mt);
%  val = find(abs(sal)<70);
%  sdt = std(sal(val(100:end-100))); 
%  try 
%    set(gca,'ylim',[mt-1.8*sdt mt+1.8*sdt])
%  catch 
%    1;
%  end 
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  yl = get(gca,'ylim');
  xl = get(gca,'xlim');
  tx=text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['mean: ',sprintf('%3.2f',mt)]); 
  set(tx,'FontWeight','bold','FontSize',11)
  
subplot (6,1,5) % U
hold off
plot (jd-jd(1), U)
    set (gca, 'Ygrid', 'on')
    set (gca, 'Xgrid', 'on')
ylabel ('U (cm/s)')
xlabel ('Time [Days]')
    meanU = mean (U)
fprintf (fidlog, 'Mean U Component = %7.4f\n', meanU);
%val = find(abs(U)<200);
%sdt = std(U(val(100:end-100))); 
%set(gca,'ylim',[mt-2.2*sdt mt+2.2*sdt ])
set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))

subplot (6,1,6) % V
hold off
plot (jd-jd(1), V)
    set (gca, 'Ygrid', 'on')
    set (gca, 'Xgrid', 'on')
ylabel ('V (cm/s)')
xlabel ('Time [Days]')
    meanV = mean (V)
fprintf (fidlog, 'Mean V Component = %7.4f\n', meanV);
%val = find(abs(V)<200);
%sdt = std(V(val(100:end-100))); 
%set(gca,'ylim',[mt-2.2*sdt mt+2.2*sdt ])
set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))

subplot(6,1,1)

title(['RCM11:  ',num2str(SerialNumber),'   Depth:  ',num2str(z),' m']) 
orient tall
eval(['print -dps ',outfile,'.ps'])


%-------------------------------------
% ---- save to rodb format
%-------------------------------------

disp(['writing data to ',outfile]) 

TIME = gregorian(jd);

%data = [TIME(:,1:3) hms2h(TIME(:,4:6)) dt(:,8) U V dt(:,11:12) dt(:,13)*100 dt(:,14:15)]; 
data = [TIME(:,1:3) hms2h(TIME(:,4:6)) dt(:,8) U V dt(:,11) cond dt(:,13)*100 dt(:,14:15)]; 

rodbsave(outfile,...
       'Latitude:Longitude:Columns:SerialNumber:Mooring:WaterDepth:Instrdepth:StartDate:StartTime:EndDate:EndTime',...
         formatheader,...
         la,lo,colheader,SerialNumber,mo,wd,z,...
         sdate,stime,edate,etime,...
         data);

fprintf(fidlog,'Instrument Target Depth[m]: %d\n',z);
fprintf(fidlog,'Start date and time: %s \n',datestr(gregorian(jd(1))));
fprintf(fidlog,'End date and time:   %s \n',datestr(gregorian(jd(end))));
sampling_rate = round(1./median(diff(jd)));
ex_samples = round((jd(end)-jd(1))*sampling_rate+1);
fprintf(fidlog,'Sampling Frequency [per jd]: %d \n',sampling_rate);
fprintf(fidlog,'Number of samples: %d; expected: %d \n',sz(1),ex_samples);

if toffset ~= 0
  fprintf(fidlog,'Offset of %8.4f days has been subtracted from recored time \n',toffset);
  fprintf(1,'Offset of %8.4f days has been subtracted from recored time \n',toffset);
end

