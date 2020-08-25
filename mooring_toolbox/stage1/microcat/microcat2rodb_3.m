% function microcat2rodb_3('infile','outfile','infofile',fidlog,[graphics],[toffset])
%
% reads ACSII output from SBE-37 MicroCAT and converts
% writes it to RODB file 
% So far the input format contain (temp,cond,day,month,year,time)  
% or .cnv format (temp,cond,press,seconds,data_flag)
% where seconds is seconds since Jan 1st 2000.
%
% input: 
%   outfile  : path and name of RODB outputfile
%   infile   : path and name of MicroCAT ASCII input file
%   infofile : path and name of mooring info.dat file  
%   fidlog   : file identifier for log file
%   graphics : 'y' = display some graphics (mooring)
%              'w' = display some graphics with whole range (rosette)
%   toffset  : offset of recorded instrument time rel to GMT (decimal days),
%              if omitted, default toffset = 0 is set 
%
%  uses: microcat_month2.m, hms2h.m  (all by T.Kanzow), rodbsave.m
%
% kanzow   26.12.2000 Xmas edition
%          28.03.2001 pressure option added 
%          20.04.2005 CSIRO seawater routine for computation of salinity 
%          24.04.2005 recorded time checked against mooring deployment time 
%          25.04.2005 input variable 'toffset' added
%          27/10/2009 - DR added functionality for .cnv format files for
%             microcat firmware 3.0 and above. and created microcat2rodb_3
%             from microcat2rodb_2_002
%          03.04.2010 - ZBS changed the method for setting 'ylim'
%             for graphics=='y' to use prctile.
%          20/11/2012 - DR added lines to tell which sample number has
%             format problems following issues with bad scans in the middle
%             of records on D382
%          30/04/2014 - DR added compatability for MicroCAT ODOs
%
% bim      2-May-2014 check added to make sure microcat is selected using id and serial number
%                     e.g. nortek and microcat share serial numbers (6805)

function microcat2rodb_3(infile,outfile,infofile,fidlog,graphics,toffset);

if nargin < 4
  disp('not enough input arguments')
elseif nargin == 4
  graphics = 'n'
  disp('graphics will not be  displayed')
end

if nargin < 5
  toffset = 0;
end

Instrument = 'MicroCAT';           % Instr. info for rodb header
cols       = 'YY:MM:DD:HH:T:C';    % column info for rodb header
colsp      = 'YY:MM:DD:HH:T:C:P'; % column info for rodb header (mc with pressure sensor)
colso      = 'YY:MM:DD:HH:T:C:P:OT:O2'; % column info for rodb header (mc with pressure sensor and oxygen sensor)
fort       = '%4.4d  %2.2d  %2.2d  %7.5f   %6.4f  %6.4f'; %data output format
fortp      = '%4.4d  %2.2d  %2.2d  %7.5f   %6.4f  %6.4f  %5.1f'; %data output format(mc with pressure sensor)  
forto      = '%4.4d  %2.2d  %2.2d  %7.5f   %6.4f  %6.4f  %5.1f %6.4f %5.2f'; %data output format(mc with pressure sensor and oxygen sensor)  

% check if infile and infofile exist

if exist(infofile) ~= 2
   disp(['infofile:  ',infofile,' does not exist'])
   pause
   return
end

if exist(infile) ~= 2
   disp(['infile:  ',infile,' does not exist'])
   pause
   return
end

if exist(outfile) == 2
   disp(['oufile:  ',outfile,' alredy exists!!'])
     overwrite =  input('Overwrite y/n  ','s');
    if overwrite ~='y'
     disp('data conversion stop')
     return
   end
end

%--------------------------------------------------------
% check for .cnv file extension
%--------------------------------------------------------
if (strcmp(infile(end-2:end),'cnv') | strcmp(infile(end-2:end),'CNV'))
    cnv=1;
else
    cnv=0;
end

%---------------------------------------------------------
% open input file, read data into string
%-------------------------------------------------------- 

SerialNumber=[]; Start_Time=[]; Start_Date=[];

fprintf(fidlog,'infile: %s\n',infile);
fprintf(fidlog,'outfile: %s\n',outfile);
fid1 = fopen(infile,'r');

if fid1 == -1
   disp(['unable to open infile:  ',infile])
   return
else
  disp(['loading ',infile]);
end 

zeile = fscanf(fid1,'%c');  %read data into string
fclose(fid1); %close file 

ret = sprintf('\n');
retx = findstr(zeile,ret);% car. return indices

if length(retx)< 10
  disp('input file does not contain data ')
  fprintf(fidlog,'input file does not contain data \n')  ;
  return
end


%---------------------------------------------------------
% find serial number inheader
%-------------------------------------------------------- 

if cnv==0
    a = findstr(zeile,'SERIAL NO.'); %find serial number
    if ~isempty(a)
     SerialNumber = str2num(zeile(a+10:a+15))  
    else
        disp('unable to find serial number\n\n')  
        input('please enter valid serial number')
        SerialNumber = str2num(input('please enter valid serial number'));
    end
else
    a=findstr(zeile,'Temperature SN'); % find serial number for .cnv files
    % currently uses Temperature sensor serial number as this is the same
    % as the instrument serial number and is unique in the header info but 
    % may run into problems if temperature serial number does not match 
    % instrument sn.
    if ~isempty(a)
        SerialNumber = str2num(zeile(a+17:a+20));  
    else
        disp('unable to find serial number\n\n')  
        input('please enter valid serial number')
        SerialNumber = str2num(input('please enter valid serial number'));
    end
end


fprintf(fidlog,'serial number: %d\n',SerialNumber);

%---------------------------------------------------------------
% get missing header variables from info.dat file
%---------------------------------------------------------------

infovar ='Mooring:Latitude:Longitude:Waterdepth:id:sn:z:StartDate:StartTime:EndDate:EndTime'; 

[mo,la,lo,wd,id,sn,z,sdate,stime,edate,etime]=rodbload(infofile,infovar); 
if isempty(id) | isnan(id)
  infovar ='Mooring:Latitude:Longitude:Waterdepth:instrument:serialnumber:z:StartDate:StartTime:EndDate:EndTime'; 
  [mo,la,lo,wd,id,sn,z,sdate,stime,edate,etime]=rodbload(infofile,infovar); 
end

if iscell(mo)
  mo = deal(mo{:}); % convert cell array
end

ii = find(SerialNumber == sn & ( id >= 333 & id <= 337 ) );
z  = z(ii)         % instrument depth


%-------------------------------------------
% get data 
%-------------------------------------------

% detect data column length

[XXX,data_length] = max(hist(diff(retx),1:300));% length of data columns 
% DR increased final term to 200 from 100 due to longer header in .cnv
% files   % changed to 300 PW

ii0 = find(diff(retx) == data_length); % data row index  
ii1 = find(diff(ii0)>1);  % 

 bg = 2;
[XXX,sede] = find(abs(diff(retx(ii0(bg):ii0(length(ii0))))-data_length) >3);

if ~isempty(XXX)
  bg = bg+1;
  [XXX,sede] = find(abs(diff(retx(ii0(bg):ii0(length(ii0))))-data_length) >3);
end

if ~isempty(XXX)
  bg = bg +1;
  [XXX,sede] = find(abs(diff(retx(ii0(4):ii0(length(ii0))))-data_length) >3);
end


if ~isempty(XXX)
  disp('severe deviation from input format')
  fprintf(fidlog,'coversion stop - severe deviation from input format: %d\n',sede+ii0(1));
  msgbox(['MC',sprintf('%4.4d','SerialNumber'),'  ',sprintf('severe deviation from input format: %d\n',sede+ii0(1))],'conversion stopped')  
 return
end

% check if there are deviations in data row length


  data_begin = retx(ii0(bg))+1;
  data_end   = retx(ii0(length(ii0)));
  disp('data begin detected')   

  if ~isempty(ii1)
    disp('warning: deviation from format')
  fprintf(fidlog,'warning: deviation from format \n');
  end



% define data stream

dt = zeile(data_begin:data_end);
dret = find(dt == ret);
comx = findstr(dt(1:dret(1)),',');
retn = length(dret);% number of data columns

ii = findstr(dt,','); % replace comma by space
dt(ii) = ' '; 
ii = findstr(dt,':'); % replace colon by space
dt(ii) = ' '; 

if cnv==0 %microcat_month2 only valid for .asc format files
    dt = microcat_month2(dt);
end

whos dt;
ii = findstr(dt,' .'); % temp. check for missing number

if ~isempty(ii)
     disp('insert dummies for missing elements')
     disp(['missing element element ',num2str(ii+1)])
       fprintf(fidlog,'warning: missing element %s \n',num2str(ii+1));
     for i = 1 : length(ii);
	dt(ii(i):ii(i)+1) = '99';
     end
     
end 


% --- convert data string

disp('Convert datastr into numbers')
%save tmp.mat dt ;


% This next line won't work if there're still ??? (carriage returns are ok) in dt
dt = str2num(dt);





% unique fix for d345 cal dip files for MCs 7347 and 7348 which have C in
% wrong format
if (strfind(infile,'d345') & strfind(infile,'7347_cal_dip'));
      check1=1;
else
      check1=0;
end
if check1|(strfind(infile,'d345') & strfind(infile,'7348_cal_dip'));
    dt(:,2)=dt(:,2)/10;
end

sz = size(dt);

% catch and locate errors in raw microcat files where have spurious data
% scans for whatever reason. Typically these come at the end of the record,
% but on D382 we found some bad data in the middle of a record that was
% difficult to track down
if isempty(dt);
    load tmp.mat % reload dt data
    k2=1;
    for k=1:length(dt)/57
        string=dt(57*(k-1)+1:57*(k-1)+55);
        dt_convert=str2num(string);
        if isempty(dt_convert)
            bad_data_line(k2)=k;
            k2=k2+1;
         end
    end
  
  msgbox(['MC',sprintf('%4.4d',SerialNumber),': format deviation at sample ' num2str(bad_data_line(1))],'conversion stopped');
  fprintf(fidlog,'conversion stopped - deviation from element number \n');
  
   return
end



% ---- determine time -------

dtl        = size(dt,1);

if cnv==0 % data is in .asc file
    if size(dt,2) ==8

      HH         = hms2h(dt(:,6),dt(:,7),dt(:,8)); % decimal hour
      Start_Date = dt(1,[5 4 3]);
      End_Date   = dt(dtl,[5 4 3]);
      jd         = julian([dt(:,[5 4 3]) HH]) - toffset;

    elseif  size(dt,2) ==9;

      HH         = hms2h(dt(:,7),dt(:,8),dt(:,9)); % decimal hour
      Start_Date = dt(1,[6 5 4]);
      End_Date   = dt(dtl,[6 5 4]);
      jd         = julian([dt(:,[6 5 4]) HH]) - toffset;

    end
else % data is in cnv file
    if size(dt,2)==5 % SMP with pressure
        cnv_secs=dt(:,4); % seconds since 1st Jan 2000.
        jd=cnv_secs/(60*60*24)+julian(2000,1,1,0) - toffset;
        disp(['toffset : ',num2str(toffset)]);
        gtime=gregorian(jd);
        HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6));
        Start_Date = gtime(1,[1 2 3]);
        End_Date = gtime(dtl,[1 2 3]);
    elseif size(dt,2)==7 % SMP-ODO with pressure
        cnv_secs=dt(:,6); % seconds since 1st Jan 2000.
        jd=cnv_secs/(60*60*24)+julian(2000,1,1,0) - toffset;
        disp(['toffset : ',num2str(toffset)]);
        gtime=gregorian(jd);
        HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6));
        Start_Date = gtime(1,[1 2 3]);
        End_Date = gtime(dtl,[1 2 3]);
    end
end

Start_Time   = HH(1);
End_Time     = HH(dtl);

bottomstart     = julian([sdate(:)' hms2h([ stime(:)' 0])]);
bottomstop      = julian([edate(:)' hms2h([ etime(:)' 0])]);

if jd(end) < bottomstop
  fprintf(1,'\n\n W A R N I N G: Record already ends before end of mooring deployment period: \n Recorded time could be wrong!! \n');
end 
if jd(1) > bottomstart
  fprintf(1,'\n\n W A R N I N G: Record only starts after beginning of mooring deployment period: \n Recorded time could be wrong!! \n');
end 

% ---- Graphics ------------------ -------


% plot data if graphics == 'y'

if graphics == 'y'  | graphics == 'w'
 
figure(33)

if size(dt,2) ==8

  dl = length(jd);
  

  subplot(3,1,1) % temperature
  hold off
  plot(jd-jd(1),dt(:,1))
     set(gca,'Ygrid','on')
     set(gca,'Xgrid','on')
  ylabel('Temp.')
   mt = median(dt(:,1));
  fprintf(fidlog,'Median Temperature %5.2f\n',mt);
   val = find(abs(dt(:,1))<40);
  sdt = std(dt(val(100:end-100),1)); 
  %set(gca,'ylim',[mt-2.2*sdt mt+2.2*sdt ])
  set(gca,'ylim',prctile(dt(:,1),[1 99])+[-1 1]*std(dt(:,1))/10)
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  if  graphics == 'w'
    set(gca,'ylim',[min(dt(:,1))  min([max(dt(:,1)) 30])])
  end
  yl = get(gca,'ylim');
  xl = get(gca,'xlim');
  tx =text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['mean: ',sprintf('%3.2f',mt)]) ;
  set(tx,'FontWeight','bold','FontSize',11)

  subplot(3,1,2)  % conductivity
  hold off
  
  
  plot(jd-jd(1),dt(:,2))
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  ylabel('Cond. [S/m]')
  xlabel('Time [Days]')
  mt = median(dt(:,2)); 
  fprintf(fidlog,'Median Conductivity %5.2f\n',mt);

  val = find(abs(dt(:,2))<7);
  sdt = std(dt(val(100:end-100),2)); 
   
  %set(gca,'ylim',[mt-2.2*sdt mt+2.2*sdt ])
  set(gca,'ylim',prctile(dt(:,2),[1 99])+[-1 1]*std(dt(:,2))/10)
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  if  graphics == 'w'
    set(gca,'ylim',[max([min(dt(:,2)) 2])  max(dt(:,2))])
  end
  yl = get(gca,'ylim');
  xl = get(gca,'xlim');

  tx=text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['mean: ',sprintf('%3.2f',mt)]); 
  set(tx,'FontWeight','bold','FontSize',11)

  subplot(3,1,3) % compute salinity from aproximate depth
  hold off
  %sal=salin78(z*ones(length(dt),1),dt(:,1),dt(:,2)*10,42.9140,0);

 % SCU JYM 20/4/05 Changed to CSIRO seawater routines
  c3515 = sw_c3515;
  
  c_ratio = (dt(:,2)*10)/c3515;
  sal = sw_salt(c_ratio,dt(:,1),z*ones(length(dt),1));
 % end SCU/JYM
  
  plot(jd-jd(1),sal)
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  ylabel('Sal.')
  xlabel('Time [Days]')
  mt = median(sal);
  val = find(abs(sal)<70& abs(sal)>30);
  sdt = std(sal(val(100:end-100))); 
  %set(gca,'ylim',[mt-3*sdt mt+3*sdt ])
  set(gca,'ylim',prctile(sal,[1 99])+[-1 1]*std(sal)/10)
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  if  graphics == 'w'
    set(gca,'ylim',[max([min(dt(:,2)) 2])  max(dt(:,2))])
  end
  yl = get(gca,'ylim');
  xl = get(gca,'xlim');

  tx=text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['mean: ',sprintf('%3.2f',mt)]); 
  set(tx,'FontWeight','bold','FontSize',11)

  subplot(3,1,1)

elseif (size(dt,2) ==9 | (size(dt,2)==5 & cnv==1))
  dl = length(jd);
  
  subplot(4,1,1)
  hold off
  plot(jd-jd(1),dt(:,1))
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  ylabel('Temp.')
  mt = median(dt(:,1));
  fprintf(fidlog,'Median Temperature %5.2f\n',mt);

  
  val = find(abs(dt(:,1))<40);

  sdt = std(dt(val(100:end-100),1)); 
  try
    %set(gca,'ylim',[mt-1.8*sdt mt+1.8*sdt ])
    set(gca,'ylim',prctile(dt(:,1),[1 99])+[-1 1]*std(dt(:,1))/10)
  catch
  1;
  end
  try
     set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  catch
    1;
  end
  if  graphics == 'w'
    try
       set(gca,'ylim',[min(dt(:,1))  min([max(dt(:,1)) 30])])
    catch
       1;
    end 
  end
 
  subplot(4,1,2)
  hold off
  
  plot(jd-jd(1),dt(:,2))
  
     set(gca,'Ygrid','on')
     set(gca,'Xgrid','on')
  ylabel('Cond. [S/m]')
  xlabel('Time [Days]')
  mt = median(dt(:,2));
  fprintf(fidlog,'Median Conductivity %5.2f\n',mt);

  val = find(abs(dt(:,2))<7);

  sdt = std(dt(val(100:end-100),2)); 
 try 
   %set(gca,'ylim',[mt-1.8*sdt mt+1.8*sdt ])
   set(gca,'ylim',prctile(dt(:,2),[1 99])+[-1 1]*std(dt(:,2))/10)
 catch
  1;
 end 
%  set(gca,'ylim',[mt-.03 mt+.03])
  try
    set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  end
  if  graphics == 'w'
    try
      set(gca,'ylim',[max([min(dt(:,2)) 2])  max(dt(:,2))])
    catch
      set(gca,'ylim',[min(dt(:,2))   max(dt(:,2))])
    end   
  end

  subplot(4,1,3)
  hold off
  plot(jd-jd(1),dt(:,3))
     set(gca,'Ygrid','on')
     set(gca,'Xgrid','on')
  ylabel('Pressure')
  xlabel('Time [Days]')
  mt = median(dt(:,3));
  fprintf(fidlog,'Median Pressure %5.1f\n',mt);

  val = find(abs(dt(:,3))<9000);

  sdt= std(dt(val(100:end-100),3)); 
  try
    %set(gca,'ylim',[mt-2.5*sdt mt+2.5*sdt])
    set(gca,'ylim',prctile(dt(:,3),[1 99])+[-1 1]*std(dt(:,3))/5)
  catch
   1;
  end
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))

  if  graphics == 'w'
% SCU
disp([max([min(dt(:,3)) 0])  max(dt(:,3))])
      set(gca,'ylim',[max([min(dt(:,3)) 0])  max(dt(:,3))])
  end

  subplot(4,1,4) % compute salinity from measured pressures
  hold off

 % sal=salin78(dt(:,3),dt(:,1),dt(:,2)*10,42.9140,0);
 % SCU JYM 20/4/05 Changed to CSIRO seawater routines
  c3515 = sw_c3515;
  
  % unique fix for MCs 7347 and 7348 from D345 cal dips - these files were
  % downloaded by RSMAS and given to us in the wrong format for C. Should
  % be in S/m in the file from the instrument but for some reason these are
  % in mS/cm
  
  c_ratio = (dt(:,2)*10)/c3515;
  
  sal = sw_salt(c_ratio,dt(:,1),dt(:,3));
 % end SCU/JYM
 
  plot(jd-jd(1),sal)
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  ylabel('Sal.')
  xlabel('Time [Days]')
  mt = median(sal);
  fprintf(fidlog,'Median Salinity %5.2f\n',mt);

  val = find(abs(sal)<70);
  sdt = std(sal(val(100:end-100))); 
  try 
    %set(gca,'ylim',[mt-1.8*sdt mt+1.8*sdt])
    set(gca,'ylim',prctile(sal,[1 99])+[-1 1]*std(sal)/10)
  catch 
    1;
  end 
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  if  graphics == 'w'
    try
      set(gca,'ylim',[max([min(dt(:,2)) 2])  max(dt(:,2))])
    catch
     1;
    end
  end
  yl = get(gca,'ylim');
  xl = get(gca,'xlim');

  tx=text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['mean: ',sprintf('%3.2f',mt)]); 
  set(tx,'FontWeight','bold','FontSize',11)
  
  subplot(4,1,1)
  
elseif (size(dt,2)==7 & cnv==1) % has ODO sensor as well as pressure and .cnv extension
  dl = length(jd);
  
  subplot(6,1,1)
  hold off
  plot(jd-jd(1),dt(:,1))
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  ylabel('Temp.')
  mt = median(dt(:,1));
  fprintf(fidlog,'Median Temperature %5.2f\n',mt);

  val = find(abs(dt(:,1))<40);

  sdt = std(dt(val(100:end-100),1)); 
  try
    set(gca,'ylim',prctile(dt(:,1),[1 99])+[-1 1]*std(dt(:,1))/10)
  catch
  1;
  end
  try
     set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  catch
    1;
  end
  if  graphics == 'w'
    try
       set(gca,'ylim',[min(dt(:,1))  min([max(dt(:,1)) 30])])
    catch
       1;
    end 
  end
 
  subplot(6,1,2)
  hold off
  
  plot(jd-jd(1),dt(:,2))
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  ylabel('Cond. [S/m]')
  xlabel('Time [Days]')
  mt = median(dt(:,2));
  fprintf(fidlog,'Median Conductivity %5.2f\n',mt);

  val = find(abs(dt(:,2))<7);

  sdt = std(dt(val(100:end-100),2)); 
 try 
   %set(gca,'ylim',[mt-1.8*sdt mt+1.8*sdt ])
   set(gca,'ylim',prctile(dt(:,2),[1 99])+[-1 1]*std(dt(:,2))/10)
 catch
  1;
 end 
%  set(gca,'ylim',[mt-.03 mt+.03])
  try
    set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  end
  if  graphics == 'w'
    try
      set(gca,'ylim',[max([min(dt(:,2)) 2])  max(dt(:,2))])
    catch
      set(gca,'ylim',[min(dt(:,2))   max(dt(:,2))])
    end   
  end

  subplot(6,1,3)
  hold off
  plot(jd-jd(1),dt(:,3))
     set(gca,'Ygrid','on')
     set(gca,'Xgrid','on')
  ylabel('Pressure')
  xlabel('Time [Days]')
  mt = median(dt(:,3));
  fprintf(fidlog,'Median Pressure %5.1f\n',mt);

  val = find(abs(dt(:,3))<9000);

  sdt= std(dt(val(100:end-100),3)); 
  try
    %set(gca,'ylim',[mt-2.5*sdt mt+2.5*sdt])
    set(gca,'ylim',prctile(dt(:,3),[1 99])+[-1 1]*std(dt(:,3))/5)
  catch
   1;
  end
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))

  if  graphics == 'w'
% SCU
disp([max([min(dt(:,3)) 0])  max(dt(:,3))])
      set(gca,'ylim',[max([min(dt(:,3)) 0])  max(dt(:,3))])
  end

  subplot(6,1,4) % compute salinity from measured pressures
  hold off

 % sal=salin78(dt(:,3),dt(:,1),dt(:,2)*10,42.9140,0);
 % SCU JYM 20/4/05 Changed to CSIRO seawater routines
  c3515 = sw_c3515;
  
  % unique fix for MCs 7347 and 7348 from D345 cal dips - these files were
  % downloaded by RSMAS and given to us in the wrong format for C. Should
  % be in S/m in the file from the instrument but for some reason these are
  % in mS/cm
  
  c_ratio = (dt(:,2)*10)/c3515;
  
  sal = sw_salt(c_ratio,dt(:,1),dt(:,3));
 % end SCU/JYM
 
  plot(jd-jd(1),sal)
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  ylabel('Sal.')
  xlabel('Time [Days]')
  mt = median(sal);
  fprintf(fidlog,'Median Salinity %5.2f\n',mt);

  val = find(abs(sal)<70);
  sdt = std(sal(val(100:end-100))); 
  try 
    %set(gca,'ylim',[mt-1.8*sdt mt+1.8*sdt])
    set(gca,'ylim',prctile(sal,[1 99])+[-1 1]*std(sal)/10)
  catch 
    1;
  end 
  set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  if  graphics == 'w'
    try
      set(gca,'ylim',[max([min(dt(:,2)) 2])  max(dt(:,2))])
    catch
     1;
    end
  end
  yl = get(gca,'ylim');
  xl = get(gca,'xlim');

  tx=text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['mean: ',sprintf('%3.2f',mt)]); 
  set(tx,'FontWeight','bold','FontSize',11)
  
  subplot(6,1,5)
  hold off
  plot(jd-jd(1),dt(:,4))
  set(gca,'Ygrid','on')
  set(gca,'Xgrid','on')
  ylabel('Oxy. Temp.')
  mt = median(dt(:,4));
  fprintf(fidlog,'Median Oxygen Temperature %5.2f\n',mt);

  val = find(abs(dt(:,4))<40);

  sdt = std(dt(val(100:end-100),1)); 
  try
    set(gca,'ylim',prctile(dt(:,4),[1 99])+[-1 1]*std(dt(:,4))/10)
  catch
  1;
  end
  try
     set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
  catch
    1;
  end
  if  graphics == 'w'
    try
       set(gca,'ylim',[min(dt(:,4))  min([max(dt(:,4)) 30])])
    catch
       1;
    end 
  end
 
    subplot(6,1,6)
    hold off
    plot(jd-jd(1),dt(:,5))
    set(gca,'Ygrid','on')
    set(gca,'Xgrid','on')
    ylabel('Oxy.')
    mo = median(dt(:,5));
    fprintf(fidlog,'Median Oxygen %5.2f\n',mo);

  
    val = find(abs(dt(:,5))<400);

    sdo = std(dt(val(100:end-100),5)); 
  try
    %set(gca,'ylim',[mt-1.8*sdt mt+1.8*sdt ])
    set(gca,'ylim',prctile(dt(:,5),[1 99])+[-1 1]*std(dt(:,5))/10)
  catch
  1;
  end
    try
        set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
    catch
        1;
    end
    if  graphics == 'w'
        try
            set(gca,'ylim',[min(dt(:,5))  min([max(dt(:,5)) 400])])
        catch
            1;
        end 
    end
  
  subplot(6,1,1)
end
  
  title(['Microcat:  ',num2str(SerialNumber),'   Depth:  ',num2str(z),' m']) 
  orient tall
  eval(['print -dpsc ',outfile,'.ps'])
end




%-------------------------------------
% ---- save to rodb format
%-------------------------------------


disp(['writing data to ',outfile]) 

TIME = gregorian(jd);


data = [TIME(:,1:3) hms2h(TIME(:,4:6)) dt(:,1) dt(:,2)*10]; 
  

if cnv==0 %.asc file
    if  size(dt,2) ==9           % with pressure
      data = [data dt(:,3)];
      cols = colsp;
      fort = fortp;   
    end
else % .cnv file
    if size(dt,2)==5
        data = [data dt(:,3)];
        cols = colsp;
        fort = fortp;
    elseif size(dt,2)==7 % ODO sensor
        data = [data dt(:,3:5)];
        cols = colso;
        fort = forto;
    end
end

rodbsave(outfile,...
       'Latitude:Longitude:Columns:SerialNumber:Mooring:WaterDepth:Instrdepth:StartDate:StartTime:EndDate:EndTime',...
         fort,...
         la,lo,cols,SerialNumber,mo,wd,z,...
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
% jym 21 april 2005: alert operator to mismatch in number of data cycles
if sz(1) ~= ex_samples
   fprintf('Number of samples: %d; expected: %d \n',sz(1),ex_samples);
   disp('Press return to continue');
   pause 
end

%tk 24/04/05
if jd(end) < bottomstop
  fprintf(fidlog,'\n\n W A R N I N G: Record already ends before end of mooring deployment period: \n Recorded time could be wrong!! \n');
end 
if jd(1) > bottomstart
  fprintf(fidlog,'\n\n W A R N I N G: Record only starts after beginning of mooring deployment period: \n Recorded time could be wrong!! \n');
end 






