% function rcm2rodb_04(moor, 'procpath', procpath, 'inpath', inpath, ...
%                      'outpath', outpath,'ctd')
%
% Reads ACSII output from Aanderaa RCM 11 (310) and converts to RODB format
%
% Required inputs:
%   moor = mooring name as string e.g. 'wb1_2_200527'
%
% Optional inputs:
%   'procpath' = path to proc directory if not using standard paths e.g.
%               '/Volumes/jrd/jrd/hydro10/rapid/data/moor/proc/'
%               Need to specify procpath if processing ctd cast as
%               proc_calib directory is seperated by cruise not mooring id.
%               e.g. '/local/users/pstar/Data/rpdmoc/rapid/data/moor/proc_calib/rb0602/cal_dip' 
%               is for the ctd casts on rb0602.
%   'inpath'   = path to raw rcm files - if not using standard raw paths
%               standard input path =
%               /local/users/pstar/Data/rpdmoc/rapid/data/moor/raw/"MOORING"/rcm/ 
%               where "MOORING" is the mooring name e.g. wb1_2_200527
%               (once this standard has been adopted which as yet it hasn't)
%   'outpath'  = path for output .raw files - default is directory function
%                run from
%   'ctd'      = include this if processing a ctd cast.
%
% Functions used:
% rodbload.m
% rcm2rodb_z02.m   - NB: This only works with Matlab v14.

%--------------------------------------------------------------------------
% By: Hao Zuo (adapted from rcm2rodb_scu.m by Aazani)
% need function rcm2rodb_z.m (adapted from rcm2rodb.m)
% moor number need to be specified\
% input files should be named with the formated of ***_data.asc and *** is
% the serial number of the instruments.
%
% CHANGES:
% 27/4/06 - create version 02 - DR modified code to prompt for inpath, outpath and moor
%         - Hao's original file is rcm2rodb_scuzuo.m
% 25/7/06 - create version 03 - change to function so that have optional/standard inputs for
%           mooring name and paths. DR.
% 11/10/06 - create version 04 - rcm2rodb_z02.m function now included in this code instead of
%            seperate .m file
% 15/12/06 - correct bug with procpath for cal_dips.
% 7/5/08 - correct pressure conversion for serial numbers > 500 (units kPa instead of MPa)
%        - create version 05 - rcm2rodb_05.m
% edited inpath to read cruise d344 rather than mooring PGW Nov 09


function rcm2rodb_04(moor,varargin)

if nargin==0
    help rcm2rodb_04
    return
end

% check for optional arguments
a=strmatch('procpath',varargin,'exact');
if a>0
    procpath=char(varargin(a+1));
else
    procpath='/local/users/pstar/Data/rpdmoc/rapid/data/moor/proc';
end

a=strmatch('inpath',varargin,'exact');
if a>0
    inpath=char(varargin(a+1));
else
    inpath=eval(['''/local/users/pstar/Data/rpdmoc/rapid/data/moor/raw/d344/rcm11/'';']);
end
disp(['inpath = ' inpath])
a=strmatch('outpath',varargin,'exact');
if a>0
    outpath=char(varargin(a+1));
else
    outpath = ['/local/users/pstar/Data/rpdmoc/rapid/data/moor/proc/' ,moor, '/rcm'];
end


% --- get moring information from infofile 
a=strmatch('ctd',varargin,'exact');
if a>0
    infofile = [procpath '/' moor 'info.dat'];
else
    infofile = [procpath '/' moor '/' moor 'info.dat'];
end
disp(['info file = ' infofile]);




[gash, operator]=system('whoami'); % This line will not work if run from a PC. May need to edit it out.
% % gash is not used, but a second variable needs to be specified for the system command


out_ext  = ['.raw'];



toffset = 0;

% vector of serial numbers

[id,sn,c1,c2]= rodbload(infofile,'id:sn:rcmc1:rcmc2');
if isempty(id) | isnan(id)  ;
[id,sn,lat,lon]= rodbload(infofile,'instrument:serialnumber:latitude:longitude');
end

ii = find(id == 310);
vec = sn(ii);           % serial numbers of RCM11
c_low = c1(ii);         % lower limit of RCM11 conductivity range [mS/cm]
c_upp = c2(ii);         % upper limit of RCM11 conductivity range [mS/cm]

fidlog = fopen([outpath,'/RCM_stage1_log'],'a');
fprintf(fidlog,'Transformation of ascii data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));

fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n',lon);




 for i = 1:length(vec),
    fprintf(fidlog,'\n\n');
 
    if findstr(moor,'cast')
        ctdcast=['_cal_dip'];
    else
        ctdcast=[''];
    end
 infile = [inpath,sprintf('%3.3d',vec(i)),ctdcast,'_data.asc'];
    
    % ......in case of other infile names.....
    % if exist(infile)~=2
    %   ......... 
    % end
    
 outfile = [outpath,'/',moor,'_',sprintf('%3.3d',vec(i)),out_ext];
    
%    infile = ['/local1/rapid/data/moor/raw/WB1_2005/rcm/test_DATA.ASC'];
%    outfile = [outpath,'test.raw'];
%    using a test input file to test the program
    
    rcm2rodb_z02(infile,outfile,infofile,fidlog,vec(i),c_low(i),c_upp(i),toffset)
    
    disp(['proceeding to next file ']);

end

comment = input('Enter additional comment to be save in Log file? ','s'); 

if ~isempty(comment)
  fprintf(fidlog,'\n COMMENT:\n %s',comment)
end

fclose(fidlog);

end


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
%       data_point.m function included in this file for editing wrapped conductivity data.
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

% 11/10/06 DR included section on correcting wrapped conductivity data.
%             also included data_point function at end.


function rcm2rodb_z02(infile,outfile,infofile,fidlog,vec,c_low,c_upp,toffset)
 
 if nargin < 7
   disp('not enough input arguments')
 end
 
 if nargin < 8
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

if isempty(id) || isnan(id)
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
disp(' ')
disp(['Calculating U and V speeds....'])
disp(' ')

for i=1:dt1;
    U(i,1) = (dt(i,9)) *cos((90-dt(i,10))* (pi/180));
    V(i,1)= (dt(i,9)) *sin((90-dt(i,10))* (pi/180));
end

% ---------converting conductivity to real value --------
disp('Converting conductivity to real values.')
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

% ================================================
% Start of Correction of Wrapped Conductivity data

fig_handle=figure;
dl = length(jd);
plot(jd-jd(1),cond,'.-')
set(gca,'Ygrid','on')
set(gca,'Xgrid','on')
ylabel('Cond. (mS/cm)')
xlabel('Time [Days]')
set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
already_plotted=1;

% ask question if want to apply correction
apply_correction=questdlg('Do you want to manually correct for wrapped conductivity data?','Wrapped Conductivity?','Yes','No','No');

% if yes launch selection routine
if strmatch('Yes',apply_correction)
    correct_again='Yes';
    while strmatch(correct_again,'Yes')
        up_or_down='Reselect Data';
        while strmatch(up_or_down,'Reselect Data')
            % modifying code to use data brushing instead of selecting
            % start and end points
            disp('Select batch of data with data brushing technique, create variable "data", and then close plot.')
            data_to_correct=data_point2(fig_handle);
            data_index=zeros(length(data_to_correct),1);
            

            for ij=1:length(data_to_correct(:,1))
                abc = find(jd-jd(1)==data_to_correct(ij,1));
                data_index(ij)=abc;
                
            end
            
            % replot graph showing data to correct
            figure(fig_handle);
            hold off
            dl = length(jd);
            plot(jd-jd(1),cond,'.-')
            hold on
            plot(jd(data_index)-jd(1),cond(data_index),'*r')
            plot([jd(1)-jd(1)-100 jd(dl)-jd(1)+100],[c_upp c_upp],'-g');
            plot([jd(1)-jd(1)-100 jd(dl)-jd(1)+100],[c_low c_low],'-g');
            set(gca,'Ygrid','on')
            set(gca,'Xgrid','on')
            ylabel('Cond. (mS/cm)')
            xlabel('Time [Days]')
            set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
            
            up_or_down=questdlg('Add or Subtract wrapping correction to data?','Add or Subtract Correction?','Add','Subtract','Reselect Data','Reselect Data');
            if strmatch(up_or_down,'Reselect Data')
                delete(data_to_correct);
            elseif strmatch(up_or_down,'Add')
                cond(data_index)=cond(data_index)+(c_upp-c_low);
            elseif strmatch(up_or_down,'Subtract')
                cond(data_index)=cond(data_index)-(c_upp-c_low);
            end
        end
%         while strmatch(up_or_down,'Reselect Data')
%             % Selection of start point to correct
%             %help_box=helpdlg('Select START point of data to apply wrapped offset, and then close plot.','Instructions');
%             disp(' ')
%             disp('Select START point of data to apply wrapped offset, and then close plot.')
%             %waitfor(help_box)
%             start_point=data_point(fig_handle,1);
% 
%             % replot graph, including start_point
%             fig_handle=figure;
%             dl = length(jd);
%             plot(jd-jd(1),cond,'.-')
%             hold on
%             plot([jd(1)-jd(1) jd(dl)-jd(1)],[c_upp c_upp],'-g');
%             plot([jd(1)-jd(1) jd(dl)-jd(1)],[c_low c_low],'-g');
%             set(gca,'Ygrid','on')
%             set(gca,'Xgrid','on')
%             ylabel('Cond. (mS/cm)')
%             xlabel('Time [Days]')
%             set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
% 
%             plot(start_point(1),start_point(2),'or')
%             text(start_point(1),start_point(2),'  Start Point','Color','r')
%             
%             % Selection of end point to correct
%             %help_box=helpdlg('Select END point of data to apply wrapped offset, and then close plot.','Instructions');
%             disp(' ')
%             disp('Select END point of data to apply wrapped offset, and then close plot.')
%             %waitfor(help_box)
%             end_point=data_point(fig_handle,2);
% 
%             start_index=find(jd-jd(1)==start_point(1));
%             end_index=find(jd-jd(1)==end_point(1));
% 
%             % replot graph, including start_point and end point and highlight data
%             % to be corrected
%             fig_handle=figure;
%             dl = length(jd);
%             plot(jd-jd(1),cond,'.-')
%             hold on
%             plot([jd(1)-jd(1) jd(dl)-jd(1)],[c_upp c_upp],'-g');
%             plot([jd(1)-jd(1) jd(dl)-jd(1)],[c_low c_low],'-g');
%             set(gca,'Ygrid','on')
%             set(gca,'Xgrid','on')
%             ylabel('Cond. (mS/cm)')
%             xlabel('Time [Days]')
%             set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
%             
%             data_to_correct=plot(jd(start_index:end_index)-jd(1),cond(start_index:end_index),'.-r');
%             text1=text(start_point(1),start_point(2),'  Start Point','Color','r');
%             text2=text(end_point(1),end_point(2),'  End Point','Color','r');
% 
%             up_or_down=questdlg('Add or Subtract wrapping correction to data?','Add or Subtract Correction?','Add','Subtract','Reselect Data','Reselect Data');
%             if strmatch(up_or_down,'Reselect Data')
%                 delete(data_to_correct); delete(text1); delete(text2);
%             elseif strmatch(up_or_down,'Add')
%                 cond(start_index:end_index)=cond(start_index:end_index)+(c_upp-c_low);
%             elseif strmatch(up_or_down,'Subtract')
%                 cond(start_index:end_index)=cond(start_index:end_index)-(c_upp-c_low);
%             end
%         end
        close(fig_handle)
        fig_handle=figure;
        plot(jd-jd(1),cond,'.-')
        hold on
        plot([jd(1)-jd(1)-100 jd(dl)-jd(1)+100],[c_upp c_upp],'-g');
        plot([jd(1)-jd(1)-100 jd(dl)-jd(1)+100],[c_low c_low],'-g');
        set(gca,'Ygrid','on')
        set(gca,'Xgrid','on')
        ylabel('Cond. (mS/cm)')
        xlabel('Time [Days]')
        set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
        
        correct_again=questdlg('Do you want to correct more data from this plot?','Correct more data?','Yes','No','No');
    end
else close(fig_handle)
end



% END OF CORRECTION OF WRAPPED CONDUCTIVITY DATA
% ==============================================

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
  
  % need to change units differently for serial numbers above 500 as
  % Aanderaa software converts to kPa and not MPa as previously.
  if vec<500
    plot(jd-jd(1),dt(:,13)*100)       % units transform from MPa to db
  else
    plot(jd-jd(1),dt(:,13)/10)       % units transform from kPa to db
  end
  
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

figure(30);

 dl = length(jd);
 
 
  
subplot (2,1,1) % U
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

subplot (2,1,2) % V
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

subplot(2,1,1)

title(['RCM11:  ',num2str(SerialNumber),'   Depth:  ',num2str(z),' m']) 
orient tall
eval(['print -dps ',outfile,'.ps'])


%-------------------------------------
% ---- save to rodb format
%-------------------------------------

disp(['writing data to ',outfile]) 

TIME = gregorian(jd);

%data = [TIME(:,1:3) hms2h(TIME(:,4:6)) dt(:,8) U V dt(:,11:12) dt(:,13)*100 dt(:,14:15)]; 
if vec<500 % if serial number <500 pressure units in MPa
    data = [TIME(:,1:3) hms2h(TIME(:,4:6)) dt(:,8) U V dt(:,11) cond dt(:,13)*100 dt(:,14:15)]; 
else  % % if serial number >500 pressure untis in kPa
    data = [TIME(:,1:3) hms2h(TIME(:,4:6)) dt(:,8) U V dt(:,11) cond dt(:,13)/10 dt(:,14:15)]; 
end

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

end

function data=data_point(fig_handle,start_or_end)
global data_value

figure(fig_handle)
datacursormode on
dcm_obj = datacursormode(fig_handle);
set(dcm_obj,'UpdateFcn',@myupdatefcn)

    function text = myupdatefcn(empt,event_obj)
        data_value = get(event_obj,'Position');
        if start_or_end==1
            text='Start Data Point';
        elseif start_or_end==2
            text='End Data Point';
        end
    end
waitfor(fig_handle)
data=data_value;
end

function data=data_point2(fig_handle)
figure(fig_handle)
brush(fig_handle);
pause
% need to export data from brushing to variable "data" for outpur from
% function
end
