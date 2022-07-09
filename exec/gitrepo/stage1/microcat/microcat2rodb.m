% function microcat2rodb_5('infile','outfile','infofile',fidlog,[graphics],[toffset])
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
% kanzow   26/12/2000 Xmas edition
%          28/03/2001 pressure option added
%          20/04/2005 CSIRO seawater routine for computation of salinity
%          24/04/2005 recorded time checked against mooring deployment time
%          25/04/2005 input variable 'toffset' added
%          27/10/2009 - DR added functionality for .cnv format files for
%             microcat firmware 3.0 and above. and created microcat2rodb_3
%             from microcat2rodb_2_002
%          03.04.2010 - ZBS changed the method for setting 'ylim'
%             for graphics=='y' to use prctile.
%          20/11/2012 - DR added lines to tell which sample number has
%             format problems following issues with bad scans in the middle
%             of records on D382
%          30/04/2014 - DR added compatability for MicroCAT ODOs
%          02/05/2014 - BIM check added to make sure microcat is selected using id and serial number
%             e.g. nortek and microcat share serial numbers (6805)
%          01/11/2015 - GDM added smart indexing for .cnv input files so
%             that if the script will search the order that that the seabird
%             datcnv program output the variables and not just assume they are
%             t:c:p:ot:ox:time
%             Minor change to allow for 5 digit serial numbers
%             Major cleanup of the plotting (90% reduction in num. of codelines)
%          05/05/2017 - LH added the possibility to read short data file
%          (length of header ~= length of data)
%                     - Add the possibility to read time as Julian day (timeJV2)  
%
function microcat2rodb(infile,outfile,infofile,fidlog,graphics,toffset)

if nargin < 4
    disp('not enough input arguments')
elseif nargin == 4
    graphics = 'n'
    disp('graphics will not be  displayed')
end

if nargin < 6
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
    
    answer = questdlg('Overwrite previous output?', ...
	'Output', ...
	'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            disp('Overwriting previous file.')
        case 'No'
            disp('Data conversion stopped.')
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
    %
    % gdm, dy039, incremented to account for new 5 digit serial numbers
    if ~isempty(a)
        SerialNumber = str2num(zeile(a+17:a+22));
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
% There is a bug here when runnning this code for older instrumnets (e.g.
% instruments with early serial numbers - 4608). The code below reads the
% number of lines in the .cnv file and estimaetes a start postion. The
% issue is that older microcats have significantly less header lines so the
% start position is well below the threshold number.I have added an elseif
% for just now as a sticking plaster - LAD 06/10/2020

% detect data column length
if length(retx)<700 %in case of short record from test (added by Loic H on DY078)
    if SerialNumber==4608
            [XXX,data_length] = max(hist(diff(retx(73:end)),1:300));% length of data columns 
    else
            [XXX,data_length] = max(hist(diff(retx(304:end)),1:300));% length of data columns 

    end
else
    [XXX,data_length] = max(hist(diff(retx),1:300));% length of data columns   
end
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
    fprintf(fidlog,'conversion stopped - severe deviation from input format: %d\n',sede+ii0(1));
    msgbox(['MC',sprintf('%4.4d','SerialNumber'),'  ',sprintf('severe deviation from input format: %d\n',sede+ii0(1))],'conversion stopped')
%    return
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

% ---- check order of variables if cnv file ----

% --- short stub inserted by gdm, dy039 to ensure the indices match the cnv file ---
% --- aditional conditions added by bim dy039, to enable DR PC conversion files to be used.

if cnv
    namei = findstr(zeile,'name');
    spani = findstr(zeile,'span');
    toti = [namei spani(1)];
    
    for k = 1:length(toti)-1
        varstr=zeile(toti(k):toti(k+1)-5);
        
        if ~isempty( strfind(varstr(10), 't') ) & ~isempty( strfind(varstr, '90') )
            
            tempi = k;
            continue;
        elseif ~isempty(strfind(varstr(10), 'c')) & ~isempty(strfind(varstr, '0S/m'))
            
            condi = k;
            continue;
        elseif ~isempty(strfind(varstr, 'pr')) |  ~isempty(strfind(varstr, 'prdM'))
            
            presi = k;
            continue;
        elseif ~isempty(strfind(varstr, 'sbeoxTC'))
            
            oxyti = k;
            continue;
        elseif ~isempty(strfind(varstr, 'sbeopoxMm/Kg'))
            
            oxyi = k;
            continue;
        elseif ~isempty(strfind(varstr, 'timeK'))|  ~isempty(strfind(varstr, 'timeS'))
            timeformat = 'sec';
            timei = k;
            continue;
        elseif ~isempty(strfind(varstr, 'timeJ'))
            timeformat = 'JD';
            timei = k;
            continue;            
        elseif ~isempty(strfind(varstr(10), 'c')) & ~isempty(strfind(varstr, '0mS/cm'))
            condi = k;
            dt(:,condi)=dt(:,condi)/10;
            
            continue;
            
        else
        end
        
    end;
end;


% % unique fix for for ar30 cal dip where jday for 3276 microcat restarted at the end of the first year of deployment 
% if (strfind(infile,'ar30') & strfind(infile,'3276_data'));
%     keyboard
%     findt(diff(dt(:,4)))
% end

% ---- determine time -------

dtl        = size(dt,1);
cnv
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
  if ~isempty(strfind(timeformat, 'sec'))
    if size(dt,2)==5 % SMP with pressure
        cnv_secs=dt(:,timei); % seconds since 1st Jan 2000.
        jd=cnv_secs/(60*60*24)+julian(2000,1,1,0) - toffset;
        disp(['toffset : ',num2str(toffset)]);
        gtime=gregorian(jd);
        HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6));
        Start_Date = gtime(1,[1 2 3]);
        End_Date = gtime(dtl,[1 2 3]);
    elseif size(dt,2)==7 % SMP-ODO with pressure
        cnv_secs=dt(:,timei); % seconds since 1st Jan 2000.
        jd=cnv_secs/(60*60*24)+julian(2000,1,1,0) - toffset;
        disp(['toffset : ',num2str(toffset)]);
        gtime=gregorian(jd);
        HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6));
        Start_Date = gtime(1,[1 2 3]);
        End_Date = gtime(dtl,[1 2 3]);
    end
  elseif ~isempty(strfind(timeformat, 'JD'))
    if size(dt,2)==5 % SMP with pressure
        cnv_days=dt(:,4); % seconds since 1st Jan 2000.
        jd=cnv_days+julian(toffset,1,1,0) - 1;
        disp(['dateoffset : ',num2str(toffset)]);
        gtime=gregorian(jd);
        HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6));
        Start_Date = gtime(1,[1 2 3]);
        End_Date = gtime(dtl,[1 2 3]);
    elseif size(dt,2)==7 % SMP-ODO with pressure
        cnv_days=dt(:,4); % seconds since 1st Jan 2000.
        jd=cnv_days+julian(toffset,1,1,0) - 1;
        disp(['dateoffset : ',num2str(toffset)]);
        gtime=gregorian(jd);
        HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6));
        Start_Date = gtime(1,[1 2 3]);
        End_Date = gtime(dtl,[1 2 3]);                
    end      
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

%%
%-------------------------------------
% ---- save to rodb format
%-------------------------------------
disp(['writing data to ',outfile])

TIME = gregorian(jd);

if cnv==0 %.asc file
    data = [TIME(:,1:3) hms2h(TIME(:,4:6)) dt(:,1) dt(:,2)*10];
    
    if  size(dt,2) ==9           % with pressure
        data = [data dt(:,3)];
        cols = colsp;
        fort = fortp;
    end
else % .cnv file
    data = [TIME(:,1:3) hms2h(TIME(:,4:6)) dt(:,tempi) dt(:,condi)*10];
    
    if size(dt,2)==5
        data = [data dt(:,presi)];
        cols = colsp;
        fort = fortp;
    elseif size(dt,2)==7 % ODO sensor
        data = [data dt(:,[presi oxyti oxyi])];
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

%% % --- New graphics section --- 
% (because why would you need 400+ lines of code to do what 40 odd will do?)
sz = size(data); % format saved by rodb, yy:mm:dd:hh:t:c:p(:ot:o2)
jd = datenum(data(:,1),data(:,2),data(:,3),data(:,4),0,0);
dl = length(jd);
nplots = sz(2) - 4; % variables other than the four time vars

for k = 1:nplots
    subplot(nplots,1,k);
    if k==1
          title(['Microcat:  ',num2str(SerialNumber),'   Depth:  ',num2str(z),' m']) 
    end
    hold off;
    plot(jd-jd(1),data(:,4+k));
    set(gca,'Ygrid','on')
    set(gca,'Xgrid','on')
    
    mt = median(data(:,4+k));
    if k == 1
        ylabel('Temp.')
        fprintf(fidlog,'Median Temperature %5.2f\n',mt);
    elseif k == 2
        ylabel('Cond.')
        fprintf(fidlog,'Median Conductivity %5.2f\n',mt);
    elseif k == 3
        ylabel('Press.')
        fprintf(fidlog,'Median Pressure %5.2f\n',mt);
    elseif k == 4
        ylabel('Oxy Temp.')
        fprintf(fidlog,'Median Oxygen Temperature %5.2f\n',mt);
    elseif k == 5
        ylabel('Oxygen')
        fprintf(fidlog,'Median Oxygen %5.2f\n',mt);
    end

    set(gca,'ylim',prctile(data(:,4+k),[1 99])+[-1 1]*std(data(:,4+k))/10)
    
    set(gca,'xlim',[jd(1) jd(dl) ]-jd(1))
    if  graphics == 'w'
        set(gca,'ylim',[min(data(:,4+k))  min([max(data(:,4+k)) 30])])
    end
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
    tx =text(xl(1) + diff(xl)*.1,yl(1) + diff(yl)*.85,['median: ',sprintf('%3.2f',mt)]) ;
    set(tx,'FontWeight','bold','FontSize',11)
end
  orient tall; grid on

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

if sz(1) ~= ex_samples
    fprintf('Number of samples: %d; expected: %d \n',sz(1),ex_samples);
    uiwait(warndlg(['Number of samples: ' num2str(sz(1)) ' expected: ' num2str(ex_samples)],'Warning'));
end

%tk 24/04/05
if jd(end) < bottomstop
    fprintf(fidlog,'\n\n W A R N I N G: Record already ends before end of mooring deployment period: \n Recorded time could be wrong!! \n');
end
if jd(1) > bottomstart
    fprintf(fidlog,'\n\n W A R N I N G: Record only starts after beginning of mooring deployment period: \n Recorded time could be wrong!! \n');
end






