% function varargout=stats_table(moor,'procpath','layout','non-verbose')
% 
% function to generate a simple stats table for data reports
% 
% required inputs:-
%   moor: complete mooring name as string. e.g. 'wb1_1_200420'
%
% optional inputs:-
%   layout: orientation of figure portrait/lanscape (default = portrait)
%           input of 'landscape' or 'portrait'
%           e.g. pressure_overlay('wb1_1_200420','layout','landscape')
%   procpath: specify exact procpath
%           e.g. pressure_overlay('wb1_1_200420','procpath','/Volumes/noc/mpoc/hydro/rpdmoc/rapid/data/moor/proc/')
%   outpath: specify exact output path
%   non-verbose: a mode to output the stats for other routines without
%           saving an ascii file. If called in this mode there will be an
%           output to the function consisiting of the Microcat serial
%           numbers and statistics
%
%   output: an ascii text file containing the summary statistics for the
%   mooring - naming convention is moor_stats.asc where moor is the mooring
%   name as input to the function
%
% functions called:-
%   rodbload, julian, auto_filt, gregorian
%   from .../exec/moor/tools and .../exec/moor/rodb paths
% 
% Routine written by Darren Rayner February 2008.
%
% Doesn't do anything with MMP or ADCP data.
%
% Loic Houpert, 4/10/16, add output path 
%
%***************************
%***************************
% MAJOR PROBLEM
%***************************
%***************************
% .use files do not have bad data changed to -9999 so statistics will be
% badly skewed by bad data such as when a microcat has a pressure sensor
% failure. Need to check with Torsten if there is a better input file to
% use.




function varargout=stats_table(moor,varargin)
if nargin <1
    help stats_table
    return
end

% check for optional arguments
a=strmatch('layout',varargin,'exact');
if a>0
    layout=char(varargin(a+1));
else
    layout='portrait';
end

a=strmatch('procpath',varargin,'exact');
if a>0
    procpath=char(varargin(a+1));
else
    data_report_tools_dir=which('data_report_tools');
    b=strfind(data_report_tools_dir,'/');
    data_report_tools_dir=data_report_tools_dir(1:b(end));
    procpath=[data_report_tools_dir '../../../moor/proc/']; %DR changed to be relative paths now that data_report_tools are on the network 19/2/12
end
b=strmatch('outpath',varargin,'exact');
if b>0
    outpath=char(varargin(b+1));
end


a=strmatch('non-verbose',varargin,'exact');
if a>0
    non_verbose=1;
else
    non_verbose=0;
end

if isunix
    procpath=[procpath '/'];
    infofile=[procpath,moor,'/',moor,'info.dat'];
elseif ispc
    procpath=[procpath '\'];
    infofile=[procpath,moor,'\',moor,'info.dat'];
end

%-----------------------------------------------------
% Load vectors of mooring information
% id instrument id, sn serial number, z nominal depth of each instrument
% s_t, e_t, s_d, e_d start and end times and dates
% lat lon mooring position, wd corrected water depth (m)
% mr mooring name
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,...
    'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');


% JULIAN Convert Gregorian date to Julian day.
% JD = JULIAN(YY,MM,DD,HH) or JD = JULIAN([YY,MM,DD,HH]) returns the Julian
% day number of the calendar date specified by year, month, day, and decimal
% hour.
% JD = JULIAN(YY,MM,DD) or JD = JULIAN([YY,MM,DD]) if decimal hour is absent,
% it is assumed to be zero.
% Although the formal definition holds that Julian days start and end at
% noon, here Julian days start and end at midnight. In this convention,
% Julian day 2440000 began at 00:00 hours, May 23, 1968.
jd_start = julian([s_d' hms2h([s_t;0]')]);
jd_end   = julian([e_d' hms2h([e_t;0]')]);

if non_verbose==0
    disp(['z : instrument id : serial number'])
    for i = 1:length(id)
        disp([z(i),id(i),sn(i)])
    end
end

% find the index number of Microcats
iiMC = find(id == 337);
vecMC = sn(iiMC);
% find the index number of RBRs
iiRBR = find(id == 330);
vecRBR = sn(iiRBR);
% find the index number of Idronauts
iiIDR = find(id == 339);
vecIDR = sn(iiIDR);
% find the index number of S4s
iiS4 = find(id == 302);
vecS4 = sn(iiS4);
% and find index number of RCM11s
iiRCM11 = find(id == 310);
vecRCM11 = sn(iiRCM11);
% and find index number of Sontek Argonauts
iiARG = find(id == 366);
vecARG = sn(iiARG);
% and find index number of Nortek Aquadopps
iiNOR = find((id == 368|id==370));
vecNOR = sn(iiNOR);
% and find index number of TRDI DVSs
iiDVS = find(id == 369);
vecDVS = sn(iiDVS);
% and find index number of both Seabird and Ixsea BPRs
iiBPR = find((id == 465)|(id==470));
vecBPR = sn(iiBPR);

depths(:,1) = id([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiDVS;iiBPR]);
depths(:,2) = z([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiDVS;iiBPR]);
depths(:,3) = sn([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiDVS;iiBPR]);
depths=sortrows(depths,2);
iiiMC=find(depths(:,1)==337);
iiiRBR=find(depths(:,1)==330);
iiiIDR=find(depths(:,1)==339);
iiiS4=find(depths(:,1)==302);  
iiiRCM11=find(depths(:,1)==310);  
iiiARG=find(depths(:,1)==366 | depths(:,1)==366337); 
iiiNOR=find(depths(:,1)==368|depths(:,1)==370);
iiiDVS=find(depths(:,1)==369);
iiiBPR=find(depths(:,1)==465 | depths(:,1)==470);
iii=[iiiS4;iiiRCM11;iiiARG;iiiNOR;iiiDVS;iiiMC;iiiRBR;iiiIDR;iiiBPR];

% -----------------------------------
% START OF READING IN INSTRUMENT DATA
% -----------------------------------
% modified to only read in MicroCat data if calling in non-verbose mode
%--------------------------------------
% Now read in Microcat data if required
%--------------------------------------
if iiMC>0
    j=1;
    
    % loop to read one file at a time

    for i=1:length(vecMC);
       serialno = vecMC(i);
       if non_verbose==0
           disp('*************************************************************')
           disp(['Reading MICROCAT - ',num2str(serialno)])
           disp('*************************************************************')
       end
       if isunix
           infile = [procpath,moor,'/microcat/',moor,'_',sprintf('%0.4d',vecMC(i)),'.use'];
       elseif ispc
           infile = [procpath,moor,'\microcat\',moor,'_',sprintf('%0.4d',vecMC(i)),'.use'];
       end

       % check if file exists
       fileopen=fopen(infile,'r');
       
       if fileopen>0
           % read data into vectors and then into structure array

           [yy,mm,dd,hh,t,c,p] = ...
               rodbload(infile,'yy:mm:dd:hh:t:c:p');
           jd=julian(yy,mm,dd,hh);

           bad_data=find(p==-9999); p(bad_data)=NaN;
           bad_data=find(t==-9999); t(bad_data)=NaN;
           bad_data=find(c==-9999); c(bad_data)=NaN;

           eval_string(iiiMC(j))={['MC_' num2str(serialno)]};

           eval([char(eval_string(iiiMC(j))) '.jd=jd;']);
           eval([char(eval_string(iiiMC(j))) '.p=p;']);
           eval([char(eval_string(iiiMC(j))) '.t=t;']);
           eval([char(eval_string(iiiMC(j))) '.c=c;']);
           
           samples=length(jd);
           start_date=gregorian(jd(1));
           end_date=gregorian(jd(end));
           eval([char(eval_string(iiiMC(j))) '.samples=samples;']);
           eval([char(eval_string(iiiMC(j))) '.start_date=start_date;']);
           eval([char(eval_string(iiiMC(j))) '.end_date=end_date;']);
           
           % calculate statistics
           tmean = nanmean(t); cmean = nanmean(c); pmean = nanmean(p);
           tstd= nanstd(t); cstd= nanstd(c); pstd= nanstd(p);
           tmax = max(t); cmax = max(c); pmax = max(p);
           tmin = min(t); cmin = min(c); pmin = min(p);
           
           eval([char(eval_string(iiiMC(j))) '.pmean=pmean;']);
           eval([char(eval_string(iiiMC(j))) '.tmean=tmean;']);
           eval([char(eval_string(iiiMC(j))) '.cmean=cmean;']);
           eval([char(eval_string(iiiMC(j))) '.pmin=pmin;']);
           eval([char(eval_string(iiiMC(j))) '.tmin=tmin;']);
           eval([char(eval_string(iiiMC(j))) '.cmin=cmin;']);
           eval([char(eval_string(iiiMC(j))) '.pmax=pmax;']);
           eval([char(eval_string(iiiMC(j))) '.tmax=tmax;']);
           eval([char(eval_string(iiiMC(j))) '.cmax=cmax;']);
           eval([char(eval_string(iiiMC(j))) '.pstd=pstd;']);
           eval([char(eval_string(iiiMC(j))) '.tstd=tstd;']);
           eval([char(eval_string(iiiMC(j))) '.cstd=cstd;']);
           
       else
           disp(['File does not exist!'])
           disp(['infile = ' infile])
           
       end
       j=j+1;
    end
    
end
if non_verbose==0
    %--------------------------------------
    % Now read in Idronaut data if required
    %--------------------------------------
    if iiIDR>0
        j=1;

        % loop to read one file at a time

        for i=1:length(vecIDR);
           serialno = vecIDR(i);
           disp('*************************************************************')
           disp(['Reading IDRONAUT - ',num2str(serialno)])
           disp('*************************************************************')

           if isunix
               infile = [procpath,moor,'/idr/',moor,'_',sprintf('%0.4d',vecIDR(i)),'.use'];
           elseif ispc
               infile = [procpath,moor,'\idr\',moor,'_',sprintf('%0.4d',vecIDR(i)),'.use'];
           end

           % read data into vectors and then into structure array

           % check if file exists
           fileopen=fopen(infile,'r');


           if fileopen>0

               [yy,mm,dd,hh,t,c,p,] = ...
                   rodbload(infile,'yy:mm:dd:hh:t:c:p');
               jd=julian(yy,mm,dd,hh);

               bad_data=find(p==-9999); p(bad_data)=NaN;
               bad_data=find(t==-9999); t(bad_data)=NaN;
               bad_data=find(c==-9999); c(bad_data)=NaN;

               eval_string(iiiIDR(j))={['IDR_' num2str(serialno)]};

               eval([char(eval_string(iiiIDR(j))) '.jd=jd;']);
               eval([char(eval_string(iiiIDR(j))) '.p=p;']);
               eval([char(eval_string(iiiIDR(j))) '.t=t;']);
               eval([char(eval_string(iiiIDR(j))) '.c=c;']);

               samples=length(jd);
               start_date=gregorian(jd(1));
               end_date=gregorian(jd(end));
               eval([char(eval_string(iiiIDR(j))) '.samples=samples;']);
               eval([char(eval_string(iiiIDR(j))) '.start_date=start_date;']);
               eval([char(eval_string(iiiIDR(j))) '.end_date=end_date;']);

               % calculate statistics
               tmean = nanmean(t); cmean = nanmean(c); pmean = nanmean(p);
               tstd= nanstd(t); cstd= nanstd(c); pstd= nanstd(p);
               tmax = max(t); cmax = max(c); pmax = max(p);
               tmin = min(t); cmin = min(c); pmin = min(p);

               eval([char(eval_string(iiiIDR(j))) '.pmean=pmean;']);
               eval([char(eval_string(iiiIDR(j))) '.tmean=tmean;']);
               eval([char(eval_string(iiiIDR(j))) '.cmean=cmean;']);
               eval([char(eval_string(iiiIDR(j))) '.pmin=pmin;']);
               eval([char(eval_string(iiiIDR(j))) '.tmin=tmin;']);
               eval([char(eval_string(iiiIDR(j))) '.cmin=cmin;']);
               eval([char(eval_string(iiiIDR(j))) '.pmax=pmax;']);
               eval([char(eval_string(iiiIDR(j))) '.tmax=tmax;']);
               eval([char(eval_string(iiiIDR(j))) '.cmax=cmax;']);
               eval([char(eval_string(iiiIDR(j))) '.pstd=pstd;']);
               eval([char(eval_string(iiiIDR(j))) '.tstd=tstd;']);
               eval([char(eval_string(iiiIDR(j))) '.cstd=cstd;']);
           else
               disp(['File does not exist!'])
               disp(['infile = ' infile])

           end
           j=j+1;
        end
    end
    %--------------------------------------
    % Now read in RBR data if required
    %--------------------------------------
    if iiRBR>0
        j=1;

        % loop to read one file at a time

        for i=1:length(vecRBR);
           serialno = vecRBR(i);
           disp('*************************************************************')
           disp(['Reading RBR - ',num2str(serialno)])
           disp('*************************************************************')

           if isunix
               infile = [procpath,moor,'/rbr/',moor,'_',sprintf('%0.4d',vecRBR(i)),'.use'];
           elseif ispc
               infile = [procpath,moor,'\rbr\',moor,'_',sprintf('%0.4d',vecRBR(i)),'.use'];
           end

           % read data into vectors and then into structure array

           % check if file exists
           fileopen=fopen(infile,'r');

           if fileopen>0

               [yy,mm,dd,hh,p,t,c] = ...
                   rodbload(infile,'yy:mm:dd:hh:p:t:c');
               jd=julian(yy,mm,dd,hh);

               bad_data=find(p==-9999); p(bad_data)=NaN;
               bad_data=find(t==-9999); t(bad_data)=NaN;
               bad_data=find(c==-9999); c(bad_data)=NaN;

               eval_string(iiiRBR(j))={['RBR_' num2str(serialno)]};

               eval([char(eval_string(iiiRBR(j))) '.jd=jd;']);
               eval([char(eval_string(iiiRBR(j))) '.p=p;']);
               eval([char(eval_string(iiiRBR(j))) '.t=t;']);
               eval([char(eval_string(iiiRBR(j))) '.c=c;']);

               samples=length(jd);
               start_date=gregorian(jd(1));
               end_date=gregorian(jd(end));
               eval([char(eval_string(iiiRBR(j))) '.samples=samples;']);
               eval([char(eval_string(iiiRBR(j))) '.start_date=start_date;']);
               eval([char(eval_string(iiiRBR(j))) '.end_date=end_date;']);

               % calculate statistics
               tmean = nanmean(t); cmean = nanmean(c); pmean = nanmean(p);
               tstd= nanstd(t); cstd= nanstd(c); pstd= nanstd(p);
               tmax = max(t); cmax = max(c); pmax = max(p);
               tmin = min(t); cmin = min(c); pmin = min(p);

               eval([char(eval_string(iiiRBR(j))) '.pmean=pmean;']);
               eval([char(eval_string(iiiRBR(j))) '.tmean=tmean;']);
               eval([char(eval_string(iiiRBR(j))) '.cmean=cmean;']);
               eval([char(eval_string(iiiRBR(j))) '.pmin=pmin;']);
               eval([char(eval_string(iiiRBR(j))) '.tmin=tmin;']);
               eval([char(eval_string(iiiRBR(j))) '.cmin=cmin;']);
               eval([char(eval_string(iiiRBR(j))) '.pmax=pmax;']);
               eval([char(eval_string(iiiRBR(j))) '.tmax=tmax;']);
               eval([char(eval_string(iiiRBR(j))) '.cmax=cmax;']);
               eval([char(eval_string(iiiRBR(j))) '.pstd=pstd;']);
               eval([char(eval_string(iiiRBR(j))) '.tstd=tstd;']);
               eval([char(eval_string(iiiRBR(j))) '.cstd=cstd;']);
           else
               disp(['File does not exist!'])
               disp(['infile = ' infile])

           end
           j=j+1;
        end
    end

    % --------------------
    % Read in S4 data if required.
    % --------------------
    if iiS4>0 

        % loop to read one file at a time
        j=1;
        for i=1:length(vecS4);
           serialno = vecS4(i);
           disp('*************************************************************')
           disp(['Reading S4 - ',num2str(serialno)])
           disp('*************************************************************')

           if isunix
               infile = [procpath,moor,'/s4/',moor,'_',sprintf('%4.4d',vecS4(i)),'.use'];
           elseif ispc
               infile = [procpath,moor,'\s4\',moor,'_',sprintf('%4.4d',vecS4(i)),'.use'];
           end

           % check if file exists
           fileopen=fopen(infile,'r');

           if fileopen>0
               % read data into vectors and then into structure array

               [yy,mm,dd,hh,u,v,t,c,p,hdg] = rodbload(infile,'yy:mm:dd:hh:u:v:t:c:p:hdg');
               jd=julian(yy,mm,dd,hh);

               bad_data=find(p==-9999); p(bad_data)=NaN;
               bad_data=find(t==-9999); t(bad_data)=NaN;
               bad_data=find(c==-9999); c(bad_data)=NaN;
               bad_data=find(u==-9999); u(bad_data)=NaN;
               bad_data=find(v==-9999); v(bad_data)=NaN;


               eval_string(iiiS4(j))={['S4_' num2str(serialno)]};

               eval([char(eval_string(iiiS4(j))) '.jd=jd;']);
               eval([char(eval_string(iiiS4(j))) '.p=p;']);
               eval([char(eval_string(iiiS4(j))) '.t=t;']);
               eval([char(eval_string(iiiS4(j))) '.c=c;']);
               eval([char(eval_string(iiiS4(j))) '.u=u;']);
               eval([char(eval_string(iiiS4(j))) '.v=v;']);

               samples=length(jd);
               start_date=gregorian(jd(1));
               end_date=gregorian(jd(end));
               eval([char(eval_string(iiiS4(j))) '.samples=samples;']);
               eval([char(eval_string(iiiS4(j))) '.start_date=start_date;']);
               eval([char(eval_string(iiiS4(j))) '.end_date=end_date;']);


               % calculate speed and direction 
               spd = sqrt(u.^2 + v.^2);
               eval([char(eval_string(iiiS4(j))) '.spd=spd;']);
               direction=atan(u./v)*180/pi;
               d = find(v<0);
               direction(d)=direction(d)+180;
               d = find(v>=0 & u<0);
               direction(d)=direction(d)+360;
               eval([char(eval_string(iiiS4(j))) '.dir=direction;']);

               % calculate statistics
               tmean = nanmean(t); cmean = nanmean(c); pmean = nanmean(p);
               tstd = nanstd(t); cstd = nanstd(c); pstd = nanstd(p);
               tmax = max(t); cmax = max(c); pmax = max(p);
               tmin = min(t); cmin = min(c); pmin = min(p);

               umin = min(u); vmin = min(v);
               umax = max(u); vmax = max(v);
               umean = nanmean(u); vmean = nanmean(v);
               ustd = nanstd(u); vstd = nanstd(v);
               spdmin = min(spd); dirmin = min(direction);
               spdmax = max(spd); dirmax = max(direction);

               eval([char(eval_string(iiiS4(j))) '.pmean=pmean;']);
               eval([char(eval_string(iiiS4(j))) '.tmean=tmean;']);
               eval([char(eval_string(iiiS4(j))) '.cmean=cmean;']);
               eval([char(eval_string(iiiS4(j))) '.umean=umean;']);
               eval([char(eval_string(iiiS4(j))) '.vmean=vmean;']);
               eval([char(eval_string(iiiS4(j))) '.pmin=pmin;']);
               eval([char(eval_string(iiiS4(j))) '.tmin=tmin;']);
               eval([char(eval_string(iiiS4(j))) '.cmin=cmin;']);
               eval([char(eval_string(iiiS4(j))) '.umin=umin;']);
               eval([char(eval_string(iiiS4(j))) '.vmin=vmin;']);
               eval([char(eval_string(iiiS4(j))) '.spdmin=spdmin;']);
               eval([char(eval_string(iiiS4(j))) '.dirmin=dirmin;']);
               eval([char(eval_string(iiiS4(j))) '.pmax=pmax;']);
               eval([char(eval_string(iiiS4(j))) '.tmax=tmax;']);
               eval([char(eval_string(iiiS4(j))) '.cmax=cmax;']);
               eval([char(eval_string(iiiS4(j))) '.umax=umax;']);
               eval([char(eval_string(iiiS4(j))) '.vmax=vmax;']);
               eval([char(eval_string(iiiS4(j))) '.spdmax=spdmax;']);
               eval([char(eval_string(iiiS4(j))) '.dirmax=dirmax;']);
               eval([char(eval_string(iiiS4(j))) '.pstd=pstd;']);
               eval([char(eval_string(iiiS4(j))) '.tstd=tstd;']);
               eval([char(eval_string(iiiS4(j))) '.cstd=cstd;']);
               eval([char(eval_string(iiiS4(j))) '.ustd=ustd;']);
               eval([char(eval_string(iiiS4(j))) '.vstd=vstd;']);

               % calculate spd and dir means and stds
               dirmean=atan(umean/vmean)*180/pi;
               if vmean<0
                   dirmean=dirmean+180;
               elseif (vmean>=0 & umean<0);
                   dirmean=dirmean+360;
               end
               eval([char(eval_string(iiiS4(j))) '.dirmean=dirmean;']);

               spdmean=sqrt(umean.^2 + vmean.^2);
               spdstd = nanstd(spd);
               eval([char(eval_string(iiiS4(j))) '.spdmean=spdmean;']);
               eval([char(eval_string(iiiS4(j))) '.spdstd=spdstd;']);

               % for dir STD convert directions to values around mean direction
               % i.e. mean direction becomes 0, and all other directions are relative to
               % that.
               if dirmean<=180
                   dir1=find(direction<=180+dirmean);
                   dir2=find(direction>180+dirmean);
                   dir_new=direction-dirmean;
                   dir_new(dir2)=dir_new(dir2)-360;
                   dirstd=nanstd(dir_new);
               elseif (dirmean>180)
                    dir1=find(direction<=dirmean-180);
                    dir2=find(direction>dirmean-180);
                    dir_new=direction-dirmean;
                    dir_new(dir1)=dir_new(dir1)+360;
                    dirstd=nanstd(dir_new);
               end
               eval([char(eval_string(iiiS4(j))) '.dirstd=dirstd;']);


           else
               disp(['File does not exist!'])
               disp(['infile = ' infile])

           end
           j=j+1;

        end
    end

    %-----------------------------------
    % Now read in RCM11 data if required
    %-----------------------------------
    if iiRCM11>0
        j=1;

        % loop to read one file at a time

        for i=1:length(vecRCM11);
           serialno = vecRCM11(i);
           disp('*************************************************************')
           disp(['Reading RCM11 - ',num2str(serialno)])
           disp('*************************************************************')

           if isunix
               infile = [procpath,moor,'/rcm/',moor,'_',sprintf('%3.3d',vecRCM11(i)),'.use'];
           elseif ispc
               infile = [procpath,moor,'\rcm\',moor,'_',sprintf('%3.3d',vecRCM11(i)),'.use'];
           end

           % check if file exists
           fileopen=fopen(infile,'r');

           if fileopen>0
               % read data into vectors and then into structure array

               [yy,mm,dd,hh,ref,u,v,t,c,p,tlt,mss] = rodbload(infile,'yy:mm:dd:hh:ref:u:v:t:c:p:tlt:mss');
               jd=julian(yy,mm,dd,hh);

               bad_data=find(p==-9999); p(bad_data)=NaN;
               bad_data=find(t==-9999); t(bad_data)=NaN;
               bad_data=find(c==-9999); c(bad_data)=NaN;
               bad_data=find(u==-9999); u(bad_data)=NaN;
               bad_data=find(v==-9999); v(bad_data)=NaN;


               eval_string(iiiRCM11(j))={['RCM11_' num2str(serialno)]};

               eval([char(eval_string(iiiRCM11(j))) '.jd=jd;']);
               eval([char(eval_string(iiiRCM11(j))) '.p=p;']);
               eval([char(eval_string(iiiRCM11(j))) '.t=t;']);
               eval([char(eval_string(iiiRCM11(j))) '.c=c;']);
               eval([char(eval_string(iiiRCM11(j))) '.u=u;']);
               eval([char(eval_string(iiiRCM11(j))) '.v=v;']);

               samples=length(jd);
               start_date=gregorian(jd(1));
               end_date=gregorian(jd(end));
               eval([char(eval_string(iiiRCM11(j))) '.samples=samples;']);
               eval([char(eval_string(iiiRCM11(j))) '.start_date=start_date;']);
               eval([char(eval_string(iiiRCM11(j))) '.end_date=end_date;']);

               % calculate speed and direction 
               spd = sqrt(u.^2 + v.^2);
               eval([char(eval_string(iiiRCM11(j))) '.spd=spd;']);
               direction=atan(u./v)*180/pi;
               d = find(v<0);
               direction(d)=direction(d)+180;
               d = find(v>=0 & u<0);
               direction(d)=direction(d)+360;
               eval([char(eval_string(iiiRCM11(j))) '.dir=direction;']);

               % calculate statistics
               tmean = nanmean(t); cmean = nanmean(c); pmean = nanmean(p);
               tstd = nanstd(t); cstd = nanstd(c); pstd = nanstd(p);
               tmax = max(t); cmax = max(c); pmax = max(p);
               tmin = min(t); cmin = min(c); pmin = min(p);

               umin = min(u); vmin = min(v);
               umax = max(u); vmax = max(v);
               umean = nanmean(u); vmean = nanmean(v);
               ustd = nanstd(u); vstd = nanstd(v);
               spdmin = min(spd); dirmin = min(direction);
               spdmax = max(spd); dirmax = max(direction);

               eval([char(eval_string(iiiRCM11(j))) '.pmean=pmean;']);
               eval([char(eval_string(iiiRCM11(j))) '.tmean=tmean;']);
               eval([char(eval_string(iiiRCM11(j))) '.cmean=cmean;']);
               eval([char(eval_string(iiiRCM11(j))) '.umean=umean;']);
               eval([char(eval_string(iiiRCM11(j))) '.vmean=vmean;']);
               eval([char(eval_string(iiiRCM11(j))) '.pmin=pmin;']);
               eval([char(eval_string(iiiRCM11(j))) '.tmin=tmin;']);
               eval([char(eval_string(iiiRCM11(j))) '.cmin=cmin;']);
               eval([char(eval_string(iiiRCM11(j))) '.umin=umin;']);
               eval([char(eval_string(iiiRCM11(j))) '.vmin=vmin;']);
               eval([char(eval_string(iiiRCM11(j))) '.spdmin=spdmin;']);
               eval([char(eval_string(iiiRCM11(j))) '.dirmin=dirmin;']);
               eval([char(eval_string(iiiRCM11(j))) '.pmax=pmax;']);
               eval([char(eval_string(iiiRCM11(j))) '.tmax=tmax;']);
               eval([char(eval_string(iiiRCM11(j))) '.cmax=cmax;']);
               eval([char(eval_string(iiiRCM11(j))) '.umax=umax;']);
               eval([char(eval_string(iiiRCM11(j))) '.vmax=vmax;']);
               eval([char(eval_string(iiiRCM11(j))) '.spdmax=spdmax;']);
               eval([char(eval_string(iiiRCM11(j))) '.dirmax=dirmax;']);
               eval([char(eval_string(iiiRCM11(j))) '.pstd=pstd;']);
               eval([char(eval_string(iiiRCM11(j))) '.tstd=tstd;']);
               eval([char(eval_string(iiiRCM11(j))) '.cstd=cstd;']);
               eval([char(eval_string(iiiRCM11(j))) '.ustd=ustd;']);
               eval([char(eval_string(iiiRCM11(j))) '.vstd=vstd;']);


               % calculate spd and dir means and std
               dirmean=atan(umean./vmean)*180/pi;
               if vmean<0
                   dirmean=dirmean+180;
               elseif (vmean>=0 & umean<0);
                   dirmean=dirmean+360;
               end
               eval([char(eval_string(iiiRCM11(j))) '.dirmean=dirmean;']);

               spdmean=sqrt(umean.^2 + vmean.^2);
               spdstd = nanstd(spd);
               eval([char(eval_string(iiiRCM11(j))) '.spdmean=spdmean;']);
               eval([char(eval_string(iiiRCM11(j))) '.spdstd=spdstd;']);

               % for dir STD convert directions to values around mean direction
               % i.e. mean direction becomes 0, and all other directions are relative to
               % that.
               if dirmean<=180
                   dir1=find(direction<=180+dirmean);
                   dir2=find(direction>180+dirmean);
                   dir_new=direction-dirmean;
                   dir_new(dir2)=dir_new(dir2)-360;
                   dirstd=nanstd(dir_new);
               elseif (dirmean>180)
                    dir1=find(direction<=dirmean-180);
                    dir2=find(direction>dirmean-180);
                    dir_new=direction-dirmean;
                    dir_new(dir1)=dir_new(dir1)+360;
                    dirstd=nanstd(dir_new);
               end
               eval([char(eval_string(iiiRCM11(j))) '.dirstd=dirstd;']);


           else
               disp(['File does not exist!'])
               disp(['infile = ' infile])

           end
           j=j+1;

        end
    end

    %--------------------------------------
    % Now read in Argonaut data if required
    %--------------------------------------
    if iiARG>0
        j=1;

        % loop to read one file at a time

        for i=1:length(vecARG);
           serialno = vecARG(i);
           disp('*************************************************************')
           disp(['Reading ARGONAUT - ',num2str(serialno)])
           disp('*************************************************************')

           if isunix
               infile = [procpath,moor,'/arg/',moor,'_',num2str(vecARG(i)),'.use'];
               if exist(infile)==0  % older Arg files had 4 digit serial number starting with zero in filename
                   infile = [procpath,moor,'/arg/',moor,'_0',num2str(vecARG(i)),'.use'];
               end
           elseif ispc
               infile = [procpath,moor,'\arg\',moor,'_',num2str(vecARG(i)),'.use'];
               if exist(infile)==0  % older Arg files had 4 digit serial number starting with zero in filename
                   infile = [procpath,moor,'\arg\',moor,'_0',num2str(vecARG(i)),'.use'];
               end
           end

           % check if file exists
           fileopen=fopen(infile,'r');

           if fileopen>0
               % read data into vectors and then into structure array

               [yy,mm,dd,hh,t,tcat,p,pcat,c,u,v,w,hdg,pit,rol,usd,vsd,wsd,uss,vss,wss,hdgsd,pitsd,rolsd,ipow] = ...
                   rodbload(infile,'yy:mm:dd:hh:t:tcat:p:pcat:c:u:v:w:hdg:pit:rol:usd:vsd:wsd:uss:vss:wss:hdgsd:pitsd:rolsd:ipow');
               jd=julian(yy,mm,dd,hh);

               bad_data=find(p==-9999); p(bad_data)=NaN;
               bad_data=find(t==-9999); t(bad_data)=NaN;
               bad_data=find(u==-9999); u(bad_data)=NaN;
               bad_data=find(v==-9999); v(bad_data)=NaN;


               eval_string(iiiARG(j))={['ARG_' num2str(serialno)]};

               eval([char(eval_string(iiiARG(j))) '.jd=jd;']);
               eval([char(eval_string(iiiARG(j))) '.p=p;']);
               eval([char(eval_string(iiiARG(j))) '.t=t;']);
               eval([char(eval_string(iiiARG(j))) '.u=u;']);
               eval([char(eval_string(iiiARG(j))) '.v=v;']);

               samples=length(jd);
               start_date=gregorian(jd(1));
               end_date=gregorian(jd(end));
               eval([char(eval_string(iiiARG(j))) '.samples=samples;']);
               eval([char(eval_string(iiiARG(j))) '.start_date=start_date;']);
               eval([char(eval_string(iiiARG(j))) '.end_date=end_date;']);

                % calculate speed and direction
               spd = sqrt(u.^2 + v.^2);
               eval([char(eval_string(iiiARG(j))) '.spd=spd;']);
               direction=atan(u./v)*180/pi;
               d = find(v<0);
               direction(d)=direction(d)+180;
               d = find(v>=0 & u<0);
               direction(d)=direction(d)+360;
               eval([char(eval_string(iiiARG(j))) '.dir=direction;']);

               % calculate statistics
               tmean = nanmean(t); pmean = nanmean(p);
               tstd = nanstd(t); pstd = nanstd(p);
               tmax = max(t); pmax = max(p);
               tmin = min(t); pmin = min(p);

               umin = min(u); vmin = min(v);
               umax = max(u); vmax = max(v);
               umean = nanmean(u); vmean = nanmean(v);
               ustd = nanstd(u); vstd = nanstd(v);
               spdmin = min(spd); dirmin = min(direction);
               spdmax = max(spd); dirmax = max(direction);

               eval([char(eval_string(iiiARG(j))) '.pmean=pmean;']);
               eval([char(eval_string(iiiARG(j))) '.tmean=tmean;']);
               eval([char(eval_string(iiiARG(j))) '.umean=umean;']);
               eval([char(eval_string(iiiARG(j))) '.vmean=vmean;']);
               eval([char(eval_string(iiiARG(j))) '.pmin=pmin;']);
               eval([char(eval_string(iiiARG(j))) '.tmin=tmin;']);
               eval([char(eval_string(iiiARG(j))) '.umin=umin;']);
               eval([char(eval_string(iiiARG(j))) '.vmin=vmin;']);
               eval([char(eval_string(iiiARG(j))) '.spdmin=spdmin;']);
               eval([char(eval_string(iiiARG(j))) '.dirmin=dirmin;']);
               eval([char(eval_string(iiiARG(j))) '.pmax=pmax;']);
               eval([char(eval_string(iiiARG(j))) '.tmax=tmax;']);
               eval([char(eval_string(iiiARG(j))) '.umax=umax;']);
               eval([char(eval_string(iiiARG(j))) '.vmax=vmax;']);
               eval([char(eval_string(iiiARG(j))) '.spdmax=spdmax;']);
               eval([char(eval_string(iiiARG(j))) '.dirmax=dirmax;']);
               eval([char(eval_string(iiiARG(j))) '.pstd=pstd;']);
               eval([char(eval_string(iiiARG(j))) '.tstd=tstd;']);
               eval([char(eval_string(iiiARG(j))) '.ustd=ustd;']);
               eval([char(eval_string(iiiARG(j))) '.vstd=vstd;']);


               % calculate spd and dir means and std
               dirmean=atan(umean./vmean)*180/pi;
               if vmean<0
                   dirmean=dirmean+180;
               elseif (vmean>=0 & umean<0);
                   dirmean=dirmean+360;
               end
               eval([char(eval_string(iiiARG(j))) '.dirmean=dirmean;']);

               spdmean=sqrt(umean.^2 + vmean.^2);
               spdstd = nanstd(spd);
               eval([char(eval_string(iiiARG(j))) '.spdmean=spdmean;']);
               eval([char(eval_string(iiiARG(j))) '.spdstd=spdstd;']);

               % for dir STD convert directions to values around mean direction
               % i.e. mean direction becomes 0, and all other directions are relative to
               % that.
               if dirmean<=180
                   dir1=find(direction<=180+dirmean);
                   dir2=find(direction>180+dirmean);
                   dir_new=direction-dirmean;
                   dir_new(dir2)=dir_new(dir2)-360;
                   dirstd=nanstd(dir_new);
               elseif (dirmean>180)
                    dir1=find(direction<=dirmean-180);
                    dir2=find(direction>dirmean-180);
                    dir_new=direction-dirmean;
                    dir_new(dir1)=dir_new(dir1)+360;
                    dirstd=nanstd(dir_new);
               end
               eval([char(eval_string(iiiARG(j))) '.dirstd=dirstd;']);


           else
               disp(['File does not exist!'])
               disp(['infile = ' infile])

           end
           j=j+1;

        end
    end
    %--------------------------------------
    % Now read in Aquadopp data if required
    %--------------------------------------
    if iiNOR>0
        j=1;

        % loop to read one file at a time

        for i=1:length(vecNOR);
           serialno = vecNOR(i);
           disp('*************************************************************')
           disp(['Reading AQUADOPP - ',num2str(serialno)])
           disp('*************************************************************')

           if isunix
               infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.use'];
           elseif ispc
               infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.use'];
           end

           % check if file exists
           fileopen=fopen(infile,'r');

           if fileopen>0
               % read data into vectors and then into structure array

               [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = ...
                   rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd');
               jd=julian(yy,mm,dd,hh);

               bad_data=find(p==-9999); p(bad_data)=NaN;
               bad_data=find(t==-9999); t(bad_data)=NaN;
               bad_data=find(u==-9999); u(bad_data)=NaN;
               bad_data=find(v==-9999); v(bad_data)=NaN;


               eval_string(iiiNOR(j))={['NOR_' num2str(serialno)]};

               eval([char(eval_string(iiiNOR(j))) '.jd=jd;']);
               eval([char(eval_string(iiiNOR(j))) '.p=p;']);
               eval([char(eval_string(iiiNOR(j))) '.t=t;']);
               eval([char(eval_string(iiiNOR(j))) '.u=u;']);
               eval([char(eval_string(iiiNOR(j))) '.v=v;']);

               samples=length(jd);
               start_date=gregorian(jd(1));
               end_date=gregorian(jd(end));
               eval([char(eval_string(iiiNOR(j))) '.samples=samples;']);
               eval([char(eval_string(iiiNOR(j))) '.start_date=start_date;']);
               eval([char(eval_string(iiiNOR(j))) '.end_date=end_date;']);

               % calculate speed and direction 
               spd = sqrt(u.^2 + v.^2);
               eval([char(eval_string(iiiNOR(j))) '.spd=spd;']);
               direction=atan(u./v)*180/pi;
               d = find(v<0);
               direction(d)=direction(d)+180;
               d = find(v>=0 & u<0);
               direction(d)=direction(d)+360;
               eval([char(eval_string(iiiNOR(j))) '.dir=direction;']);

               % calculate statistics
               tmean = nanmean(t); pmean = nanmean(p);
               tstd = nanstd(t); pstd = nanstd(p);
               tmax = max(t); pmax = max(p);
               tmin = min(t); pmin = min(p);

               umin = min(u); vmin = min(v);
               umax = max(u); vmax = max(v);
               umean = nanmean(u); vmean = nanmean(v);
               ustd = nanstd(u); vstd = nanstd(v);
               spdmin = min(spd); dirmin = min(direction);
               spdmax = max(spd); dirmax = max(direction);

               eval([char(eval_string(iiiNOR(j))) '.pmean=pmean;']);
               eval([char(eval_string(iiiNOR(j))) '.tmean=tmean;']);
               eval([char(eval_string(iiiNOR(j))) '.umean=umean;']);
               eval([char(eval_string(iiiNOR(j))) '.vmean=vmean;']);
               eval([char(eval_string(iiiNOR(j))) '.pmin=pmin;']);
               eval([char(eval_string(iiiNOR(j))) '.tmin=tmin;']);
               eval([char(eval_string(iiiNOR(j))) '.umin=umin;']);
               eval([char(eval_string(iiiNOR(j))) '.vmin=vmin;']);
               eval([char(eval_string(iiiNOR(j))) '.spdmin=spdmin;']);
               eval([char(eval_string(iiiNOR(j))) '.dirmin=dirmin;']);
               eval([char(eval_string(iiiNOR(j))) '.pmax=pmax;']);
               eval([char(eval_string(iiiNOR(j))) '.tmax=tmax;']);
               eval([char(eval_string(iiiNOR(j))) '.umax=umax;']);
               eval([char(eval_string(iiiNOR(j))) '.vmax=vmax;']);
               eval([char(eval_string(iiiNOR(j))) '.spdmax=spdmax;']);
               eval([char(eval_string(iiiNOR(j))) '.dirmax=dirmax;']);
               eval([char(eval_string(iiiNOR(j))) '.pstd=pstd;']);
               eval([char(eval_string(iiiNOR(j))) '.tstd=tstd;']);
               eval([char(eval_string(iiiNOR(j))) '.ustd=ustd;']);
               eval([char(eval_string(iiiNOR(j))) '.vstd=vstd;']);


               % calculate spd and dir means and std
               dirmean=atan(umean./vmean)*180/pi;
               if vmean<0
                   dirmean=dirmean+180;
               elseif (vmean>=0 & umean<0);
                   dirmean=dirmean+360;
               end
               eval([char(eval_string(iiiNOR(j))) '.dirmean=dirmean;']);

               spdmean=sqrt(umean.^2 + vmean.^2);
               spdstd = nanstd(spd);
               eval([char(eval_string(iiiNOR(j))) '.spdmean=spdmean;']);
               eval([char(eval_string(iiiNOR(j))) '.spdstd=spdstd;']);

               % for dir STD convert directions to values around mean direction
               % i.e. mean direction becomes 0, and all other directions are relative to
               % that.
               if dirmean<=180
                   dir1=find(direction<=180+dirmean);
                   dir2=find(direction>180+dirmean);
                   dir_new=direction-dirmean;
                   dir_new(dir2)=dir_new(dir2)-360;
                   dirstd=nanstd(dir_new);
               elseif (dirmean>180)
                    dir1=find(direction<=dirmean-180);
                    dir2=find(direction>dirmean-180);
                    dir_new=direction-dirmean;
                    dir_new(dir1)=dir_new(dir1)+360;
                    dirstd=nanstd(dir_new);
               end
               eval([char(eval_string(iiiNOR(j))) '.dirstd=dirstd;']);

           else
               disp(['File does not exist!'])
               disp(['infile = ' infile])

           end
           j=j+1;

        end
    end
    %--------------------------------------
    % Now read in DVS data if required
    %--------------------------------------
    if iiDVS>0
        j=1;

        % loop to read one file at a time

        for i=1:length(vecDVS);
           serialno = vecDVS(i);
           disp('*************************************************************')
           disp(['Reading DVS - ',num2str(serialno)])
           disp('currently only reading Bin 2 data - will need to update code!')
           disp('*************************************************************')

           if isunix
               infile = [procpath,moor,'/dvs/',moor,'_',sprintf('%5.5d',vecDVS(i)),'_bin2.use'];
           elseif ispc
               infile = [procpath,moor,'\dvs\',moor,'_',sprintf('%5.5d',vecDVS(i)),'_bin2.use'];
           end

           % check if file exists
           fileopen=fopen(infile,'r');

           if fileopen>0
               % read data into vectors and then into structure array

               [yy,mm,dd,hh,t,u,v,w,hdg,pit,rol,cs,cd] = rodbload(infile,'yy:mm:dd:hh:t:u:v:w:hdg:pit:rol:cs:cd');
               jd=julian(yy,mm,dd,hh);

               bad_data=find(t==-9999); t(bad_data)=NaN;
               bad_data=find(u==-9999); u(bad_data)=NaN;
               bad_data=find(v==-9999); v(bad_data)=NaN;



               eval_string(iiiDVS(j))={['DVS_' num2str(serialno)]};

               eval([char(eval_string(iiiDVS(j))) '.jd=jd;']);
               eval([char(eval_string(iiiDVS(j))) '.t=t;']);
               eval([char(eval_string(iiiDVS(j))) '.u=u;']);
               eval([char(eval_string(iiiDVS(j))) '.v=v;']);


               samples=length(jd);
               start_date=gregorian(jd(1));
               end_date=gregorian(jd(end));
               eval([char(eval_string(iiiDVS(j))) '.samples=samples;']);
               eval([char(eval_string(iiiDVS(j))) '.start_date=start_date;']);
               eval([char(eval_string(iiiDVS(j))) '.end_date=end_date;']);

               % calculate speed and direction 
               spd = sqrt(u.^2 + v.^2);
               eval([char(eval_string(iiiDVS(j))) '.spd=spd;']);
               direction=atan(u./v)*180/pi;
               d = find(v<0);
               direction(d)=direction(d)+180;
               d = find(v>=0 & u<0);
               direction(d)=direction(d)+360;
               eval([char(eval_string(iiiDVS(j))) '.dir=direction;']);

               % calculate statistics
               tmean = nanmean(t);
               tstd = nanstd(t); 
               tmax = max(t); 
               tmin = min(t); 

               umin = min(u); vmin = min(v);
               umax = max(u); vmax = max(v);
               umean = nanmean(u); vmean = nanmean(v);
               ustd = nanstd(u); vstd = nanstd(v);
               spdmin = min(spd); dirmin = min(direction);
               spdmax = max(spd); dirmax = max(direction);


               eval([char(eval_string(iiiDVS(j))) '.tmean=tmean;']);
               eval([char(eval_string(iiiDVS(j))) '.umean=umean;']);
               eval([char(eval_string(iiiDVS(j))) '.vmean=vmean;']);
               eval([char(eval_string(iiiDVS(j))) '.tmin=tmin;']);
               eval([char(eval_string(iiiDVS(j))) '.umin=umin;']);
               eval([char(eval_string(iiiDVS(j))) '.vmin=vmin;']);
               eval([char(eval_string(iiiDVS(j))) '.spdmin=spdmin;']);
               eval([char(eval_string(iiiDVS(j))) '.dirmin=dirmin;']);
               eval([char(eval_string(iiiDVS(j))) '.tmax=tmax;']);
               eval([char(eval_string(iiiDVS(j))) '.umax=umax;']);
               eval([char(eval_string(iiiDVS(j))) '.vmax=vmax;']);
               eval([char(eval_string(iiiDVS(j))) '.spdmax=spdmax;']);
               eval([char(eval_string(iiiDVS(j))) '.dirmax=dirmax;']);
               eval([char(eval_string(iiiDVS(j))) '.tstd=tstd;']);
               eval([char(eval_string(iiiDVS(j))) '.ustd=ustd;']);
               eval([char(eval_string(iiiDVS(j))) '.vstd=vstd;']);


               % calculate spd and dir means and std
               dirmean=atan(umean./vmean)*180/pi;
               if vmean<0
                   dirmean=dirmean+180;
               elseif (vmean>=0 & umean<0);
                   dirmean=dirmean+360;
               end
               eval([char(eval_string(iiiDVS(j))) '.dirmean=dirmean;']);

               spdmean=sqrt(umean.^2 + vmean.^2);
               spdstd = nanstd(spd);
               eval([char(eval_string(iiiDVS(j))) '.spdmean=spdmean;']);
               eval([char(eval_string(iiiDVS(j))) '.spdstd=spdstd;']);

               % for dir STD convert directions to values around mean direction
               % i.e. mean direction becomes 0, and all other directions are relative to
               % that.
               if dirmean<=180
                   dir1=find(direction<=180+dirmean);
                   dir2=find(direction>180+dirmean);
                   dir_new=direction-dirmean;
                   dir_new(dir2)=dir_new(dir2)-360;
                   dirstd=nanstd(dir_new);
               elseif (dirmean>180)
                    dir1=find(direction<=dirmean-180);
                    dir2=find(direction>dirmean-180);
                    dir_new=direction-dirmean;
                    dir_new(dir1)=dir_new(dir1)+360;
                    dirstd=nanstd(dir_new);
               end
               eval([char(eval_string(iiiDVS(j))) '.dirstd=dirstd;']);

           else
               disp(['File does not exist!'])
               disp(['infile = ' infile])

           end
           j=j+1;

        end
    end
    %---------------------------------------------------
    % Now read in Seabird and Ixsea BPR data if required
    %---------------------------------------------------
    if iiBPR>0
        j=1;

        % loop to read one file at a time

        for i=1:length(vecBPR);
           serialno = vecBPR(i);
           disp('*************************************************************')
           disp(['Reading BPR - ',num2str(serialno)])
           disp('*************************************************************')

           if id(iiBPR)==465
               if isunix
                   infile = [procpath,moor,'/seagauge/',moor,'_',sprintf('%04d',vecBPR(i)),'.use'];
               elseif ispc
                   infile = [procpath,moor,'\seagauge\',moor,'_',sprintf('%04d',vecBPR(i)),'.use'];
               end

           elseif id(iiBPR)==470
               if isunix
                   infile = [procpath,moor,'/ixsbpr/',moor,'_',sprintf('%04d',vecBPR(i)),'.use'];
               elseif ispc
                   infile = [procpath,moor,'\ixsbpr\',moor,'_',sprintf('%04d',vecBPR(i)),'.use'];
               end
           end
           % check if file exists
           fileopen=fopen(infile,'r');

           if fileopen>0
               % read data into vectors and then into structure array

               [yy,mm,dd,hh,t,p] = ...
                   rodbload(infile,'yy:mm:dd:hh:t:p');
               jd=julian(yy,mm,dd,hh);

               bad_data=find(p==-9999); p(bad_data)=NaN;
               bad_data=find(t==-9999); t(bad_data)=NaN;

               eval_string(iiiBPR(j))={['BPR_' num2str(serialno)]};

               eval([char(eval_string(iiiBPR(j))) '.jd=jd;']);
               eval([char(eval_string(iiiBPR(j))) '.p=p;']);
               eval([char(eval_string(iiiBPR(j))) '.t=t;']);

               samples=length(jd);
               start_date=gregorian(jd(1));
               end_date=gregorian(jd(end));
               eval([char(eval_string(iiiBPR(j))) '.samples=samples;']);
               eval([char(eval_string(iiiBPR(j))) '.start_date=start_date;']);
               eval([char(eval_string(iiiBPR(j))) '.end_date=end_date;']);

               % calculate statistics
               tmean = nanmean(t); pmean = nanmean(p);
               tstd= nanstd(t); pstd= nanstd(p);
               tmax = max(t); pmax = max(p);
               tmin = min(t); pmin = min(p);

               eval([char(eval_string(iiiBPR(j))) '.pmean=pmean;']);
               eval([char(eval_string(iiiBPR(j))) '.tmean=tmean;']);
               eval([char(eval_string(iiiBPR(j))) '.pmin=pmin;']);
               eval([char(eval_string(iiiBPR(j))) '.tmin=tmin;']);
               eval([char(eval_string(iiiBPR(j))) '.pmax=pmax;']);
               eval([char(eval_string(iiiBPR(j))) '.tmax=tmax;']);
               eval([char(eval_string(iiiBPR(j))) '.pstd=pstd;']);
               eval([char(eval_string(iiiBPR(j))) '.tstd=tstd;']);
           else
               disp(['File does not exist!'])
               disp(['infile = ' infile])

           end
           j=j+1;
        end

    end

    % ============================================
    % START OF OUPUTTING SUMMARY DATA TO TEXT FILE
    % ============================================
    outfile=[outpath,moor '_stats.asc'];
    check=0;
    while check==0
        if exist(outfile,'file')
            disp(['Outfile ' outfile ' already exists!'])
            overwrite=input('Do you wish to overwrite the file? y/n (default=n):- ','s');
            if overwrite~='y';
                outfile=input('Please enter alternative outfile name:- ','s');
            else
                check=1;
            end
        else
            check=1;
        end
    end    
    fid=fopen(outfile,'w');

    fprintf(fid,'%s \n%s %s \n','OSNAP Mooring Array.','Simple Statisctics for Mooring:- ',moor);
    fprintf(fid,'%s %0.2i/%0.2i/%i %0.2i:%0.2i\n','Mooring deployment - start: ',s_d(3),s_d(2),s_d(1),s_t(1),s_t(2));
    fprintf(fid,'%s %0.2i/%0.2i/%i %0.2i:%0.2i\n\n','                       end: ',e_d(3),e_d(2),e_d(1),e_t(1),e_t(2));
    fprintf(fid,'%s\n','-----------------------------------------------------------------------------------------');
    fprintf(fid,'%s\n','       SN  var       first            last        valid    mean   stdev     min     max');
    fprintf(fid,'%s\n','                    record          record      records');
    fprintf(fid,'%s\n','-----------------------------------------------------------------------------------------');

    % remove double occurence of Sonteks with SMPs
    % and remove MMPs as this routine is too simple for them
    % also remove releases as some info.dat files have them in

    %a=find((depths(:,1)~=366337)&(depths(:,1)~=380)&(depths(:,1)~=311));
    sn2=depths(:,3);



    % NOTE this routine assumes that all current meters have pressure and
    % temperature sensors, and that the RCM11 and S4 also have conductivity
    % this is likely to change in future as new CMs are unlikely to have cond
    % and possibly not temp.

    CTD=[iiiMC;iiiRBR;iiiIDR];
    CMCTD=[iiiRCM11;iiiS4];
    CM=[iiiARG;iiiNOR];
    BPR=[iiiBPR];
    CM_NO_P=[iiiDVS];

    ctd_or_cm=zeros(size(sn2));
    ctd_or_cm(CTD)=1;
    ctd_or_cm(CMCTD)=2;
    ctd_or_cm(CM)=3;
    ctd_or_cm(BPR)=4;
    ctd_or_cm(CM_NO_P)=5;

    for i=1:size(sn2)
        variables=[' p ';' t ';' c ';' u ';' v ';'spd';'dir']; % pad single character variables with spaces
        if ctd_or_cm(i)==1 % i.e. is a ctd
            stats_rows=3;
            variables=variables(1:3,:);
        elseif ctd_or_cm(i)==2 % i.e. is a cm with ctd
            stats_rows=7;
            %variables unchanged
        elseif ctd_or_cm(i)==3 % i.e. is a cm with t and p but not c (NOR or ARG)
            stats_rows=6;
            variables=[variables(1:2,:);variables(4:7,:)]; % remove c variable
        elseif ctd_or_cm(i)==4 % i.e. is a SBE BPR
            stats_rows=2;
            variables=variables(1:2,:);
        elseif ctd_or_cm(i)==5 % current meter with t but not c or p
            stats_rows=5;
            variables=[variables(2,:);variables(4:7,:)]; % remove c and p variables
        end
        clear SN
        for k=1:stats_rows
            SN{k}=' ';
        end
        SN{1}=num2str(sn2(i));
        if find(sn2(i)==vecMC)
            eval_string=['MC_' num2str(sn2(i))];
        elseif find(sn2(i)==vecRCM11)
            eval_string=['RCM11_' num2str(sn2(i))];
        elseif find(sn2(i)==vecARG)
            eval_string=['ARG_' num2str(sn2(i))];
        elseif find(sn2(i)==vecIDR)
            eval_string=['IDR_' num2str(sn2(i))];
        elseif find(sn2(i)==vecNOR)
            eval_string=['NOR_' num2str(sn2(i))];
        elseif find(sn2(i)==vecS4)
            eval_string=['S4_' num2str(sn2(i))];
        elseif find(sn2(i)==vecDVS)
            eval_string=['DVS_' num2str(sn2(i))];
        elseif find(sn2(i)==vecRBR)
            eval_string=['RBR_' num2str(sn2(i))];
        elseif find(sn2(i)==vecBPR)
            eval_string=['BPR_' num2str(sn2(i))];
        end



        % Use automatic detection of first and last useable record
        % NB: this works on the assumption that bad data has been replaced with
        % -9999 values - which are then replaced with NaNs during loading
        j2=size(variables);
        j2=j2(1);
        for j=1:j2
            eval(['validRecs=find(~isnan(' char(eval_string) '.' strtrim(variables(j,:)) '));']);
            if length(validRecs) > 0
              firstRec(j)=validRecs(1);
              lastRec(j)=validRecs(end);
              numRecs(j)=length(validRecs);

              eval(['firstRecJD=' char(eval_string) '.jd(firstRec(j));']);
              eval(['lastRecJD=' char(eval_string) '.jd(lastRec(j));']);
              firstRecGREG(j,:)=gregorian(firstRecJD);
              lastRecGREG(j,:)=gregorian(lastRecJD);

            % combine mean, max, min and std into columns for all variables per
            % instrument
              eval(['meanRec(j)=' char(eval_string) '.' strtrim(variables(j,:)) 'mean;']);
              eval(['stdRec(j)=' char(eval_string) '.' strtrim(variables(j,:)) 'std;']);
              eval(['minRec(j)=' char(eval_string) '.' strtrim(variables(j,:)) 'min;']);
              eval(['maxRec(j)=' char(eval_string) '.' strtrim(variables(j,:)) 'max;']);

              short_year1=num2str(firstRecGREG(j,1));
              short_year2=num2str(lastRecGREG(j,1));
              short_year1=short_year1(3:4);
              short_year2=short_year2(3:4);
            % write data to file
              fprintf(fid,' %8s  %s  %02.0f/%02.0f/%s %02.0f:%02.0f  %02.0f/%02.0f/%s %02.0f:%02.0f  %7.0f  %6.1f  %6.1f  %6.1f  %6.1f\n',...
              char(SN{j}),variables(j,:),firstRecGREG(j,3),firstRecGREG(j,2), ...
              short_year1,firstRecGREG(j,4),firstRecGREG(j,5),...
              lastRecGREG(j,3),lastRecGREG(j,2),short_year2,lastRecGREG(j,4),lastRecGREG(j,5),...
              numRecs(j),meanRec(j),stdRec(j),minRec(j),maxRec(j));
           else
              fprintf(fid,' %8s  %s  No valid data \n',char(SN{j}),variables(j,:));
           end
        end




        % print dividing line between instruments
        fprintf(fid,'%s\n','-----------------------------------------------------------------------------------------');
    end
    fclose(fid);
end

if non_verbose==1
    a=who('MC*');
    for i=1:length(a)
        eval([a{i} '.sn=a{i}(4:end)']);
        deletefields={'p','t','c','jd'};
        
        eval(['rmfield(' a{i} ',deletefields)'])
        %for j=1:length(deletefields)
        %    eval(['rmfield(' a{i} ',deletefields)'])
        %end
        
        data(i)=eval(['{' char(a(i)) '}']);
                
    end
    varargout{1}=data;
    %keyboard
end
