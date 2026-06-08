function vel_vs_temp_corr(moor,varargin)
%
% Function for calculating and plotting the correlation betwen northward 
% velocity and temperautre. - Suggested by Harry for checking if any
% correlation for heat transport estimates
%
% NB: currently using temp and not potemp.
%
% function vel_vs_temp_corr('moor','layout','procpath')
%
% required inputs:-
%   moor: complete mooring name as string. e.g. 'wb1_1_200420'
%
% optional inputs:-
%   layout: orientation of figure portrait/lanscape (default = portrait)
%           input of 1 = landscape, 0 = portrait
%   procpath: can specify exact path to proc directory if not using 
%             standard data paths. 
%             e.g. '/Volumes/noc/mpoc/rpdmoc/rapid/data/moor/proc/'
%   unfiltered:  plots data in unfiltered format - useful for more detailed
%                inspection
%
% functions called:-
%   rodbload, julian, auto_filt
%   from .../exec/moor/tools and .../exec/moor/rodb paths
% 
% Routine written by Darren Rayner Feb 2010 - adapted from currents_overlay.m.


if nargin <1
    help vel_vs_temp_corr
    return
end

varargin_string=varargin;
for i=1:length(varargin) % need to change numeric values in varargin to string so can search with strmatch below
    if isnumeric(varargin{i})
        varargin_string{i}=num2str(varargin{i});
    end
end

% check for optional arguments
a=strmatch('layout',varargin_string,'exact');
if a>0
    layout=str2num(varargin_string{a+1});
else
    layout='portrait';
end

a=strmatch('procpath',varargin_string,'exact');
if a>0
    procpath=char(varargin{a+1});
else
    data_report_tools_dir=which('data_report_tools');
    b=strfind(data_report_tools_dir,'/');
    data_report_tools_dir=data_report_tools_dir(1:b(end));
    procpath=[data_report_tools_dir '../../../moor/proc/']; %DR changed to be relative paths now that data_report_tools are on the network 19/2/12
end

a=strmatch('unfiltered',varargin_string,'exact');
if a>0
    filt_or_unfilt=1;
else
    filt_or_unfilt=0; 
end


if isunix
    infofile=[procpath,moor,'/',moor,'info.dat'];
elseif ispc
    infofile=[procpath,moor,'\',moor,'info.dat'];
end


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
jd_start = julian([s_d' hms2h([s_t;0]')']);
jd_end   = julian([e_d' hms2h([e_t;0]')']);

disp(['z : instrument id : serial number'])
for i = 1:length(id)
    disp([z(i),id(i),sn(i)])
end





% find the index number of S4s
iiS4 = find(id == 302);
vecS4 = sn(iiS4);
% and find index number of RCM11s
iiRCM11 = find(id == 310);
vecRCM11 = sn(iiRCM11);
% and find index number of Sontek Argonauts
iiARG = find(id == 366);
vecARG = sn(iiARG);
% and find index number of Norteks
iiNOR = find((id == 368|370));
vecNOR = sn(iiNOR);
% and find index number of MicroCATs - important if have near a current meter.
% Can use more accurate MicroCAT rather than current
% meter temp sensor
iiMC = find(id == 337);
vecMC = sn(iiMC);

depths(:,1) = id([iiS4;iiRCM11;iiARG;iiNOR;iiMC]);
depths(:,2) = z([iiS4;iiRCM11;iiARG;iiNOR;iiMC]);
depths=sortrows(depths,2);
iiiS4=find(depths(:,1)==302);  
iiiRCM11=find(depths(:,1)==310);  
iiiARG=find(depths(:,1)==366 | depths(:,1)==366337);  
iiiNOR=find(depths(:,1)==368|depths(:,1)==370);
iiiMC=find(depths(:,1)==337);
iii=[iiiS4;iiiRCM11;iiiARG;iiiNOR]; %Microcats not included

%% set up figure window
% set figure size on screen for better viewing
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

%set print area of figure
pos1  = [1/8*scnsize(3),8*bdwidth,1/2*scnsize(3),(scnsize(4) - 30*bdwidth)];
correlation_plot=figure('Position',pos1);
    
set(correlation_plot,'PaperUnits','centimeters');
set(correlation_plot, 'PaperType', 'A4');
set(correlation_plot, 'PaperOrientation',layout);
papersize = get(correlation_plot,'PaperSize');
width=17; height=26; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
figuresize = [left, bottom, width, height];
set(correlation_plot, 'PaperPosition', figuresize);

%%
plot_string={};

MAX_v=0;
MAX_t=0;

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

       % read data into vectors and then into structure array

       [yy,mm,dd,hh,u,v,t,c,p,hdg] = rodbload(infile,'yy:mm:dd:hh:u:v:t:c:p:hdg');
       jd=julian(yy,mm,dd,hh);

       bad_data=find(u==-9999); u(bad_data)=NaN;
       bad_data=find(v==-9999); v(bad_data)=NaN;
       bad_data=find(t==-9999); t(bad_data)=NaN;
       bad_data=find(c==-9999); c(bad_data)=NaN;
       bad_data=find(p==-9999); p(bad_data)=NaN;
       bad_data=find(hdg==-9999); hdg(bad_data)=NaN;

       eval_string(iiiS4(j))={['S4_' num2str(serialno)]};
       
       eval([char(eval_string(iiiS4(j))) '.jd=jd;']);
       eval([char(eval_string(iiiS4(j))) '.u=u;']);
       eval([char(eval_string(iiiS4(j))) '.v=v;']);
       eval([char(eval_string(iiiS4(j))) '.t=t;']);
       eval([char(eval_string(iiiS4(j))) '.c=c;']);
       eval([char(eval_string(iiiS4(j))) '.p=p;']);
       eval([char(eval_string(iiiS4(j))) '.hdg=hdg;']);

       sampling_rate = 1/median(diff(jd));

       if filt_or_unfilt==0 % i.e. want to filter the data prior to plotting
           % Apply a butterworth filter to the data using auto_filt and use for
           % plots
           
           ii = eval(['find(~isnan(' char(eval_string(iiiS4(j))) '.u));']); 
           eval([char(eval_string(iiiS4(j))) '.u(ii)=auto_filt(' char(eval_string(iiiS4(j)))...
                 '.u(ii), sampling_rate, 1/2,''low'',4);']);
           ii = eval(['find(~isnan(' char(eval_string(iiiS4(j))) '.v));']); 
           eval([char(eval_string(iiiS4(j))) '.v(ii)=auto_filt(' char(eval_string(iiiS4(j)))...
                 '.v(ii), sampling_rate, 1/2,''low'',4);']);
           ii = eval(['find(~isnan(' char(eval_string(iiiS4(j))) '.t));']);
           eval([char(eval_string(iiiS4(j))) '.t(ii)=auto_filt(' char(eval_string(iiiS4(j)))...
                 '.t(ii), sampling_rate, 1/2,''low'',4);']);
       end

       % calculate spd and dir from u and v - may not be needed but don't
       % delete yet
       current_speed = eval(['sqrt((' char(eval_string(iiiS4(j))) '.u).^2 + (' char(eval_string(iiiS4(j))) '.v) .^2);']);
       eval([char(eval_string(iiiS4(j))) '.spd=current_speed;']);
       current_direction=eval(['atan(' char(eval_string(iiiS4(j))) '.u ./' char(eval_string(iiiS4(j))) '.v)*180/pi;']);
       d = eval(['find(' char(eval_string(iiiS4(j))) '.v<0);']);
       current_direction(d)=current_direction(d)+180;
       d = eval(['find(' char(eval_string(iiiS4(j))) '.v>=0 & ' char(eval_string(iiiS4(j))) '.u<0);']);
       current_direction(d)=current_direction(d)+360;
       eval([char(eval_string(iiiS4(j))) '.dir=current_direction;']);

       
       % determine max v and max t for setting axes after plotting
       max_v=max(abs(v));
       if max_v > MAX_v
           MAX_v = max_v;
       end
       max_t=max(t);
       if max_t > MAX_t
           MAX_t = max_t;
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
       

       % read data into vectors and then into structure array

       [yy,mm,dd,hh,ref,u,v,t,c,p,tlt,mss] = rodbload(infile,'yy:mm:dd:hh:ref:u:v:t:c:p:tlt:mss');
       jd=julian(yy,mm,dd,hh);

       bad_data=find(t==-9999); t(bad_data)=NaN;
       bad_data=find(p==-9999); p(bad_data)=NaN;
       bad_data=find(u==-9999); u(bad_data)=NaN;
       bad_data=find(v==-9999); v(bad_data)=NaN;
       bad_data=find(c==-9999); c(bad_data)=NaN;
       bad_data=find(tlt==-9999); tlt(bad_data)=NaN;
       bad_data=find(mss==-9999); mss(bad_data)=NaN;

       eval_string(iiiRCM11(j))={['RCM11_' num2str(serialno)]};
       
       eval([char(eval_string(iiiRCM11(j))) '.jd=jd;']);
       eval([char(eval_string(iiiRCM11(j))) '.t=t;']);
       eval([char(eval_string(iiiRCM11(j))) '.p=p;']);
       eval([char(eval_string(iiiRCM11(j))) '.u=u;']);
       eval([char(eval_string(iiiRCM11(j))) '.v=v;']);
       eval([char(eval_string(iiiRCM11(j))) '.c=c;']);
       eval([char(eval_string(iiiRCM11(j))) '.tlt=tlt;']);
       eval([char(eval_string(iiiRCM11(j))) '.mss=mss;']);

       sampling_rate = 1/median(diff(jd));

       if filt_or_unfilt==0 % i.e. want to filter the data prior to plotting
           % Apply a butterworth filter to the data using auto_filt and use for
           % plots
           
           ii = eval(['find(~isnan(' char(eval_string(iiiRCM11(j))) '.u));']);
           eval([char(eval_string(iiiRCM11(j))) '.u(ii)=auto_filt(' char(eval_string(iiiRCM11(j)))...
                 '.u(ii), sampling_rate, 1/2,''low'',4);']);
           ii = eval(['find(~isnan(' char(eval_string(iiiRCM11(j))) '.v));']); 
           eval([char(eval_string(iiiRCM11(j))) '.v(ii)=auto_filt(' char(eval_string(iiiRCM11(j)))...
                 '.v(ii), sampling_rate, 1/2,''low'',4);']);
           ii = eval(['find(~isnan(' char(eval_string(iiiRCM11(j))) '.t));']); 
           eval([char(eval_string(iiiRCM11(j))) '.t(ii)=auto_filt(' char(eval_string(iiiRCM11(j)))...
                 '.t(ii), sampling_rate, 1/2,''low'',4);']);
       end
       
       current_speed = eval(['sqrt((' char(eval_string(iiiRCM11(j))) '.u).^2 + (' char(eval_string(iiiRCM11(j))) '.v).^2);']);
       eval([char(eval_string(iiiRCM11(j))) '.spd=current_speed;']);
       current_direction=eval(['atan(' char(eval_string(iiiRCM11(j))) '.u ./' char(eval_string(iiiRCM11(j))) '.v)*180/pi;']);
       d = eval(['find(' char(eval_string(iiiRCM11(j))) '.v<0);']);
       current_direction(d)=current_direction(d)+180;
       d = eval(['find(' char(eval_string(iiiRCM11(j))) '.v>=0 & ' char(eval_string(iiiRCM11(j))) '.u<0);']);
       current_direction(d)=current_direction(d)+360;
       eval([char(eval_string(iiiRCM11(j))) '.dir=current_direction;']);
       
       % determine max v and t for setting axes after plotting
       max_v=max(abs(v));
       if max_v > MAX_v
           MAX_v = max_v;
       end
       max_t=max(abs(t));
       if max_t > MAX_t
           MAX_t = max_t;
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

       % read data into vectors and then into structure array

       [yy,mm,dd,hh,t,tcat,p,pcat,c,u,v,w,hdg,pit,rol,usd,vsd,wsd,uss,vss,wss,hdgsd,pitsd,rolsd,ipow] = ...
           rodbload(infile,'yy:mm:dd:hh:t:tcat:p:pcat:c:u:v:w:hdg:pit:rol:usd:vsd:wsd:uss:vss:wss:hdgsd:pitsd:rolsd:ipow');
       jd=julian(yy,mm,dd,hh);

       bad_data=find(t==-9999); t(bad_data)=NaN;
       bad_data=find(p==-9999); p(bad_data)=NaN;
       bad_data=find(u==-9999); u(bad_data)=NaN;
       bad_data=find(v==-9999); v(bad_data)=NaN;
       bad_data=find(w==-9999); w(bad_data)=NaN;
       bad_data=find(hdg==-9999); hdg(bad_data)=NaN;
       bad_data=find(pit==-9999); pit(bad_data)=NaN;
       bad_data=find(rol==-9999); rol(bad_data)=NaN;
       bad_data=find(usd==-9999); usd(bad_data)=NaN;
       bad_data=find(vsd==-9999); vsd(bad_data)=NaN;
       bad_data=find(wsd==-9999); wsd(bad_data)=NaN;
       bad_data=find(uss==-9999); uss(bad_data)=NaN;
       bad_data=find(vss==-9999); vss(bad_data)=NaN;
       bad_data=find(wss==-9999); wss(bad_data)=NaN;
       bad_data=find(hdgsd==-9999); hdgsd(bad_data)=NaN;
       bad_data=find(pitsd==-9999); pitsd(bad_data)=NaN;
       bad_data=find(rolsd==-9999); rolsd(bad_data)=NaN;
       bad_data=find(ipow==-9999); ipow(bad_data)=NaN;
       
       eval_string(iiiARG(j))={['ARG_' num2str(serialno)]};
       
       eval([char(eval_string(iiiARG(j))) '.jd=jd;']);
       eval([char(eval_string(iiiARG(j))) '.t=t;']);
       eval([char(eval_string(iiiARG(j))) '.p=p;']);
       eval([char(eval_string(iiiARG(j))) '.u=u;']);
       eval([char(eval_string(iiiARG(j))) '.v=v;']);
       eval([char(eval_string(iiiARG(j))) '.w=w;']);
       eval([char(eval_string(iiiARG(j))) '.hdg=hdg;']);
       eval([char(eval_string(iiiARG(j))) '.pit=pit;']);
       eval([char(eval_string(iiiARG(j))) '.rol=rol;']);
       eval([char(eval_string(iiiARG(j))) '.usd=usd;']);
       eval([char(eval_string(iiiARG(j))) '.vsd=vsd;']);
       eval([char(eval_string(iiiARG(j))) '.wsd=wsd;']);
       eval([char(eval_string(iiiARG(j))) '.uss=uss;']);
       eval([char(eval_string(iiiARG(j))) '.vss=vss;']);
       eval([char(eval_string(iiiARG(j))) '.wss=wss;']);
       eval([char(eval_string(iiiARG(j))) '.hdgsd=hdgsd;']);
       eval([char(eval_string(iiiARG(j))) '.pitsd=pitsd;']);
       eval([char(eval_string(iiiARG(j))) '.rolsd=rolsd;']);
       eval([char(eval_string(iiiARG(j))) '.ipow=ipow;']);
       
       sampling_rate = 1/median(diff(jd));

       if filt_or_unfilt==0 % i.e. want to filter the data prior to plotting
           % Apply a butterworth filter to the data using auto_filt and use for
           % plots
           
           ii = eval(['find(~isnan(' char(eval_string(iiiARG(j))) '.u));']); 
           eval([char(eval_string(iiiARG(j))) '.u(ii)=auto_filt(' char(eval_string(iiiARG(j)))...
                 '.u(ii), sampling_rate, 1/2,''low'',4);']);
           ii = eval(['find(~isnan(' char(eval_string(iiiARG(j))) '.v));']); 
           eval([char(eval_string(iiiARG(j))) '.v(ii)=auto_filt(' char(eval_string(iiiARG(j)))...
                 '.v(ii), sampling_rate, 1/2,''low'',4);']);
           ii = eval(['find(~isnan(' char(eval_string(iiiARG(j))) '.t));']); 
           eval([char(eval_string(iiiARG(j))) '.t(ii)=auto_filt(' char(eval_string(iiiARG(j)))...
                 '.t(ii), sampling_rate, 1/2,''low'',4);']);
       end
       
       current_speed = eval(['sqrt((' char(eval_string(iiiARG(j))) '.u).^2 + (' char(eval_string(iiiARG(j))) '.v).^2);']);
       eval([char(eval_string(iiiARG(j))) '.spd=current_speed;']);
       current_direction=eval(['atan(' char(eval_string(iiiARG(j))) '.u ./' char(eval_string(iiiARG(j))) '.v)*180/pi;']);
       d = eval(['find(' char(eval_string(iiiARG(j))) '.v<0);']);
       current_direction(d)=current_direction(d)+180;
       d = eval(['find(' char(eval_string(iiiARG(j))) '.v>=0 & ' char(eval_string(iiiARG(j))) '.u<0);']);
       current_direction(d)=current_direction(d)+360;
       eval([char(eval_string(iiiARG(j))) '.dir=current_direction;']);
       
       % determine current_speed for setting y-axis after plotting
       max_v=max(abs(v));
       if max_v > MAX_v
           MAX_v = max_v;
       end
       if max_t > MAX_t
           MAX_t = max_t;
       end
       j=j+1;
    end
end

%------------------------------------
% Now read in NORTEK data if required
%------------------------------------
if iiNOR>0
    j=1;
    
    % loop to read one file at a time

    for i=1:length(vecNOR);
       serialno = vecNOR(i);
       disp('*************************************************************')
       disp(['Reading NORTEK - ',num2str(serialno)])
       disp('*************************************************************')

       if isunix
           infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.use'];
       elseif ispc
           infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.use'];
       end
       

       % read data into vectors and then into structure array

       [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd');
       jd=julian(yy,mm,dd,hh);

       bad_data=find(t==-9999); t(bad_data)=NaN;
       bad_data=find(p==-9999); p(bad_data)=NaN;
       bad_data=find(u==-9999); u(bad_data)=NaN;
       bad_data=find(v==-9999); v(bad_data)=NaN;
       bad_data=find(w==-9999); w(bad_data)=NaN;
       bad_data=find(hdg==-9999); hdg(bad_data)=NaN;
       bad_data=find(pit==-9999); pit(bad_data)=NaN;
       bad_data=find(rol==-9999); rol(bad_data)=NaN;
       bad_data=find(uss==-9999); uss(bad_data)=NaN;
       bad_data=find(vss==-9999); vss(bad_data)=NaN;
       bad_data=find(wss==-9999); wss(bad_data)=NaN;
       bad_data=find(ipow==-9999); ipow(bad_data)=NaN;
       bad_data=find(cs==-9999); cs(bad_data)=NaN;
       bad_data=find(cd==-9999); cd(bad_data)=NaN;

       eval_string(iiiNOR(j))={['NOR_' num2str(serialno)]};
       
       eval([char(eval_string(iiiNOR(j))) '.jd=jd;']);
       eval([char(eval_string(iiiNOR(j))) '.t=t;']);
       eval([char(eval_string(iiiNOR(j))) '.p=p;']);
       eval([char(eval_string(iiiNOR(j))) '.u=u;']);
       eval([char(eval_string(iiiNOR(j))) '.v=v;']);
       eval([char(eval_string(iiiNOR(j))) '.w=w;']);
       eval([char(eval_string(iiiNOR(j))) '.hdg=hdg;']);
       eval([char(eval_string(iiiNOR(j))) '.pit=pit;']);
       eval([char(eval_string(iiiNOR(j))) '.rol=rol;']);
       eval([char(eval_string(iiiNOR(j))) '.uss=uss;']);
       eval([char(eval_string(iiiNOR(j))) '.vss=vss;']);
       eval([char(eval_string(iiiNOR(j))) '.wss=wss;']);
       eval([char(eval_string(iiiNOR(j))) '.ipow=ipow;']);
       eval([char(eval_string(iiiNOR(j))) '.cs=cs;']);
       eval([char(eval_string(iiiNOR(j))) '.cd=cd;']);

       sampling_rate = 1/median(diff(jd));

       if filt_or_unfilt==0 % i.e. if want to filter the data prior to plotting
           % Apply a butterworth filter to the data using auto_filt and use for
           % plots
           ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.u));']);
           eval([char(eval_string(iiiNOR(j))) '.u(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
                 '.u(ii), sampling_rate, 1/2,''low'',4);']);
           ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.v));']); 
           eval([char(eval_string(iiiNOR(j))) '.v(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
                 '.v(ii), sampling_rate, 1/2,''low'',4);']);
           ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.t));']); 
           eval([char(eval_string(iiiNOR(j))) '.t(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
                 '.t(ii), sampling_rate, 1/2,''low'',4);']);
       end
       
       current_speed = eval(['sqrt((' char(eval_string(iiiNOR(j))) '.u).^2 + (' char(eval_string(iiiNOR(j))) '.v).^2);']);
       eval([char(eval_string(iiiNOR(j))) '.spd=current_speed;']);
       current_direction=eval(['atan(' char(eval_string(iiiNOR(j))) '.u ./' char(eval_string(iiiNOR(j))) '.v)*180/pi;']);
       d = eval(['find(' char(eval_string(iiiNOR(j))) '.v<0);']);
       current_direction(d)=current_direction(d)+180;
       d = eval(['find(' char(eval_string(iiiNOR(j))) '.v>=0 & ' char(eval_string(iiiNOR(j))) '.u<0);']);
       current_direction(d)=current_direction(d)+360;
       eval([char(eval_string(iiiNOR(j))) '.dir=current_direction;']);
       
       % determine current_speed for setting y-axis after plotting
       max_v=max(abs(v));
       if max_v > MAX_v
           MAX_v = max_v;
       end
       max_v=max(t);
       if max_t > MAX_t
           MAX_t = max_t;
       end
       j=j+1;
    end
end
%--------------------------------------
% Now read in Microcat data if required
%--------------------------------------
if iiMC>0
    j=1;
    
    % loop to read one file at a time

    for i=1:length(vecMC);
       serialno = vecMC(i);
       disp('*************************************************************')
       disp(['Reading MICROCAT - ',num2str(serialno)])
       disp('*************************************************************')

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

           bad_data=find(c==-9999); c(bad_data)=NaN;

           eval_string(iiiMC(j))={['MC_' num2str(serialno)]};

           eval([char(eval_string(iiiMC(j))) '.jd=jd;']);
           eval([char(eval_string(iiiMC(j))) '.c=c;']);
           eval([char(eval_string(iiiMC(j))) '.t=t;']);

           % Apply a butterworth filter to the data using auto_filt and use for
           % plots
           sampling_rate = 1/median(diff(jd));

           if filt_or_unfilt==0 % i.e. if want to filter the data prior to plotting
               % Apply a butterworth filter to the data using auto_filt and use for
               % plots
               ii = eval(['find(~isnan(' char(eval_string(iiiMC(j))) '.t));']); 
               eval([char(eval_string(iiiMC(j))) '.t(ii)=auto_filt(' char(eval_string(iiiMC(j)))...
                     '.t(ii), sampling_rate, 1/2,''low'',4);']);
           end
           
       end
       j=j+1;
    end
end

%%
% Section to find if Microcats are within certain depth of the current meter
matching_depth=10; % set the depth at which instruments are close enough to each other to be considered at the same depth
matching=find(diff(depths(:,2))<=matching_depth); %produces index of matching instrument depths, but only pointing to first instrument of the pair
a=size(matching);
i=1; j=1;
matching_insts=zeros(a(1)*2,3); % pre-allocate array with zeros
while i<=a(1)*2 % while loop creates array of matching instrument types and depths and expands from "matching" index to include both instruments in each pair
    matching_insts(i,2:3)=depths(matching(j),:);
    matching_insts(i,1)=sn(matching(j));
    matching_insts(i+1,2:3)=depths(matching(j)+1,:);
    matching_insts(i+1,1)=sn(matching(j)+1);
    i=i+2;
    j=j+1;
end 


%%
% ------------------------------
% Plotting section
% ------------------------------
% number of current meters given by iii (which doesn't include Microcats)

%set up number of panels for subplot
if length(iii)==1
    m=1;
elseif length(iii)==2 | length(iii)==4
    m=2;
elseif length(iii)==3 | length(iii)==5 | length(iii)==6 | length(iii)==9
    m=3;
elseif length(iii)==7 | length(iii)==8 | length(iii)==10 | length(iii)==11 | length(iii)==12
    m=4;
end

if length(iii)==1 | length(iii)==2 | length(iii)==3
    n=1;
elseif length(iii)>=4 & length(iii)<=8
    n=2;
elseif length(iii)>=9
    n=3;
end

s = regexprep(moor,'_','\\_'); % change mooring text to make compatible with plotting title
eval_string2=eval_string(iii);
for i=1:length(iii)
    figure(correlation_plot);
    subplot(m,n,i)
    hold on
    eval(['plot(' char(eval_string2(i)) '.v,' char(eval_string2(i)) '.t,''.'');'])
    title(['Current meter at ' num2str(depths(iii(i),2)) 'm.']);
    xlabel('v velocity component (cm/s)')
    ylabel('temperature (deg C)')
    
    % now overlay Microcats if close by.
    disp('speak to Harry first to see if worth continuing with this.')
    disp('Or if better to use gridded data.')
    
    
end



if filt_or_unfilt==0 % i.e. filtered data
    suptitle_string=['2-day lowpass filtered temperature vs northward velocity for mooring ' s '.'];
else
    suptitle_string=['Unfiltered temperature vs northward velocity for mooring ' s '.'];
end
suptitle(suptitle_string)
pause
end
