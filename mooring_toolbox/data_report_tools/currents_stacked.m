function currents_stacked(moor,varargin)
%
% Function for plotting u,v, speed and direction plots from RCM11s, S4s and
% Argonauts. Composite plot by mooring.
%
% function currents_stacked('moor','proclvl','layout','plot_interval','procpath','unfiltered')
%
% required inputs:-
%   moor: complete mooring name as string. e.g. 'wb1_1_200420'
%
% optional inputs:-
%   layout: orientation of figure portrait/lanscape (default = portrait)
%           input of 1 = landscape, 0 = portrait
%   plot_interval: matrix of start and end dates for plot
%                  e.g. [2004 02 01 00; 2005 06 01 00]
%                  dates are:- yyyy mm dd hh
%   procpath: can specify exact path to proc directory if not using 
%           standard data paths. 
%           e.g. '/Volumes/noc/mpoc/hydro/rpdmoc/rapid/data/moor/proc/'
%   proclvl: can specify level of processing of the data to plot. 
%           e.g. 'proclvl','2': will plot the .use file ; 'proclvl','3' will plot the .microcat and .edt files
%
% functions called:-
%   rodbload, julian, auto_filt
%   from .../exec/moor/tools and .../exec/moor/rodb paths
%   suptitle.m downloaded from Mathworks file exchange website 
%
% Routine written by Darren Rayner January 2008 - adapted from stick_plot.
%
% 15/4/2011 - aboard cruise kn200-4: added Seaguard functionality 
% 05/10/16 - Loic Houpert: added option to process lvl 3 data (.microcat and .edt files for nortek) and save plot
%

if nargin <1
    help currents_stacked
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

a=strmatch('proclvl',varargin,'exact');
if a>0
    proclvlstr0=char(varargin(a+1));
    proclvl   = str2num(proclvlstr0);
else
    proclvl=2;
    proclvlstr0 = num2str(proclvl);
end

a=strmatch('plot_interval',varargin,'exact');
if a>0
    plot_interval=eval(varargin{a+1});
else
    plot_interval=0;
end

a=strmatch('unfiltered',varargin,'exact');
if a>0
    unfilt=1;
    proclvlstr = [proclvlstr0 '_unfilt'];
else
    unfilt=0;
    proclvlstr = [proclvlstr0 '_lpfilt'];    
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
jd_start = julian([s_d' hms2h([s_t;0]')]);
jd_end   = julian([e_d' hms2h([e_t;0]')]);

disp(['z : instrument id : serial number'])
for i = 1:length(id)
    disp([z(i),id(i),sn(i)])
end

% ------------------------------------------------
% Determine plot_interval if not input to function
if nargin == 1
    plot_interval = zeros(2,4);
    plot_interval(1,1) = s_d(1); plot_interval(1,2) = s_d(2); plot_interval(1,3) = 1; plot_interval(1,4) = 0;
    plot_interval(2,1) = e_d(1); plot_interval(2,2) = e_d(2)+1; plot_interval(2,3) = 1; plot_interval(2,4) = 0;
    if plot_interval(2,2)==13
        plot_interval(2,2)=1; plot_interval(2,1)=plot_interval(2,1)+1;
    end
end

% ------------------------------------------------
% Determine plot_interval if not input to function
if plot_interval==0
    plot_interval = zeros(2,4);
    plot_interval(1,1) = s_d(1); plot_interval(1,2) = s_d(2); plot_interval(1,3) = 1; plot_interval(1,4) = 0;
    plot_interval(2,1) = e_d(1); plot_interval(2,2) = e_d(2)+1; plot_interval(2,3) = 1; plot_interval(2,4) = 0;
    if plot_interval(2,2)==13
        plot_interval(2,2)=1; plot_interval(2,1)=plot_interval(2,1)+1;
    end
end


% create xtick spacings based on start of months     
check=0;
i=2;
xticks(1,:)=plot_interval(1,:);
while check~=1
    xticks(i,:)=xticks(i-1,:);
    if xticks(i,2)<12
        xticks(i,2)=xticks(i-1,2)+1;
    else
        xticks(i,2)=1;
        xticks(i,1)=xticks(i-1,1)+1;
    end
    if xticks(i,:)>=plot_interval(2,:)
        check = 1;
    end
    i=i+1;
end

if i<4
   jdxticks =julian(plot_interval(1,:)):(julian(plot_interval(2,:))-julian(plot_interval(1,:)))/5:julian(plot_interval(2,:));
   gxticks = gregorian(jdxticks);
   xticks = gxticks(:,1:4);
   xticklabels = datestr(gxticks,'dd mmm');
else
	jdxticks=julian(xticks);
	% create xticklabels from xticks
	months=['Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'];
	xticklabels=months(xticks(:,2),1:3);
end

% cannot have multi-line xticklabels so have to use manual label command
% this is not really a problem as only want to display years on bottom plot
year_indexes =[];
for i=1:length(xticklabels)
    if find(strfind(xticklabels(i,1:3),'Jan'))
        year_indexes=[year_indexes; i];
    end
end
% use year_indexes later for plotting on bottom graph

jd1 = julian(plot_interval(1,:));
jd2 = julian(plot_interval(2,:)); 


% find the index number of S4s
iiS4 = find(id == 302);
vecS4 = sn(iiS4);
% and find index number of RCM11s
iiRCM11 = find(id == 310);
vecRCM11 = sn(iiRCM11);
% and find index number of Sontek Argonauts
iiARG = find(id == 366 | id == 366337);
vecARG = sn(iiARG);
% and find index number of Nortek Aquadopps
iiNOR = find((id ==368|id==370));
vecNOR = sn(iiNOR);
% and find index number of Nortek Aquadopps
iiSG = find(id ==301);
vecSG = sn(iiSG);

depths(:,1) = id([iiS4;iiRCM11;iiARG;iiNOR;iiSG]);
depths(:,2) = z([iiS4;iiRCM11;iiARG;iiNOR;iiSG]);
depths=sortrows(depths,2);
iiiS4=find(depths(:,1)==302);  
iiiRCM11=find(depths(:,1)==310);  
iiiARG=find(depths(:,1)==366 | depths(:,1)==366337);  
iiiNOR=find(depths(:,1)==368|depths(:,1)==370);
iiiSG=find(depths(:,1)==301);
iii=[iiiS4;iiiRCM11;iiiARG;iiiNOR;iiiSG];

%set figure size on screen for better viewing
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

%set print area of figure
pos1  = [1/8*scnsize(3),8*bdwidth,1/2*scnsize(3),(scnsize(4) - 30*bdwidth)];
hordirection_plot=figure('Position',pos1);
horspeed_plot=figure('Position',pos1);
u_plot=figure('Position',pos1);
v_plot=figure('Position',pos1);
w_plot=figure('Position',pos1);

set(hordirection_plot,'PaperUnits','centimeters');
set(hordirection_plot, 'PaperType', 'A4');
set(hordirection_plot, 'PaperOrientation',layout);
set(horspeed_plot,'PaperUnits','centimeters');
set(horspeed_plot, 'PaperType', 'A4');
set(horspeed_plot, 'PaperOrientation',layout);
set(u_plot,'PaperUnits','centimeters');
set(u_plot, 'PaperType', 'A4');
set(u_plot, 'PaperOrientation',layout);
set(v_plot,'PaperUnits','centimeters');
set(v_plot, 'PaperType', 'A4');
set(v_plot, 'PaperOrientation',layout);
set(w_plot,'PaperUnits','centimeters');
set(w_plot, 'PaperType', 'A4');
set(w_plot, 'PaperOrientation',layout);

papersize = get(hordirection_plot,'PaperSize');
width=17; height=26; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
figuresize = [left, bottom, width, height];
set(hordirection_plot, 'PaperPosition', figuresize);
set(horspeed_plot, 'PaperPosition', figuresize);
set(u_plot, 'PaperPosition', figuresize);
set(v_plot, 'PaperPosition', figuresize);
set(w_plot, 'PaperPosition', figuresize);

MAX_u=0;
MAX_v=0;
MAX_www=0;
MAX_spd=0;
plot_string={};

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

       % Apply a butterworth filter to the data using auto_filt and use for
       % plots
       sampling_rate = 1/median(diff(jd));
       
       uf=[];
       vf=[];
       
       ii = eval(['find(~isnan(' char(eval_string(iiiS4(j))) '.u));']); 
       eval(['uf=auto_filt(' char(eval_string(iiiS4(j)))...
             '.u(ii), sampling_rate, 1/2,''low'',4);']);
       
       ii = eval(['find(~isnan(' char(eval_string(iiiS4(j))) '.v));']); 
       eval(['vf=auto_filt(' char(eval_string(iiiS4(j)))...
             '.v(ii), sampling_rate, 1/2,''low'',4);']);
       
       if unfilt==0
           eval([char(eval_string(iiiS4(j))) '.v(ii)=vf;']);
           eval([char(eval_string(iiiS4(j))) '.u(ii)=uf;']);
       else
           % do not reassign u and v as filtered version
       end
       
       % calculate and filter speed and direction
       current_speed=[];
       current_direction=[];
       
       current_speed = eval(['sqrt((' char(eval_string(iiiS4(j))) '.u).^2 + (' char(eval_string(iiiS4(j))) '.v) .^2);']);
       eval([char(eval_string(iiiS4(j))) '.spd=current_speed;']);
       current_direction=eval(['atan(' char(eval_string(iiiS4(j))) '.u ./' char(eval_string(iiiS4(j))) '.v)*180/pi;']);
       d = eval(['find(' char(eval_string(iiiS4(j))) '.v<0);']);
       current_direction(d)=current_direction(d)+180;
       d = eval(['find(' char(eval_string(iiiS4(j))) '.v>=0 & ' char(eval_string(iiiS4(j))) '.u<0);']);
       current_direction(d)=current_direction(d)+360;
       eval([char(eval_string(iiiS4(j))) '.dir=current_direction;']);
      
       % determine current_speed for setting y-axis after plotting
       max_spd=max(current_speed);
       if max_spd > MAX_spd
           MAX_spd = max_spd;
       end
       
       % determine maximum u and v components for setting y-axis after plotting
       if unfilt==0
           max_u=max(abs(uf));
           max_v=max(abs(vf));
       else
           max_u=max(abs(u));
           max_v=max(abs(v));
       end
       if max_u > MAX_u
           MAX_u = max_u;
       end
       if max_v > MAX_v
           MAX_v = max_v;
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

       
       % Apply a butterworth filter to the data using auto_filt and use for
       % plots
       sampling_rate = 1/median(diff(jd));

       %ii = eval(['find(~isnan(' char(eval_string(iiiRCM11(j))) '.u));']); 
       %eval([char(eval_string(iiiRCM11(j))) '.u(ii)=auto_filt(' char(eval_string(iiiRCM11(j)))...
       %      '.u(ii), sampling_rate, 1/2,''low'',4);']);
       %ii = eval(['find(~isnan(' char(eval_string(iiiRCM11(j))) '.v));']); 
       %eval([char(eval_string(iiiRCM11(j))) '.v(ii)=auto_filt(' char(eval_string(iiiRCM11(j)))...
       %      '.v(ii), sampling_rate, 1/2,''low'',4);']);
       
       uf=[];
       vf=[];
       
       ii = eval(['find(~isnan(' char(eval_string(iiiRCM11(j))) '.u));']); 
       eval(['uf=auto_filt(' char(eval_string(iiiRCM11(j)))...
             '.u(ii), sampling_rate, 1/2,''low'',4);']);
       
       ii = eval(['find(~isnan(' char(eval_string(iiiRCM11(j))) '.v));']); 
       eval(['vf=auto_filt(' char(eval_string(iiiRCM11(j)))...
             '.v(ii), sampling_rate, 1/2,''low'',4);']);
       if unfilt==0
           eval([char(eval_string(iiiRCM11(j))) '.u(ii)=uf;']);
           eval([char(eval_string(iiiRCM11(j))) '.v(ii)=vf;']);
       else
           % do no reassign u and v as filtered versions
       end
             
             
       % calculate speed and direction
       current_speed=[];
       current_direction=[];
       
       current_speed = eval(['sqrt((' char(eval_string(iiiRCM11(j))) '.u).^2 + (' char(eval_string(iiiRCM11(j))) '.v) .^2);']);
       eval([char(eval_string(iiiRCM11(j))) '.spd=current_speed;']);
       current_direction=eval(['atan(' char(eval_string(iiiRCM11(j))) '.u ./' char(eval_string(iiiRCM11(j))) '.v)*180/pi;']);
       d = eval(['find(' char(eval_string(iiiRCM11(j))) '.v<0);']);
       current_direction(d)=current_direction(d)+180;
       d = eval(['find(' char(eval_string(iiiRCM11(j))) '.v>=0 & ' char(eval_string(iiiRCM11(j))) '.u<0);']);
       current_direction(d)=current_direction(d)+360;
       eval([char(eval_string(iiiRCM11(j))) '.dir=current_direction;']);
      
       % determine current_speed for setting y-axis after plotting
       max_spd=max(current_speed);
       if max_spd > MAX_spd
           MAX_spd = max_spd;
       end
       
       % determine maximum u and v components for setting y-axis after plotting
       if unfilt==0
           max_u=max(abs(uf));
           max_v=max(abs(vf));
       else
           max_u=max(abs(u));
           max_v=max(abs(v));
       end
       if max_u > MAX_u
           MAX_u = max_u;
       end
       if max_v > MAX_v
           MAX_v = max_v;
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
       
       % Apply a butterworth filter to the data using auto_filt and use for
       % plots
       sampling_rate = 1/median(diff(jd));

       %ii = eval(['find(~isnan(' char(eval_string(iiiARG(j))) '.u));']); 
       %eval([char(eval_string(iiiARG(j))) '.u(ii)=auto_filt(' char(eval_string(iiiARG(j)))...
       %      '.u(ii), sampling_rate, 1/2,''low'',4);']);
       %ii = eval(['find(~isnan(' char(eval_string(iiiARG(j))) '.v));']); 
       %eval([char(eval_string(iiiARG(j))) '.v(ii)=auto_filt(' char(eval_string(iiiARG(j)))...
       %      '.v(ii), sampling_rate, 1/2,''low'',4);']);
       
       uf=[];
       vf=[];
       
       ii = eval(['find(~isnan(' char(eval_string(iiiARG(j))) '.u));']); 
       eval(['uf=auto_filt(' char(eval_string(iiiARG(j)))...
             '.u(ii), sampling_rate, 1/2,''low'',4);']);
       
       ii = eval(['find(~isnan(' char(eval_string(iiiARG(j))) '.v));']); 
       eval(['vf=auto_filt(' char(eval_string(iiiARG(j)))...
             '.v(ii), sampling_rate, 1/2,''low'',4);']);
       if unfilt==0
           eval([char(eval_string(iiiARG(j))) '.u(ii)=uf;']);
           eval([char(eval_string(iiiARG(j))) '.v(ii)=vf;']);
       else
           % do not reassign u and v as filtered versions
       end
       
       % calculate speed and direction
       current_speed=[];
       current_direction=[];
       
       current_speed = eval(['sqrt((' char(eval_string(iiiARG(j))) '.u).^2 + (' char(eval_string(iiiARG(j))) '.v) .^2);']);
       eval([char(eval_string(iiiARG(j))) '.spd=current_speed;']);
       current_direction=eval(['atan(' char(eval_string(iiiARG(j))) '.u ./' char(eval_string(iiiARG(j))) '.v)*180/pi;']);
       d = eval(['find(' char(eval_string(iiiARG(j))) '.v<0);']);
       current_direction(d)=current_direction(d)+180;
       d = eval(['find(' char(eval_string(iiiARG(j))) '.v>=0 & ' char(eval_string(iiiARG(j))) '.u<0);']);
       current_direction(d)=current_direction(d)+360;
       eval([char(eval_string(iiiARG(j))) '.dir=current_direction;']);
      
       % determine current_speed for setting y-axis after plotting
       max_spd=max(current_speed);
       if max_spd > MAX_spd
           MAX_spd = max_spd;
       end
       
       % determine maximum u and v components for setting y-axis after plotting
       if unfilt==0
           max_u=max(abs(uf));
           max_v=max(abs(vf));
       else
           max_u=max(abs(u));
           max_v=max(abs(v));
       end
       
       if max_u > MAX_u
           MAX_u = max_u;
       end
       if max_v > MAX_v
           MAX_v = max_v;
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

	if proclvl==2
	       if isunix
        	   infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.use'];
       		elseif ispc
        	   infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.use'];
       		end
	elseif proclvl==3
	       if isunix
        	   infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
       		elseif ispc
        	   infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
       		end	
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

       
       %ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.u));']); 
       %eval([char(eval_string(iiiNOR(j))) '.u(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
       %      '.u(ii), sampling_rate, 1/2,''low'',4);']);
       %ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.v));']); 
       %eval([char(eval_string(iiiNOR(j))) '.v(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
       %      '.v(ii), sampling_rate, 1/2,''low'',4);']);

       % Apply a butterworth filter to the data using auto_filt and use for
       % plots
       
       sampling_rate = 1/median(diff(jd));

       
       if unfilt==0
           
       uf=[];
       vf=[];          
           if ~isempty(eval([ char(eval_string(iiiNOR(j))) '.u']))
           ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.u));']); 
           eval(['uf=auto_filt(' char(eval_string(iiiNOR(j)))...
                 '.u(ii), sampling_rate, 1/2,''low'',4);']);

           ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.v));']); 
           eval(['vf=auto_filt(' char(eval_string(iiiNOR(j)))...
                 '.v(ii), sampling_rate, 1/2,''low'',4);']);           
               eval([char(eval_string(iiiNOR(j))) '.u(ii)=uf;']);
               eval([char(eval_string(iiiNOR(j))) '.v(ii)=vf;']);
           else
           end
       else
           % do not reassign u and v with filtered versions
       end
       
       % calculate speed and direction
       current_speed=[];
       current_direction=[];
       
       current_speed = eval(['sqrt((' char(eval_string(iiiNOR(j))) '.u).^2 + (' char(eval_string(iiiNOR(j))) '.v) .^2);']);
       eval([char(eval_string(iiiNOR(j))) '.spd=current_speed;']);
       vert_speed = eval([ char(eval_string(iiiNOR(j))) '.w;']);
       eval([char(eval_string(iiiNOR(j))) '.www=vert_speed;']);
       current_direction=eval(['atan(' char(eval_string(iiiNOR(j))) '.u ./' char(eval_string(iiiNOR(j))) '.v)*180/pi;']);
       d = eval(['find(' char(eval_string(iiiNOR(j))) '.v<0);']);
       current_direction(d)=current_direction(d)+180;
       d = eval(['find(' char(eval_string(iiiNOR(j))) '.v>=0 & ' char(eval_string(iiiNOR(j))) '.u<0);']);
       current_direction(d)=current_direction(d)+360;
       eval([char(eval_string(iiiNOR(j))) '.dir=current_direction;']);
      
       % determine current_speed for setting y-axis after plotting
       max_spd=max(current_speed);
       if max_spd > MAX_spd
           MAX_spd = max_spd;
       end
       
       % determine maximum u and v components for setting y-axis after plotting
       if unfilt==0
           max_u=max(abs(uf));
           max_v=max(abs(vf));
       else
           max_u=max(abs(u));
           max_v=max(abs(v));
       end
       if max_u > MAX_u
           MAX_u = max_u;
       end
       if max_v > MAX_v
           MAX_v = max_v;
       end
         max_www=max(vert_speed);
       if max_www > MAX_www
           MAX_www = max_www;
       end        
       j=j+1;
    end
end
%--------------------------------------
% Now read in Seaguard data if required
%--------------------------------------
if iiSG>0
    j=1;
    
    % loop to read one file at a time

    for i=1:length(vecSG);
       serialno = vecSG(i);
       disp('*************************************************************')
       disp(['Reading SG - ',num2str(serialno)])
       disp('*************************************************************')

       if isunix
           infile = [procpath,moor,'/seaguard/',moor,'_',sprintf('%3.3d',vecSG(i)),'.use'];
       elseif ispc
           infile = [procpath,moor,'\seaguard\',moor,'_',sprintf('%3.3d',vecSG(i)),'.use'];
       end
       

       % read data into vectors and then into structure array

       [yy,mm,dd,hh,u,v,cs,cd,cssd,mss,hdg,pit,rol,t,c,tc,p,tp,ipow] = ...
            rodbload(infile,'YY:MM:DD:HH:U:V:CS:CD:CSSD:MSS:HDG:PIT:ROL:T:C:TC:P:TP:IPOW');
       jd=julian(yy,mm,dd,hh);

       bad_data=find(t==-9999); t(bad_data)=NaN;
       bad_data=find(p==-9999); p(bad_data)=NaN;
       bad_data=find(u==-9999); u(bad_data)=NaN;
       bad_data=find(v==-9999); v(bad_data)=NaN;
       bad_data=find(c==-9999); c(bad_data)=NaN;
       bad_data=find(hdg==-9999); hdg(bad_data)=NaN;
       bad_data=find(pit==-9999); pit(bad_data)=NaN;
       bad_data=find(rol==-9999); rol(bad_data)=NaN;

       eval_string(iiiSG(j))={['SG_' num2str(serialno)]};
       
       eval([char(eval_string(iiiSG(j))) '.jd=jd;']);
       eval([char(eval_string(iiiSG(j))) '.t=t;']);
       eval([char(eval_string(iiiSG(j))) '.p=p;']);
       eval([char(eval_string(iiiSG(j))) '.u=u;']);
       eval([char(eval_string(iiiSG(j))) '.v=v;']);
       eval([char(eval_string(iiiSG(j))) '.c=c;']);
       eval([char(eval_string(iiiSG(j))) '.hdg=hdg;']);
       eval([char(eval_string(iiiSG(j))) '.pit=pit;']);
       eval([char(eval_string(iiiSG(j))) '.rol=rol;']);

       sampling_rate = 1/median(diff(jd));

       uf=[];
       vf=[];
       
       ii = eval(['find(~isnan(' char(eval_string(iiiSG(j))) '.u));']); 
       eval(['uf=auto_filt(' char(eval_string(iiiSG(j)))...
             '.u(ii), sampling_rate, 1/2,''low'',4);']);
       
       ii = eval(['find(~isnan(' char(eval_string(iiiSG(j))) '.v));']); 
       eval(['vf=auto_filt(' char(eval_string(iiiSG(j)))...
             '.v(ii), sampling_rate, 1/2,''low'',4);']);
       if unfilt==0
           eval([char(eval_string(iiiSG(j))) '.u(ii)=uf;']);
           eval([char(eval_string(iiiSG(j))) '.v(ii)=vf;']);
       else
           % do no reassign u and v as filtered versions
       end
       
       current_speed = eval(['sqrt((' char(eval_string(iiiSG(j))) '.u).^2 + (' char(eval_string(iiiSG(j))) '.v).^2);']);
       eval([char(eval_string(iiiSG(j))) '.spd=current_speed;']);
       current_direction=eval(['atan(' char(eval_string(iiiSG(j))) '.u ./' char(eval_string(iiiSG(j))) '.v)*180/pi;']);
       d = eval(['find(' char(eval_string(iiiSG(j))) '.v<0);']);
       current_direction(d)=current_direction(d)+180;
       d = eval(['find(' char(eval_string(iiiSG(j))) '.v>=0 & ' char(eval_string(iiiSG(j))) '.u<0);']);
       current_direction(d)=current_direction(d)+360;
       eval([char(eval_string(iiiSG(j))) '.dir=current_direction;']);
       
       % determine current_speed for setting y-axis after plotting
       max_spd=max(current_speed);
       if max_spd > MAX_spd
           MAX_spd = max_spd;
       end
       
       % determine maximum u and v components for setting y-axis after plotting
       if unfilt==0
           max_u=max(abs(uf));
           max_v=max(abs(vf));
       else
           max_u=max(abs(u));
           max_v=max(abs(v));
       end
       if max_u > MAX_u
           MAX_u = max_u;
       end
       if max_v > MAX_v
           MAX_v = max_v;
       end
       
       j=j+1;
    end
end

% ------------------------------
% Plotting section
% ------------------------------
% calculate limits of u and v axes using MAX_u and MAX_v
MAX_u=ceil(MAX_u/10)*10;
MAX_v=ceil(MAX_v/10)*10;
% calculate limits of spd axes using MAX_spd
MAX_spd=ceil(MAX_spd/10)*10;

plot_handle={'horspeed_plot','hordirection_plot','u_plot','v_plot','w_plot'};
ylabels={'speed (cm/s)','direction (deg)','velocity (cm/s)','velocity (cm/s)','velocity (cm/s)'};
yvalue={'spd', 'dir', 'u', 'v', 'w'};
suptitle_label={'hor. current speed','hor. current direction','u-velocity component','v-velocity component','w-velocity component'};
ylimits={[0 MAX_spd],[0 360],[-MAX_u MAX_u],[-MAX_v MAX_v],[-MAX_www MAX_www]};

for k=1:length(plot_handle)
    eval(['figure(' plot_handle{k} ')']);
    % setup subplot number of panels
    subplot(length(depths(:,1)),1,1);

    % create axes
    for i=1:length(iii)
        % create axes
        axes_string=['axes' num2str(iii(i))];
        eval([axes_string '=subplot(length(depths(:,1)),1,iii(i));']);
        
        % plot data
        eval(['plot((' eval_string{i} '.jd-jd1),' eval_string{i} '.' yvalue{k} ');']);
        
        eval(['ylabel(''' ylabels{k} ''');']);
        xlim([0 jd2-jd1]);
        set(gca,'YMinorTick','on');
        set(gca,'xTickLabel',xticklabels);
        set(gca,'XTick',jdxticks-jd1);
        
        Y_limits=(ylimits{k});
        ylim(Y_limits);
        X_limits=xlim;
        set(gca,'YMinorTick','on');
        set(gca,'xTickLabel',xticklabels);
        set(gca,'XTick',(jdxticks-jd1));
        % Not using timeaxis function.
        
        if k==2
            set(gca,'ytick',[0 90 180 270 360])
            set(gca,'yminortick','off')
        end
        
        % draw zero axis
        hold on
        plot(X_limits,[0 0],'k')
        
        % Label plots with serial numbers and depths
        % need to insert backslash into text string to prevent subscript in
        % labels
        s = regexprep(eval_string(iii(i)), '_', '\nsn:');
        text((X_limits(2)-X_limits(1))*1.005+X_limits(1), (Y_limits(2)-Y_limits(1))*0.75+Y_limits(1),s,'FontSize',8);

        text((X_limits(2)-X_limits(1))*1.005+X_limits(1), (Y_limits(2)-Y_limits(1))*0.5+Y_limits(1),[num2str(depths(iii(i),2)) 'm'],'FontSize',8);

    end

    
    % Display year labels on bottom graph
    eval(['axes(axes' num2str(i) ');']);
    a=(Y_limits(1)-Y_limits(2))*1.5+Y_limits(2);
    text(X_limits(1),a,num2str(xticks(1,1)),'FontSize',10);
    for i=1:length(year_indexes)
        text((jd2-jd1)*(year_indexes(i)-1)/(length(xticklabels)-1),a,num2str(xticks(year_indexes(i),1)),'FontSize',10);
    end

    % Add title to top graph
    s = regexprep(moor,'_','\\_');
    % suptitle is function to place a title over subplots - download from the
    % Mathworks file exchange.
    if unfilt==0
        suptitle(['Low-pass filtered ' suptitle_label{k} ' from mooring: ' s]);
    else
        suptitle(['Unfiltered ' suptitle_label{k} ' from mooring: ' s]);
    end
    
    print('-dpng',[moor '_currents_stacked_' plot_handle{k} '_proclvl_' proclvlstr])
    savefig([moor '_currents_stacked_' plot_handle{k} '_proclvl_' proclvlstr])

end
