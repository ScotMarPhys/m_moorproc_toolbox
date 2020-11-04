
function currents_overlay(moor,varargin)
%
% Function for plotting directions of a mooring overlayed on the same axes
% and speeds of a mooring overlayed on a second pair of axes.
%
% function currents_overlay('moor','proclvl','layout','plot_interval','procpath')
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
%             standard data paths. 
%             e.g. '/Volumes/noc/mpoc/rpdmoc/rapid/data/moor/proc/'
%   proclvl: can specify level of processing of the data to plot. 
%           e.g. 'proclvl','2': will plot the .use file ; 'proclvl','3' will plot the .microcat and .edt files
%   num_to_plot: number of samples per day to plot. Default is 2 but for
%                short time intervals in plot_interval this may not be 
%                appropriate - NB: cannot be higher than the sampling rate,
%                i.e. not greater number per day than the data actually has
%   unfiltered:  plots data in unfiltered format - useful for more detailed
%                inspection
%
% functions called:-
%   rodbload, julian, auto_filt
%   from .../exec/moor/tools and .../exec/moor/rodb paths
% 
% Routine written by Darren Rayner July 2006.
%
% 25/3/07 - added Nortek capability and made PC compatible
% 15/4/11 - onboard KN200-4: added Seaguard capability
% 05/10/16 - Loic Houpert: added option to process lvl 3 data (.microcat and .edt files for nortek) and save plot
%
if nargin <1
    help currents_overlay
    return
end

varargin_string=varargin;
for i=1:length(varargin) % need to change numeric values in varargin to string so can search with strmatch below
    if isnumeric(varargin{i})
        varargin_string{i}=num2str(varargin{i});
    end
end

% check for optional arguments
a=strmatch('plot_interval',varargin,'exact');
if a>0
    plot_interval=eval(varargin{a+1});
else
    plot_interval=0;
end

a=strmatch('layout',varargin_string,'exact');
if a>0
    %layout=char(varargin{a+1});
    width=26; height=17;
    layout='portrait';    
else
    layout='portrait';
    width=17; height=26;
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

a=strmatch('proclvl',varargin,'exact');
if a>0
    proclvlstr0=char(varargin(a+1));
    proclvl   = str2num(proclvlstr0);
else
    proclvl=2;
    proclvlstr0 = num2str(proclvl);
end


a=strmatch('num_to_plot',varargin_string,'exact');
if a>0
    num_to_plot=varargin{a+1};
else
    num_to_plot=2; 
end

a=strmatch('unfiltered',varargin_string,'exact');
if a>0
    filt_or_unfilt=1;
    proclvlstr = [proclvlstr0 '_unfilt'];    
else
    filt_or_unfilt=0; 
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
if plot_interval==0
    plot_interval = zeros(2,4);
    plot_interval(1,1) = s_d(1); plot_interval(1,2) = s_d(2)-1; plot_interval(1,3) = 1; plot_interval(1,4) = 0;
    plot_interval(2,1) = e_d(1); plot_interval(2,2) = e_d(2)+1; plot_interval(2,3) = 1; plot_interval(2,4) = 0;
    if plot_interval(1,2)==0
        plot_interval(1,2)=12; plot_interval(1,1)=plot_interval(1,1)-1;
    end
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
        plot_x_labels=1; % toggle value to be used later
    elseif i<3
        if (plot_interval(i,1)==plot_interval(i-1,1)) & (plot_interval(i,2)==plot_interval(i-1,2))
            check = 1;
            disp('short plot_interval so xticks may be limited')
            plot_x_labels=0; % toggle value to be used later
        end
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

if plot_x_labels>0

    % cannot have multi-line xticklabels so have to use manual label command
    % this is not really a problem as only want to display years on bottom plot
    year_indexes =[];
    for i=1:length(xticklabels)
        if find(strfind(xticklabels(i,1:3),'Jan'))
            year_indexes=[year_indexes; i];
        end
    end
    % use year_indexes later for plotting on bottom graph
end

jd1 = julian(plot_interval(1,:));
jd2 = julian(plot_interval(2,:)); 



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
iiNOR = find((id == 368|id==370));
vecNOR = sn(iiNOR);
% and find index number of Seaguard RCMs
iiSG = find(id == 301);
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
direction_plot=figure('Position',pos1);
speed_plot=figure('Position',pos1);
www_plot=figure('Position',pos1);   

set(direction_plot,'PaperUnits','centimeters');
set(direction_plot, 'PaperType', 'A4');
set(direction_plot, 'PaperOrientation',layout);
set(speed_plot,'PaperUnits','centimeters');
set(speed_plot, 'PaperType', 'A4');
set(speed_plot, 'PaperOrientation',layout);
set(www_plot,'PaperUnits','centimeters');
set(www_plot, 'PaperType', 'A4');
set(www_plot, 'PaperOrientation',layout);

papersize = get(direction_plot,'PaperSize');
left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
figuresize = [left, bottom, width, height];
set(direction_plot, 'PaperPosition', figuresize);
set(speed_plot, 'PaperPosition', figuresize);
set(www_plot, 'PaperPosition', figuresize);

%set print area of figure
pos1  = [1/8*scnsize(3),8*bdwidth,1/2*scnsize(3),(scnsize(4) - 30*bdwidth)];
pressure_plot=figure('Position',pos1);
set(pressure_plot,'PaperUnits','centimeters');
set(pressure_plot, 'PaperType', 'A4');
set(pressure_plot, 'PaperOrientation',layout);
papersize = get(pressure_plot,'PaperSize');
left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
figuresize = [left, bottom, width, height];
set(pressure_plot, 'PaperPosition', figuresize);

plot_string={};

MAX_spd=0;
MAX_www=0;
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
       end
       
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
       end
       
       current_speed = eval(['sqrt((' char(eval_string(iiiRCM11(j))) '.u).^2 + (' char(eval_string(iiiRCM11(j))) '.v).^2);']);
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
       max_spd=max(current_speed);
       if max_spd > MAX_spd
           MAX_spd = max_spd;
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
           ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.w));']); 
           eval([char(eval_string(iiiNOR(j))) '.w(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
                 '.w(ii), sampling_rate, 1/2,''low'',4);']);            
       end
       
       current_speed = eval(['sqrt((' char(eval_string(iiiNOR(j))) '.u).^2 + (' char(eval_string(iiiNOR(j))) '.v).^2);']);
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

       if filt_or_unfilt==0 % i.e. want to filter the data prior to plotting
           % Apply a butterworth filter to the data using auto_filt and use for
           % plots
           
           ii = eval(['find(~isnan(' char(eval_string(iiiSG(j))) '.u));']);
           eval([char(eval_string(iiiSG(j))) '.u(ii)=auto_filt(' char(eval_string(iiiSG(j)))...
                 '.u(ii), sampling_rate, 1/2,''low'',4);']);
           ii = eval(['find(~isnan(' char(eval_string(iiiSG(j))) '.v));']); 
           eval([char(eval_string(iiiSG(j))) '.v(ii)=auto_filt(' char(eval_string(iiiSG(j)))...
                 '.v(ii), sampling_rate, 1/2,''low'',4);']);
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
       
       j=j+1;
    end
end

% ------------------------------
% Plotting section
% ------------------------------
colours = 'brgcmybrgcmybrgcmy';
for i=1:length(iii)
    if i<7 
        line_style = '-';
    elseif i>12
        line_style = '.';
    else
        line_style = 'o';
    end
    figure(direction_plot); hold on
    a=1;%sampling_rate/num_to_plot; %value to decimate data by to ease plot visualisation doesn't actually average it, just decimates
    % if sampling rate equals num_to_plot entered in function then no decimation occurs.
    % Assumes same sampling rate for all current meters on the mooring.
    %a = 5; 
    b = eval(['length(' char(eval_string(i)) '.jd);']);
    eval(['plot(' char(eval_string(i)) '.jd(1:a:b)-jd1,' char(eval_string(i)) '.dir(1:a:b),''' colours(i) line_style ''',''MarkerSize'',3);']);
    figure(speed_plot); hold on
    eval(['plot(' char(eval_string(i)) '.jd(1:a:b)-jd1,' char(eval_string(i)) '.spd(1:a:b),''' colours(i) line_style ''',''MarkerSize'',3);']);
    figure(www_plot); hold on
    eval(['plot(' char(eval_string(i)) '.jd(1:a:b)-jd1,' char(eval_string(i)) '.www(1:a:b),''' colours(i) line_style ''',''MarkerSize'',3);']);    
end


figure(direction_plot);
ylabel('current direction (deg M)');
xlim([0 jd2-jd1]);
ylim([0 360]);
set(gca,'YMinorTick','on');
set(gca,'xTickLabel',xticklabels);
set(gca,'XTick',jdxticks-jd1);
set(gca,'ytick',[0 90 180 270 360]);
s = regexprep(moor,'_','\\_');
if filt_or_unfilt==0
    title(['Current directions of low pass filtered currents at mooring ' s '.'])
else
    title(['Current directions of unfiltered currents at mooring ' s '.'])
end
for i=1:length(eval_string)
    legend_text(i)={[num2str(depths(i,2)) 'm']};
end
legend(legend_text);
text(1,-20,num2str(xticks(1,1)),'FontSize',10);

if plot_x_labels>0
    for i=1:length(year_indexes)
        text((jd2-jd1)*(year_indexes(i)-1)/(length(xticklabels)-1),-20,num2str(xticks(year_indexes(i),1)),'FontSize',10);
    end
end
print('-dpng',[moor '_horcurrents_overlay_dir_proclvl_' proclvlstr])
savefig([moor '_horcurrents_overlay_dir_proclvl_' proclvlstr])

figure(speed_plot);
ylabel('current speed (cm/s)');
xlim([0 jd2-jd1]);
MAX_spd=ceil(MAX_spd/10)*10;
ylim([0 MAX_spd]);
set(gca,'YMinorTick','on');
set(gca,'xTickLabel',xticklabels);
set(gca,'XTick',jdxticks-jd1);
s = regexprep(moor,'_','\\_');
if filt_or_unfilt==0
    title(['Hor. current speed of low pass filtered currents at mooring ' s '.'])
else
    title(['Hor. current speed of ufiltered currents at mooring ' s '.'])
end
legend(legend_text);
text(1,-5,num2str(xticks(1,1)),'FontSize',10);

if plot_x_labels>0
    for i=1:length(year_indexes)
        text((jd2-jd1)*(year_indexes(i)-1)/(length(xticklabels)-1),-5,num2str(xticks(year_indexes(i),1)),'FontSize',10);
    end
end


print('-dpng',[moor '_horcurrents_overlay_speed_proclvl_' proclvlstr])
savefig([moor '_horcurrents_overlay_speed_proclvl_' proclvlstr])


figure(www_plot);
ylabel('current speed (cm/s)');
xlim([0 jd2-jd1]);
MAX_www=ceil(MAX_www/10)*10;
ylim([0 MAX_www]);
set(gca,'YMinorTick','on');
set(gca,'xTickLabel',xticklabels);
set(gca,'XTick',jdxticks-jd1);
s = regexprep(moor,'_','\\_');
if filt_or_unfilt==0
    title(['Vert. current speed of low pass filtered currents at mooring ' s '.'])
else
    title(['Vert. current speed of ufiltered currents at mooring ' s '.'])
end
legend(legend_text);
text(1,-5,num2str(xticks(1,1)),'FontSize',10);

if plot_x_labels>0
    for i=1:length(year_indexes)
        text((jd2-jd1)*(year_indexes(i)-1)/(length(xticklabels)-1),-5,num2str(xticks(year_indexes(i),1)),'FontSize',10);
    end
end


print('-dpng',[moor '_horcurrents_overlay_speed_proclvl_' proclvlstr])
savefig([moor '_horcurrents_overlay_speed_proclvl_' proclvlstr])
