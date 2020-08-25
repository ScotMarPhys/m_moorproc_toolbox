function progressive_vector(moor,varargin)
%
% Function for plotting progressive vector plots of velocity from RCM11s, S4s and
% Argonauts. Composite plot by mooring.
%
% function progressive_vector('moor','proclvl','layout','procpath')
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
%   rodbload, julian, auto_filt, suptitle (super title for subplot figures)
%   from .../exec/moor/tools and .../exec/moor/rodb paths
% 
% Routine written by Darren Rayner March 2007.
% 15/4/11 - aboard KN200-4 added Seaguard compatability
% 05/10/16 - Loic Houpert: added option to process lvl 3 data (.microcat
% and .edt files for nortek) and save plot; add plot_interval option
%
if nargin <1
    help progressive_vector
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
    proclvlstr=char(varargin(a+1));
    proclvl   = str2num(proclvlstr);
else
    proclvl=2;
    proclvlstr = num2str(proclvl);
end

a=strmatch('plot_interval',varargin,'exact');
if a>0
    plot_interval=eval(varargin{a+1});
else
    plot_interval=0;
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

% Determine plot_interval if not input to function
if plot_interval==0
    plot_interval = zeros(2,4);
    plot_interval(1,1) = s_d(1); plot_interval(1,2) = s_d(2); plot_interval(1,3) = 1; plot_interval(1,4) = 0;
    plot_interval(2,1) = e_d(1); plot_interval(2,2) = e_d(2)+1; plot_interval(2,3) = 1; plot_interval(2,4) = 0;
    if plot_interval(2,2)==13
        plot_interval(2,2)=1; plot_interval(2,1)=plot_interval(2,1)+1;
    end
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
iiARG = find(id == 366 | id == 366337);
vecARG = sn(iiARG);
% and find index number of Nortek Aquadopps
iiNOR = find((id==368|id==370));
vecNOR = sn(iiNOR);
% and find index number of Seaguard RCMs
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
prog_plot=figure('Position',pos1);
set(prog_plot,'PaperUnits','centimeters');
set(prog_plot, 'PaperType', 'A4');
set(prog_plot, 'PaperOrientation',layout);
papersize = get(prog_plot,'PaperSize');
width=17; height=26; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
figuresize = [left, bottom, width, height];
set(prog_plot, 'PaperPosition', figuresize);

% setup subplot number of panels
if length(depths)<=3
    across=1;
    down=length(depths);
elseif length(depths)==4
    across=2; down=2;
elseif 4<length(depths) & length(depths)<=7
    across=2; down=3;
elseif 6<length(depths) & length(depths)<=9
    across=3;down=3;
elseif 9<length(depths) & length(depths)<=12
    across=3;down=4;
end
plot_string={};

% % --------------------
% % Read in S4 data if required.
% % --------------------
% if iiS4>0 
%     
%     % loop to read one file at a time
%     j=1;
%     for i=1:length(vecS4);
%        serialno = vecS4(i);
%        disp('*************************************************************')
%        disp(['Reading S4 - ',num2str(serialno)])
%        disp('*************************************************************')
% 
%        if isunix
%            infile = [procpath,moor,'/s4/',moor,'_',sprintf('%4.4d',vecS4(i)),'.use'];
%        elseif ispc
%            infile = [procpath,moor,'\s4\',moor,'_',sprintf('%4.4d',vecS4(i)),'.use'];
%        end
% 
%        % read data into vectors and then into structure array
% 
%        [yy,mm,dd,hh,u,v,t,c,p,hdg] = rodbload(infile,'yy:mm:dd:hh:u:v:t:c:p:hdg');
%        jd=julian(yy,mm,dd,hh);
% 
%        bad_data=find(u==-9999); u(bad_data)=NaN;
%        bad_data=find(v==-9999); v(bad_data)=NaN;
%        bad_data=find(t==-9999); t(bad_data)=NaN;
%        bad_data=find(c==-9999); c(bad_data)=NaN;
%        bad_data=find(p==-9999); p(bad_data)=NaN;
%        bad_data=find(hdg==-9999); hdg(bad_data)=NaN;
% 
%        eval_string(iiiS4(j))={['S4_' num2str(serialno)]};
%        
%        eval([char(eval_string(iiiS4(j))) '.jd=jd;']);
%        eval([char(eval_string(iiiS4(j))) '.u=u;']);
%        eval([char(eval_string(iiiS4(j))) '.v=v;']);
%        eval([char(eval_string(iiiS4(j))) '.t=t;']);
%        eval([char(eval_string(iiiS4(j))) '.c=c;']);
%        eval([char(eval_string(iiiS4(j))) '.p=p;']);
%        eval([char(eval_string(iiiS4(j))) '.hdg=hdg;']);
% 
%        % Apply a butterworth filter to the data using auto_filt and use for
%        % plots
%        sampling_rate = 1/median(diff(jd));
% 
%        ii = eval(['find(~isnan(' char(eval_string(iiiS4(j))) '.u));']); 
%        eval([char(eval_string(iiiS4(j))) '.u(ii)=auto_filt(' char(eval_string(iiiS4(j)))...
%              '.u(ii), sampling_rate, 1/2,''low'',4);']);
%        ii = eval(['find(~isnan(' char(eval_string(iiiS4(j))) '.v));']); 
%        eval([char(eval_string(iiiS4(j))) '.v(ii)=auto_filt(' char(eval_string(iiiS4(j)))...
%              '.v(ii), sampling_rate, 1/2,''low'',4);']);
% 
%        
%        % Calculate cumulative vectors
%        clear prog_u prog_v;
%        prog_u(1)=0;
%        prog_v(1)=0;
%        for k=2:length(u)
%            prog_u(k)=u(k-1)+prog_u(k-1);
%            prog_v(k)=v(k-1)+prog_v(k-1);
%        end
%        eval([char(eval_string(iiiS4(j))) '.prog_u=prog_u;']);
%        eval([char(eval_string(iiiS4(j))) '.prog_v=prog_v;']);
%        
%        % decimate to 1 record every ten days to use as time stamp for plot
%        b=round(sampling_rate)*10;
%        clear prog_u2; clear prog_v2; clear jd_2;
%        prog_u2=prog_u(1:b:length(prog_u));
%        prog_v2=prog_v(1:b:length(prog_v));
%        jd_2=jd(1:b:length(jd));
%        eval([char(eval_string(iiiS4(j))) '.prog_u2=prog_u2;']);
%        eval([char(eval_string(iiiS4(j))) '.prog_v2=prog_v2;']);
%        eval([char(eval_string(iiiS4(j))) '.jd_2=jd_2;']);
%        
%        j=j+1;
%        
%     end
% end
% 
% %-----------------------------------
% % Now read in RCM11 data if required
% %-----------------------------------
% if iiRCM11>0
%     j=1;
%     
%     % loop to read one file at a time
% 
%     for i=1:length(vecRCM11);
%        serialno = vecRCM11(i);
%        disp('*************************************************************')
%        disp(['Reading RCM11 - ',num2str(serialno)])
%        disp('*************************************************************')
%        
%        if isunix
%            infile = [procpath,moor,'/rcm/',moor,'_',sprintf('%3.3d',vecRCM11(i)),'.use'];
%        elseif ispc
%            infile = [procpath,moor,'\rcm\',moor,'_',sprintf('%3.3d',vecRCM11(i)),'.use'];
%        end
%        
%        % read data into vectors and then into structure array
% 
%        [yy,mm,dd,hh,ref,u,v,t,c,p,tlt,mss] = rodbload(infile,'yy:mm:dd:hh:ref:u:v:t:c:p:tlt:mss');
%        jd=julian(yy,mm,dd,hh);
% 
%        bad_data=find(t==-9999); t(bad_data)=NaN;
%        bad_data=find(p==-9999); p(bad_data)=NaN;
%        bad_data=find(u==-9999); u(bad_data)=NaN;
%        bad_data=find(v==-9999); v(bad_data)=NaN;
%        bad_data=find(c==-9999); c(bad_data)=NaN;
%        bad_data=find(tlt==-9999); tlt(bad_data)=NaN;
%        bad_data=find(mss==-9999); mss(bad_data)=NaN;
% 
%        eval_string(iiiRCM11(j))={['RCM11_' num2str(serialno)]};
%        
%        eval([char(eval_string(iiiRCM11(j))) '.jd=jd;']);
%        eval([char(eval_string(iiiRCM11(j))) '.t=t;']);
%        eval([char(eval_string(iiiRCM11(j))) '.p=p;']);
%        eval([char(eval_string(iiiRCM11(j))) '.u=u;']);
%        eval([char(eval_string(iiiRCM11(j))) '.v=v;']);
%        eval([char(eval_string(iiiRCM11(j))) '.c=c;']);
%        eval([char(eval_string(iiiRCM11(j))) '.tlt=tlt;']);
%        eval([char(eval_string(iiiRCM11(j))) '.mss=mss;']);
% 
%        
%        % Apply a butterworth filter to the data using auto_filt and use for
%        % plots
%        sampling_rate = 1/median(diff(jd));
% 
%        ii = eval(['find(~isnan(' char(eval_string(iiiRCM11(j))) '.u));']); 
%        eval([char(eval_string(iiiRCM11(j))) '.u(ii)=auto_filt(' char(eval_string(iiiRCM11(j)))...
%              '.u(ii), sampling_rate, 1/2,''low'',4);']);
%        ii = eval(['find(~isnan(' char(eval_string(iiiRCM11(j))) '.v));']); 
%        eval([char(eval_string(iiiRCM11(j))) '.v(ii)=auto_filt(' char(eval_string(iiiRCM11(j)))...
%              '.v(ii), sampling_rate, 1/2,''low'',4);']);
%        
%        % Calculate cumulative vectors
%        clear prog_u prog_v;
%        prog_u(1)=0;
%        prog_v(1)=0;
%        for k=2:length(u)
%            prog_u(k)=u(k-1)+prog_u(k-1);
%            prog_v(k)=v(k-1)+prog_v(k-1);
%        end
%        eval([char(eval_string(iiiRCM11(j))) '.prog_u=prog_u;']);
%        eval([char(eval_string(iiiRCM11(j))) '.prog_v=prog_v;']);
%               
%        % decimate to 1 record every ten days to use as time stamp for plot
%        b=round(sampling_rate)*10;
%        clear prog_u2; clear prog_v2; clear jd_2;
%        prog_u2=prog_u(1:b:length(prog_u));
%        prog_v2=prog_v(1:b:length(prog_v));
%        jd_2=jd(1:b:length(jd));
%        eval([char(eval_string(iiiRCM11(j))) '.prog_u2=prog_u2;']);
%        eval([char(eval_string(iiiRCM11(j))) '.prog_v2=prog_v2;']);
%        eval([char(eval_string(iiiRCM11(j))) '.jd_2=jd_2;']);
%        
%        j=j+1;
%     end
% end
% 
% %--------------------------------------
% % Now read in Argonaut data if required
% %--------------------------------------
% if iiARG>0
%     j=1;
%     
%     % loop to read one file at a time
% 
%     for i=1:length(vecARG);
%        serialno = vecARG(i);
%        disp('*************************************************************')
%        disp(['Reading ARGONAUT - ',num2str(serialno)])
%        disp('*************************************************************')
%        
%        if isunix
%            infile = [procpath,moor,'/arg/',moor,'_',num2str(vecARG(i)),'.use'];
%            if exist(infile)==0  % older Arg files had 4 digit serial number starting with zero in filename
%                infile = [procpath,moor,'/arg/',moor,'_0',num2str(vecARG(i)),'.use'];
%            end
%        elseif ispc
%            infile = [procpath,moor,'\arg\',moor,'_',num2str(vecARG(i)),'.use'];
%            if exist(infile)==0  % older Arg files had 4 digit serial number starting with zero in filename
%                infile = [procpath,moor,'\arg\',moor,'_0',num2str(vecARG(i)),'.use'];
%            end
%        end
%        
% 
%        
%        % read data into vectors and then into structure array
% 
%        % This line will need changing if have microcat data stored on the
%        % Argonaut, as per older Sontek files. The format is different.
%        [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,usd,vsd,wsd,uss,vss,wss,hdgsd,pitsd,rolsd,ipow,unoise,vnoise,wnoise,usnr,vsnr,wsnr,pgp] = ...
%            rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:usd:vsd:wsd:uss:vss:wss:hdgsd:pitsd:rolsd:ipow:unoise:vnoise:wnoise:usnr:vsnr:wsnr:pgp');
%        jd=julian(yy,mm,dd,hh);
% 
%        bad_data=find(t==-9999); t(bad_data)=NaN;
%        bad_data=find(p==-9999); p(bad_data)=NaN;
%        bad_data=find(u==-9999); u(bad_data)=NaN;
%        bad_data=find(v==-9999); v(bad_data)=NaN;
%        bad_data=find(w==-9999); w(bad_data)=NaN;
%        bad_data=find(hdg==-9999); hdg(bad_data)=NaN;
%        bad_data=find(pit==-9999); pit(bad_data)=NaN;
%        bad_data=find(rol==-9999); rol(bad_data)=NaN;
%        bad_data=find(usd==-9999); usd(bad_data)=NaN;
%        bad_data=find(vsd==-9999); vsd(bad_data)=NaN;
%        bad_data=find(wsd==-9999); wsd(bad_data)=NaN;
%        bad_data=find(uss==-9999); uss(bad_data)=NaN;
%        bad_data=find(vss==-9999); vss(bad_data)=NaN;
%        bad_data=find(wss==-9999); wss(bad_data)=NaN;
%        bad_data=find(hdgsd==-9999); hdgsd(bad_data)=NaN;
%        bad_data=find(pitsd==-9999); pitsd(bad_data)=NaN;
%        bad_data=find(rolsd==-9999); rolsd(bad_data)=NaN;
%        bad_data=find(ipow==-9999); ipow(bad_data)=NaN;
%        
%        eval_string(iiiARG(j))={['ARG_' num2str(serialno)]};
%        
%        eval([char(eval_string(iiiARG(j))) '.jd=jd;']);
%        eval([char(eval_string(iiiARG(j))) '.t=t;']);
%        eval([char(eval_string(iiiARG(j))) '.p=p;']);
%        eval([char(eval_string(iiiARG(j))) '.u=u;']);
%        eval([char(eval_string(iiiARG(j))) '.v=v;']);
%        eval([char(eval_string(iiiARG(j))) '.w=w;']);
%        eval([char(eval_string(iiiARG(j))) '.hdg=hdg;']);
%        eval([char(eval_string(iiiARG(j))) '.pit=pit;']);
%        eval([char(eval_string(iiiARG(j))) '.rol=rol;']);
%        eval([char(eval_string(iiiARG(j))) '.usd=usd;']);
%        eval([char(eval_string(iiiARG(j))) '.vsd=vsd;']);
%        eval([char(eval_string(iiiARG(j))) '.wsd=wsd;']);
%        eval([char(eval_string(iiiARG(j))) '.uss=uss;']);
%        eval([char(eval_string(iiiARG(j))) '.vss=vss;']);
%        eval([char(eval_string(iiiARG(j))) '.wss=wss;']);
%        eval([char(eval_string(iiiARG(j))) '.hdgsd=hdgsd;']);
%        eval([char(eval_string(iiiARG(j))) '.pitsd=pitsd;']);
%        eval([char(eval_string(iiiARG(j))) '.rolsd=rolsd;']);
%        eval([char(eval_string(iiiARG(j))) '.ipow=ipow;']);
%        
%        % Apply a butterworth filter to the data using auto_filt and use for
%        % plots
%        sampling_rate = 1/median(diff(jd));
% 
%        ii = eval(['find(~isnan(' char(eval_string(iiiARG(j))) '.u));']); 
%        eval([char(eval_string(iiiARG(j))) '.u(ii)=auto_filt(' char(eval_string(iiiARG(j)))...
%              '.u(ii), sampling_rate, 1/2,''low'',4);']);
%        ii = eval(['find(~isnan(' char(eval_string(iiiARG(j))) '.v));']); 
%        eval([char(eval_string(iiiARG(j))) '.v(ii)=auto_filt(' char(eval_string(iiiARG(j)))...
%              '.v(ii), sampling_rate, 1/2,''low'',4);']);
%        
%        % Calculate cumulative vectors
%        clear prog_u prog_v;
%        prog_u(1)=0;
%        prog_v(1)=0;
%        for k=2:length(u)
%            prog_u(k)=u(k-1)+prog_u(k-1);
%            prog_v(k)=v(k-1)+prog_v(k-1);
%        end
%        eval([char(eval_string(iiiARG(j))) '.prog_u=prog_u;']);
%        eval([char(eval_string(iiiARG(j))) '.prog_v=prog_v;']);
%        
%        % decimate to 1 record every ten days to use as time stamp for plot
%        b=round(sampling_rate)*10;
%        clear prog_u2; clear prog_v2; clear jd_2;
%        prog_u2=prog_u(1:b:length(prog_u));
%        prog_v2=prog_v(1:b:length(prog_v));
%        jd_2=jd(1:b:length(jd));
%        eval([char(eval_string(iiiARG(j))) '.prog_u2=prog_u2;']);
%        eval([char(eval_string(iiiARG(j))) '.prog_v2=prog_v2;']);
%        eval([char(eval_string(iiiARG(j))) '.jd_2=jd_2;']);
%               
%        j=j+1;
%     end
% end

%---------------------------------------------
% Now read in Nortek Aquadopp data if required
%---------------------------------------------
if iiNOR>0
    j=1;
    
    % loop to read one file at a time

    for i=1:length(vecNOR);
       serialno = vecNOR(i);
       disp('*************************************************************')
       disp(['Reading AQUADOPP - ',num2str(serialno)])
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

       [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = ...
           rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd');
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
       
       %remove data if there is some nan in the beginning of the timeseries
       % otherwise the cumsum did not work
       if length(u)>1
           ibad=find(isnan(u(1:10)));
           if ~isempty(ibad)
               t(1:ibad(end))=[];
               p(1:ibad(end))=[];           
               u(1:ibad(end))=[];
               v(1:ibad(end))=[];           
               w(1:ibad(end))=[];      
               jd(1:ibad(end))=[];                
           end
           %remove data outside the plot_interval    
           ibad=find(jd<jd1 | jd>jd2);
           if ~isempty(ibad)        
               u(ibad)=[];
               v(ibad)=[];             
               jd(ibad)=[];                
           end       
       end
           
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
       
       % Apply a butterworth filter to the data using auto_filt and use for
       % plots
       sampling_rate = 1/median(diff(jd));
       if length(u)>1
           ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.u));']); 
           eval([char(eval_string(iiiNOR(j))) '.u(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
                 '.u(ii), sampling_rate, 1/2,''low'',4);']);
           ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.v));']); 
           eval([char(eval_string(iiiNOR(j))) '.v(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
                 '.v(ii), sampling_rate, 1/2,''low'',4);']);
       end
         
       % Calculate cumulative vectors
       clear prog_u prog_v;
       prog_u(1)=0;
       prog_v(1)=0;
       for k=2:length(u)
           prog_u(k)=u(k-1)+prog_u(k-1);
           prog_v(k)=v(k-1)+prog_v(k-1);
       end
       eval([char(eval_string(iiiNOR(j))) '.prog_u=prog_u;']);
       eval([char(eval_string(iiiNOR(j))) '.prog_v=prog_v;']);
       
       % decimate to 1 record every ten days to use as time stamp for plot
       b=round(sampling_rate)*10;
       clear prog_u2; clear prog_v2; clear jd_2;
       if length(u)>1      
        prog_u2=prog_u(1:b:length(prog_u));
        prog_v2=prog_v(1:b:length(prog_v));
        jd_2=jd(1:b:length(jd));
       else
           prog_u2 = [];
           prog_v2 = [];
           jd_2    = [];
       end
       eval([char(eval_string(iiiNOR(j))) '.prog_u2=prog_u2;']);
       eval([char(eval_string(iiiNOR(j))) '.prog_v2=prog_v2;']);
       eval([char(eval_string(iiiNOR(j))) '.jd_2=jd_2;']);
              
       j=j+1;
    end
end
% %--------------------------------------
% % Now read in Seaguard data if required
% %--------------------------------------
% if iiSG>0
%     j=1;
%     
%     % loop to read one file at a time
% 
%     for i=1:length(vecSG);
%        serialno = vecSG(i);
%        disp('*************************************************************')
%        disp(['Reading SG - ',num2str(serialno)])
%        disp('*************************************************************')
% 
%        if isunix
%            infile = [procpath,moor,'/seaguard/',moor,'_',sprintf('%3.3d',vecSG(i)),'.use'];
%        elseif ispc
%            infile = [procpath,moor,'\seaguard\',moor,'_',sprintf('%3.3d',vecSG(i)),'.use'];
%        end
%        
% 
%        % read data into vectors and then into structure array
% 
%        [yy,mm,dd,hh,u,v,cs,cd,cssd,mss,hdg,pit,rol,t,c,tc,p,tp,ipow] = ...
%             rodbload(infile,'YY:MM:DD:HH:U:V:CS:CD:CSSD:MSS:HDG:PIT:ROL:T:C:TC:P:TP:IPOW');
%        jd=julian(yy,mm,dd,hh);
% 
%        bad_data=find(t==-9999); t(bad_data)=NaN;
%        bad_data=find(p==-9999); p(bad_data)=NaN;
%        bad_data=find(u==-9999); u(bad_data)=NaN;
%        bad_data=find(v==-9999); v(bad_data)=NaN;
%        bad_data=find(c==-9999); c(bad_data)=NaN;
%        bad_data=find(hdg==-9999); hdg(bad_data)=NaN;
%        bad_data=find(pit==-9999); pit(bad_data)=NaN;
%        bad_data=find(rol==-9999); rol(bad_data)=NaN;
% 
%        eval_string(iiiSG(j))={['SG_' num2str(serialno)]};
%        
%        eval([char(eval_string(iiiSG(j))) '.jd=jd;']);
%        eval([char(eval_string(iiiSG(j))) '.t=t;']);
%        eval([char(eval_string(iiiSG(j))) '.p=p;']);
%        eval([char(eval_string(iiiSG(j))) '.u=u;']);
%        eval([char(eval_string(iiiSG(j))) '.v=v;']);
%        eval([char(eval_string(iiiSG(j))) '.c=c;']);
%        eval([char(eval_string(iiiSG(j))) '.hdg=hdg;']);
%        eval([char(eval_string(iiiSG(j))) '.pit=pit;']);
%        eval([char(eval_string(iiiSG(j))) '.rol=rol;']);
% 
%        sampling_rate = 1/median(diff(jd));
%         
%        % Apply a butterworth filter to the data using auto_filt and use for
%            % plots
%            
%        ii = eval(['find(~isnan(' char(eval_string(iiiSG(j))) '.u));']);
%        eval([char(eval_string(iiiSG(j))) '.u(ii)=auto_filt(' char(eval_string(iiiSG(j)))...
%              '.u(ii), sampling_rate, 1/2,''low'',4);']);
%        ii = eval(['find(~isnan(' char(eval_string(iiiSG(j))) '.v));']); 
%        eval([char(eval_string(iiiSG(j))) '.v(ii)=auto_filt(' char(eval_string(iiiSG(j)))...
%              '.v(ii), sampling_rate, 1/2,''low'',4);']);
%        
%        
%        % Calculate cumulative vectors
%        clear prog_u prog_v;
%        prog_u(1)=0;
%        prog_v(1)=0;
%        for k=2:length(u)
%            prog_u(k)=u(k-1)+prog_u(k-1);
%            prog_v(k)=v(k-1)+prog_v(k-1);
%        end
%        eval([char(eval_string(iiiSG(j))) '.prog_u=prog_u;']);
%        eval([char(eval_string(iiiSG(j))) '.prog_v=prog_v;']);
%        
%        % decimate to 1 record every ten days to use as time stamp for plot
%        b=round(sampling_rate)*10;
%        clear prog_u2; clear prog_v2; clear jd_2;
%        prog_u2=prog_u(1:b:length(prog_u));
%        prog_v2=prog_v(1:b:length(prog_v));
%        jd_2=jd(1:b:length(jd));
%        eval([char(eval_string(iiiSG(j))) '.prog_u2=prog_u2;']);
%        eval([char(eval_string(iiiSG(j))) '.prog_v2=prog_v2;']);
%        eval([char(eval_string(iiiSG(j))) '.jd_2=jd_2;']);
%        
%        j=j+1;
%     end
% end
% ------------------------------
% Plotting section
% ------------------------------

% create axes

for i=1:length(iii)
    subplot_handle(iii(i))=subplot(down,across,iii(i),'align');
    eval(['plot(' char(eval_string(iii(i))) '.prog_u,' char(eval_string(iii(i))) '.prog_v,''k-'');']);
    hold on
    eval(['plot(' char(eval_string(iii(i))) '.prog_u2,' char(eval_string(iii(i))) '.prog_v2,''r.'');']);
    if length(eval([char(eval_string(iii(i))) '.prog_u2'])) > 1
        eval(['plot(' char(eval_string(iii(i))) '.prog_u2(1),' char(eval_string(iii(i))) '.prog_v2(1),''b*'');']);
    end
    ylabel('North (m)','FontSize',8);
    xlabel('East (m)','FontSize',8);
    axis equal
    set(gca,'YMinorTick','on');
    set(gca,'XMinorTick','on');
    set(gca,'FontSize',8);
    s = regexprep(eval_string(iii(i)), '_', ' ');
    
    title(strcat(s,' :',' ', num2str(depths(iii(i),2)), 'm'),'FontSize',8);
    
end

% rescale plots to match largest excursion
xlim2=[0 0]; ylim2=[0 0];
for i=1:length(iii)
    subplot(subplot_handle(i))
    xlim1=get(gca,'xlim');
    ylim1=get(gca,'ylim');
    if xlim1(2)>xlim2(2)
        xlim2(2)=xlim1(2);
    end
    if xlim1(1)<xlim2(1)
        xlim2(1)=xlim1(1);
    end
    if ylim1(2)>ylim2(2)
        ylim2(2)=ylim1(2);
    end
    if ylim1(1)<ylim2(1)
        ylim2(1)=ylim1(1);
    end
end
for i=1:length(iii)
    subplot(subplot_handle(i))
    xlim(xlim2); ylim(ylim2);
end
s2 = regexprep(moor,'_','\\_');
% suptitle is function to place a title over subplots - download from the
% Mathworks file exchange.
suptitle(strcat('Progressive Vector Plot from mooring:  ',' ',s2));
subplot(subplot_handle(1))
legend('Progession','10 day marks','start point')

print('-dpng',[moor '_pro_vec_proclvl_' proclvlstr])
savefig([moor '_pro_vec_proclvl_' proclvlstr])

