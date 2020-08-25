function plot_stacked(moor,varargin)
%
% Function for plotting  pressure, temp, conductivity, salinity, potential density for microcat.
%Composite plot by mooring.
%
% function plot_stacked('moor','proclvl','layout','plot_interval','procpath','unfiltered')
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
% 07/10/16 - Loic Houpert: cpy and adapted from currents_stacked.m

if nargin <1
    help  plot_stacked
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


%set figure size on screen for better viewing
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

%set print area of figure
pos1  = [1/8*scnsize(3),8*bdwidth,1/2*scnsize(3),(scnsize(4) - 30*bdwidth)];

plot_handle={'pressure_plot','temperature_plot','conductivity_plot','salinity_plot','pden1000_plot','ptmp1000_plot'};
ylabels={'db','deg','mS/cm','','kg/m^{-3}','deg'};
yvalue={'p','t','c','s','pden1000','ptmp1000'};
suptitle_label={'pressure','temperature','conductivity','salinity','pot. dens. rel 1000','pot. temp. rel 1000'};

for iplot=1:length(plot_handle)
        eval([ plot_handle{iplot} '=figure(' '''Position'',pos1)']);
        eval(['set(' plot_handle{iplot} ',''PaperUnits'',''centimeters'',''PaperType'', ''A4'', ''PaperOrientation'',layout);' ]);
        eval(['papersize = get('  plot_handle{iplot}  ',''PaperSize'');']);
        width=20; height=30; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
        figuresize = [left, bottom, width, height];
        eval(['set(' plot_handle{iplot} ',''PaperPosition'', figuresize);']);
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
% and find index number of AADI Seaguards
iiSG = find(id == 301);
vecSG = sn(iiSG);

% depths(:,1) = id([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG]);
% depths(:,2) = z([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG]);
depths(:,1) = id(iiMC);
depths(:,2) = z(iiMC);
depths=sortrows(depths,2);
iiiMC=find(depths(:,1)==337);
%iii=[iiiMC;iiiRBR;iiiIDR;iiiS4;iiiRCM11;iiiARG;iiiNOR;iiiSG];
iii=[iiiMC];


plot_string={};


% -----------------------------------
% START OF READING IN INSTRUMENT DATA
% -----------------------------------
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
	if proclvl==2
	       if isunix
        	   infile = [procpath,moor,'/microcat/',moor,'_',sprintf('%0.4d',vecMC(i)),'.use'];
       		elseif ispc
        	   infile = [procpath,moor,'\microcat\',moor,'_',sprintf('%0.4d',vecMC(i)),'.use'];
       		end
	elseif proclvl==3
	       if isunix
        	   infile = [procpath,moor,'/microcat/',moor,'_',sprintf('%0.3d',i),'.microcat'];
       		elseif ispc
        	   infile = [procpath,moor,'\microcat\',moor,'_',sprintf('%0.3d',i),'.microcat'];
       		end	
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
           
           c3515=sw_c3515;
           c_ratio=c/c3515;
           s=sw_salt(c_ratio,t,p); 
           pden1000=sw_pden(s,t,p,1000);
           ptmp1000=sw_ptmp(s,t,p,1000);
           
           eval_string(iiiMC(j))={['MC_' num2str(serialno)]};

           eval([char(eval_string(iiiMC(j))) '.jd=jd;']);
           eval([char(eval_string(iiiMC(j))) '.p=p;']);
           eval([char(eval_string(iiiMC(j))) '.c=c;']);    
           eval([char(eval_string(iiiMC(j))) '.t=t;']);        
           eval([char(eval_string(iiiMC(j))) '.s=s;']);               
           eval([char(eval_string(iiiMC(j))) '.pden1000=pden1000;']);
           eval([char(eval_string(iiiMC(j))) '.ptmp1000=ptmp1000;']);          

           if unfilt==0
           % Apply a butterworth filter to the data using auto_filt and use for
           % plots
           sampling_rate = 1/median(diff(jd));

            ii = eval(['find(~isnan(' char(eval_string(iiiMC(j))) '.p));']); 
            eval([char(eval_string(iiiMC(j))) '.p(ii)=auto_filt(' char(eval_string(iiiMC(j)))...
                  '.p(ii), sampling_rate, 1/2,''low'',4);']);
            ii = eval(['find(~isnan(' char(eval_string(iiiMC(j))) '.t));']); 
            eval([char(eval_string(iiiMC(j))) '.t(ii)=auto_filt(' char(eval_string(iiiMC(j)))...
                  '.t(ii), sampling_rate, 1/2,''low'',4);']);             
            ii = eval(['find(~isnan(' char(eval_string(iiiMC(j))) '.c));']); 
            eval([char(eval_string(iiiMC(j))) '.c(ii)=auto_filt(' char(eval_string(iiiMC(j)))...
                  '.c(ii), sampling_rate, 1/2,''low'',4);']);           
            ii = eval(['find(~isnan(' char(eval_string(iiiMC(j))) '.s));']); 
            eval([char(eval_string(iiiMC(j))) '.s(ii)=auto_filt(' char(eval_string(iiiMC(j)))...
                  '.s(ii), sampling_rate, 1/2,''low'',4);']);
            ii = eval(['find(~isnan(' char(eval_string(iiiMC(j))) '.pden1000));']); 
            eval([char(eval_string(iiiMC(j))) '.pden1000(ii)=auto_filt(' char(eval_string(iiiMC(j)))...
                  '.pden1000(ii), sampling_rate, 1/2,''low'',4);']);             
            ii = eval(['find(~isnan(' char(eval_string(iiiMC(j))) '.ptmp1000));']); 
            eval([char(eval_string(iiiMC(j))) '.ptmp1000(ii)=auto_filt(' char(eval_string(iiiMC(j)))...
                  '.ptmp1000(ii), sampling_rate, 1/2,''low'',4);']);           
                            
           end
        end
       j=j+1;
    end
end

% %------------------------------------
% % Now read in NORTEK data if required
% %------------------------------------
% if iiNOR>0
%     j=1;
%     
%     % loop to read one file at a time
% 
%     for i=1:length(vecNOR);
%        serialno = vecNOR(i);
%        disp('*************************************************************')
%        disp(['Reading NORTEK - ',num2str(serialno)])
%        disp('*************************************************************')
% 
% 	if proclvl==2
% 	       if isunix
%         	   infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.use'];
%        		elseif ispc
%         	   infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.use'];
%        		end
% 	elseif proclvl==3
% 	       if isunix
%         	   infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
%        		elseif ispc
%         	   infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
%        		end	
% 	end       
%        
% 
%        % read data into vectors and then into structure array
% 
%        [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd');
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
%        bad_data=find(uss==-9999); uss(bad_data)=NaN;
%        bad_data=find(vss==-9999); vss(bad_data)=NaN;
%        bad_data=find(wss==-9999); wss(bad_data)=NaN;
%        bad_data=find(ipow==-9999); ipow(bad_data)=NaN;
%        bad_data=find(cs==-9999); cs(bad_data)=NaN;
%        bad_data=find(cd==-9999); cd(bad_data)=NaN;
% 
%        eval_string(iiiNOR(j))={['NOR_' num2str(serialno)]};
%        
%        eval([char(eval_string(iiiNOR(j))) '.jd=jd;']);
%        eval([char(eval_string(iiiNOR(j))) '.t=t;']);
%        eval([char(eval_string(iiiNOR(j))) '.p=p;']);
%            eval([char(eval_string(iiiNOR(j))) '.c=t*nan;']);           
%            eval([char(eval_string(iiiNOR(j))) '.s=t*nan;']);               
%            eval([char(eval_string(iiiNOR(j))) '.pden1000=t*nan;']);
%            eval([char(eval_string(iiiNOR(j))) '.ptmp1000=t*nan;']);             
%            if unfilt==0
%                % Apply a butterworth filter to the data using auto_filt and use for
%                % plots
%                sampling_rate = 1/median(diff(jd));
% 
%                 ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.p));']); 
%                 eval([char(eval_string(iiiNOR(j))) '.p(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
%                       '.p(ii), sampling_rate, 1/2,''low'',4);']);
%                 ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.t));']); 
%                 eval([char(eval_string(iiiNOR(j))) '.t(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
%                       '.t(ii), sampling_rate, 1/2,''low'',4);']);                 
%            end
%                
%        j=j+1;
%     end
% end


% ------------------------------
% Plotting section
% ------------------------------
% calculate limits of u and v axes using MAX_u and MAX_v
emptycellcheck=cellfun('isempty',eval_string);
eval_string(emptycellcheck)=[];    
    
for k=1:length(plot_handle)
    eval(['figure(' plot_handle{k} ')']);
    
    % setup subplot number of panels
    subplot(length(depths(:,1)),1,1);

    % create axes
    for i=1:length(eval_string)
        % create axes
        axes_string=['axes' num2str(i)];
        eval([axes_string '=subplot(length(depths(:,1)),1,i);']);
        
        % plot data
        eval(['plot((' eval_string{i} '.jd-jd1),' eval_string{i} '.' yvalue{k} ');']);
        
        eval(['ylabel(''' ylabels{k} ''');']);
        xlim([0 jd2-jd1]);
        set(gca,'YMinorTick','on');
        set(gca,'xTickLabel',xticklabels);
        set(gca,'XTick',jdxticks-jd1);
        
        Y_limits=get(gca,'ylim');
%         Y_limits=(ylimits{k});
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
        s = regexprep(eval_string(i), '_', '\nsn:');
        text((X_limits(2)-X_limits(1))*1.005+X_limits(1), (Y_limits(2)-Y_limits(1))*0.75+Y_limits(1),s,'FontSize',8);

        text((X_limits(2)-X_limits(1))*1.005+X_limits(1), (Y_limits(2)-Y_limits(1))*0.5+Y_limits(1),[num2str(depths(i,2)) 'm'],'FontSize',8);

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
    
    print('-dpng',[moor '_plot_stacked_' plot_handle{k} '_proclvl_' proclvlstr])
    savefig([moor '_plot_stacked_' plot_handle{k} '_proclvl_' proclvlstr])

end
