% function overlay_plot_severalyr(moorlist,vartoplot,'procpath','proclvl','layout','plot_interval','unfiltered')
%
% Function for plotting pressure, temp, conductivity, salinity, potential density of a mooring overlayed on the same axes
%
% required inputs:-
%   moorlist: list of complete mooring name as string. e.g.
%   {'nocm5_01_2014','nocm5_02_2015'}
%   vartoplot: variable to plot
%
% optional inputs:-
%   layout: orientation of figure portrait/lanscape (default = portrait)
%           input of 'landscape' or 'portrait'
%           e.g. pressure_overlay('wb1_1_200420','layout','landscape')
%   plot_interval: matrix of start and end dates for plot
%           e.g. pressure_overlay('wb1_1_200420','plot_interval',[2004 02 01 00; 2005 06 01 00])
%           dates are:- yyyy mm dd hh. Default is calculated automatically
%   procpath: can specify exact procpath if not using standard data paths. 
%           e.g. pressure_overlay('wb1_1_200420','inpath','/Volumes/noc/mpoc/hydro/rpdmoc/rapid/data/moor/proc/')
%   proclvl: can specify level of processing of the data to plot. 
%           e.g. 'proclvl','2': will plot the .use file ; 'proclvl','3' will plot the .microcat and .edt files
%   unfiltered: plot data in unfiltered format - input 'unfiltered'
%
% functions called:-
%   rodbload, julian, auto_filt
%   from .../exec/moor/tools and .../exec/moor/rodb paths
% 
% Routine written by Darren Rayner July 2006.

% 09/01/2008 - added Nortek Aquadopp capability
% 09/01/2008 - added PC compatability
% 15/04/2011 - aboard KN200-4: added SeaGuard capability
% 05/10/16 - Loic Houpert: added option to process lvl 3 data (.microcat and .edt files for nortek) and save plot
% 06/10/16 - copy from pressure_overlay to plot several deployment on the
% same plot and adapt to plot different variable
%
function overlay_plot_severalyr(moorlist,varargin)
if nargin <1
    help overlay_plot_severalyr
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
disp('The proc path needed to be specify')
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
        width=17; height=26; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
        figuresize = [left, bottom, width, height];
        eval(['set(' plot_handle{iplot} ',''PaperPosition'', figuresize);']);
end

timemax=0;
timemin=999999999999999;
% Need to be define
% xlim(plotxlim); %[0 jd2-jd1]);
% set(gca,'xTickLabel',plotxticklabels);%xticklabels);
% set(gca,'XTick',plotxtick);%jdxticks-jd1);
for imoor =1:length(moorlist)

    moor = moorlist{imoor};
    
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

disp('z : instrument id : serial number')
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


jd1 = julian(plot_interval(1,:));
if julian(plot_interval(1,:))< timemin
    timemin = julian(plot_interval(1,:)); %jd2
end
if julian(plot_interval(2,:)) >timemax
    timemax = julian(plot_interval(2,:)); %jd2
end

num_to_plot=2; % num_to_plot is the number of samples per day to plot and can be adjusted accordingly

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

depths(:,1) = id([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG]);
depths(:,2) = z([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG]);
depths=sortrows(depths,2);
iiiMC=find(depths(:,1)==337);
iiiRBR=find(depths(:,1)==330);
iiiIDR=find(depths(:,1)==339);
iiiS4=find(depths(:,1)==302);  
iiiRCM11=find(depths(:,1)==310);  
iiiARG=find(depths(:,1)==366 | depths(:,1)==366337); 
iiiNOR=find(depths(:,1)==368|depths(:,1)==370);
iiiSG=find(depths(:,1)==301);
iii=[iiiMC;iiiRBR;iiiIDR;iiiS4;iiiRCM11;iiiARG;iiiNOR;iiiSG];


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
       
       % check if file exists
       fileopen=fopen(infile,'r');
      
       if fileopen>0
           % read data into vectors and then into structure array

           [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = ...
               rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd');
           jd=julian(yy,mm,dd,hh);

           bad_data=find(p==-9999); p(bad_data)=NaN;
           bad_data=find(t==-9999); t(bad_data)=NaN;
           
           eval_string(iiiNOR(j))={['NOR_' num2str(serialno)]};

           eval([char(eval_string(iiiNOR(j))) '.jd=jd;']);
           eval([char(eval_string(iiiNOR(j))) '.p=p;']);
           eval([char(eval_string(iiiNOR(j))) '.t=t;']);       
           eval([char(eval_string(iiiNOR(j))) '.c=t*nan;']);           
           eval([char(eval_string(iiiNOR(j))) '.s=t*nan;']);               
           eval([char(eval_string(iiiNOR(j))) '.pden1000=t*nan;']);
           eval([char(eval_string(iiiNOR(j))) '.ptmp1000=t*nan;']);             
           if unfilt==0
               % Apply a butterworth filter to the data using auto_filt and use for
               % plots
               sampling_rate = 1/median(diff(jd));

                ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.p));']); 
                eval([char(eval_string(iiiNOR(j))) '.p(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
                      '.p(ii), sampling_rate, 1/2,''low'',4);']);
                ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.t));']); 
                eval([char(eval_string(iiiNOR(j))) '.t(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
                      '.t(ii), sampling_rate, 1/2,''low'',4);']);                 
           end
       end
       
       j=j+1;
    end
end


% ------------------------------
% Plotting section
% ------------------------------
colours = 'bkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbk';

emptycellcheck=cellfun('isempty',eval_string);
eval_string(emptycellcheck)=[];    
    
for iplot=1:length(plot_handle)
for i=1:length(eval_string)
    eval(['figure(' plot_handle{iplot} ')']); hold on
    eval(['plot(' char(eval_string(i)) '.jd-timemin,' char(eval_string(i)) '.' yvalue{iplot} ',''' colours(i) ''');']);
end


% Display year labels on bottom graph
Y_limits=ylim;   X_limits=xlim;
% % label data with nominal instrument depth
if imoor==1
    label_x_positions=X_limits(1);   
else
    label_x_positions=X_limits(2)*0.97; %(X_limits(2)-X_limits(1))*1.005+X_limits(1);
end
% % line below is opposite to currents_stacked and temperature_overlay etc
% % because y axis is inverted in this pressure_overlay function.
% a=(Y_limits(2)-Y_limits(1))*1.05+Y_limits(1);
% text(X_limits(1),a,num2str(xticks(1,1)),'FontSize',10);
% for i=1:length(year_indexes)
%     text((timemax-timemin)*(year_indexes(i)-1)/(length(xticklabels)-1),a,num2str(xticks(year_indexes(i),1)),'FontSize',10);
% end


for i=1:length(eval_string)
    if imoor==1
        label_y_positions=eval([ char(eval_string{i})  '.' yvalue{iplot} '(1)']); %label_y_positions=eval(['nanmedian(' char(eval_string{i})  '.' yvalue{iplot} ');']); 
    else
        label_y_positions=eval([ char(eval_string{i})  '.' yvalue{iplot} '(end);']);        
    end
    eval(['text(label_x_positions, label_y_positions,[num2str(depths(i,2)) ''m''],''FontSize'',8,''color'',''' colours(i) ''');'])
end

plot_interval=0;

end
end


for iplot=1:length(plot_handle)
    eval(['figure(' plot_handle{iplot} ')']);
    plotxlim = [0 timemax-timemin];
    plotxtick = timemin:round((timemax-timemin)/6):timemax;
    plotxticklabels = datestr(datenum(gregorian(plotxtick)),'mmm yy');
    ylabel(ylabels{iplot});
    xlim(plotxlim); %[0 jd2-jd1]);
%     ylimits=get(gca,'ylim');
%     %ylimits(1)=0;
%     set(gca,'ylim',ylimits)
    set(gca,'YMinorTick','on');
    if ~isempty(strfind(plot_handle{iplot},'pressure'))
        set(gca,'YDir','reverse');
    end
    set(gca,'xTickLabel',plotxticklabels);%xticklabels);
    set(gca,'XTick',plotxtick-timemin);%jdxticks-jd1);
    s = regexp(moor,'_','split');

    if unfilt==0
        title(['Low-pass filtered ' suptitle_label{iplot} ' from mooring: ' s{1}])
    else
        title(['Unfiltered ' suptitle_label{iplot} ' from mooring: ' s{1}])
    end


    %keyboard
    print('-dpng',[s{1} '_' plot_handle{iplot} '_overlay_proclvl_' proclvlstr])
    savefig([s{1} '_'  plot_handle{iplot} '_overlay_proclvl_' proclvlstr])
end

close all

