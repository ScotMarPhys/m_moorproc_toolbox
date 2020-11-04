% function currents_stacked_severalyr(moorlist,vartoplot,'procpath','proclvl','layout','plot_interval','unfiltered')
% 
% Function for plotting u,v, speed and direction plots from RCM11s, S4s and
% Argonauts. Composite plot by mooring.
%
% function currents_stacked('moor','proclvl','layout','plot_interval','procpath','unfiltered')
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
% Routine written by Darren Rayner January 2008 - adapted from stick_plot.
%
% 15/4/2011 - aboard cruise kn200-4: added Seaguard functionality 
% 05/10/16 - Loic Houpert: added option to process lvl 3 data (.microcat and .edt files for nortek) and save plot
% 12/02/18 - Loic Houpert: created and adapted from
% overlay_plot_severalyr.m

function currents_stacked_severalyr(moorlist,varargin)
if nargin <1
    help  currents_stacked_severalyr
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

plot_handle={'horspeed_plot','hordirection_plot','u_plot','v_plot','w_plot'};
ylabels={'speed (cm/s)','direction (deg)','velocity (cm/s)','velocity (cm/s)','velocity (cm/s)'};
yvalue={'spd', 'dir', 'u', 'v', 'w'};
suptitle_label={'hor. current speed','hor. current direction','u-velocity component','v-velocity component','w-velocity component'};
vertshift ={'50','400','40','40','15'};

for iplot=1:length(plot_handle)
        eval([ plot_handle{iplot} '=figure(' '''Position'',pos1)']);
        eval(['set(' plot_handle{iplot} ',''PaperUnits'',''centimeters'',''PaperType'', ''A4'', ''PaperOrientation'',layout);' ]);
        eval(['papersize = get('  plot_handle{iplot}  ',''PaperSize'');']);
        width=17; height=26; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
        figuresize = [left, bottom, width, height];
        eval(['set(' plot_handle{iplot} ',''PaperPosition'', figuresize);']);
end



MAX_u=0;
MAX_v=0;
MAX_www=0;
MAX_spd=0;

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

plot_string={};

for iplot=1:length(plot_handle)
    for i=1:length(iii)
    eval(['figure(' plot_handle{iplot} ')']); hold on           
    eval(['plot([0 timemax-timemin],[0 timemax-timemin]*0 - (i-1)*' vertshift{iplot}  ',''--k'');']); 
    end
end

% -----------------------------------
% START OF READING IN INSTRUMENT DATA
% -----------------------------------


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
           if unfilt==1      
               if isunix
                   infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
                elseif ispc
                   infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
               end
           end                
           if unfilt==0
               if isunix
                   infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.lowedt'];
               elseif ispc
                   infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.lowedt'];
               end	       
           end        
    end
       
       % check if file exists
       fileopen=fopen(infile,'r');
      
       if fileopen>0
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


%            if unfilt==0
% 
%            uf=[];
%            vf=[];          
%                if ~isempty(eval([ char(eval_string(iiiNOR(j))) '.u']))
%                ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.u));']); 
%                eval(['uf=auto_filt(' char(eval_string(iiiNOR(j)))...
%                      '.u(ii), sampling_rate, 1/2,''low'',4);']);
% 
%                ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.v));']); 
%                eval(['vf=auto_filt(' char(eval_string(iiiNOR(j)))...
%                      '.v(ii), sampling_rate, 1/2,''low'',4);']);           
%                    eval([char(eval_string(iiiNOR(j))) '.u(ii)=uf;']);
%                    eval([char(eval_string(iiiNOR(j))) '.v(ii)=vf;']);
%                else
%                end
%            else
%                % do not reassign u and v with filtered versions
%            end

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

               max_u=max(abs(u));
               max_v=max(abs(v));

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


% ------------------------------
% Plotting section
% ------------------------------
colours = 'bkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbkbk';

emptycellcheck=cellfun('isempty',eval_string);
eval_string(emptycellcheck)=[];    
    
for iplot=1:length(plot_handle)
for i=1:length(eval_string)
    eval(['figure(' plot_handle{iplot} ')']); hold on   
    eval(['plot(' char(eval_string(i)) '.jd-timemin,' char(eval_string(i)) '.' yvalue{iplot} ' - (i-1)*' vertshift{iplot}  ',''' colours(i) ''');']);
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
    eval(['text(label_x_positions, label_y_positions  - (i-1)*' vertshift{iplot}  ',[num2str(depths(i,2)) ''m''],''FontSize'',8,''color'',''' colours(i) ''');'])
end

plot_interval=0;

end
end


% ylimits={[-MAX_u MAX_u],[-MAX_v MAX_v],[-MAX_www MAX_www],[0 MAX_spd],[0 360]};

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
    print('-dpng',[s{1} '_' plot_handle{iplot} '_stacked_proclvl_' proclvlstr])
    savefig([s{1} '_'  plot_handle{iplot} '_stacked_proclvl_' proclvlstr])
end

close all



