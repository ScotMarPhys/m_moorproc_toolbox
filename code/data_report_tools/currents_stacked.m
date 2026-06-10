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
% 08/02/26 - Y Firing put repeated code into function, reduce use of eval

if nargin <1
    help currents_stacked
    return
end

layout = 'portrait';
ii = find(strcmp(varargin,'layout')); 
if ~isempty(ii); layout = varargin{ii+1}; end
proclvl = 2;
ii = find(strcmp(varargin,'proclvl')); 
if ~isempty(ii); proclvl = varargin{ii+1}; end
plot_interval = 0;
ii = find(strcmp(varargin,'plot_interval')); 
if ~isempty(ii); plot_interval = varargin{ii+1}; end
unfilt = 0;
ii = find(strcmp(varargin,'unfiltered')); 
if ~isempty(ii); unfilt = varargin{ii+1}; end

% Load vectors of mooring information
% id instrument id, sn serial number, z nominal depth of each instrument
% s_t, e_t, s_d, e_d start and end times and dates
% lat lon mooring position, wd corrected water depth (m)
% mr mooring name
pd = moor_inoutpaths('reports',moor);
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(pd.infofile,...
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


%get information about all the instruments
id_z_sn = all_inst_table(id, z, sn);
%these rows have variables defined for currents_stacked, keep them
    mc = cellfun(@(x) length(x)>1, id_z_sn.vars_cs); 
if ~sum(mc)
    warning('no current data from %s; skipping',moor)
    return
end
id_z_sn = id_z_sn(mc,:);
ninst = length(id_z_sn.id);
%sort by z
id_z_sn = sortrows(id_z_sn,2); 
disp(['z : instrument id : serial number for stacked currents plots'])
disp([id_z_sn(:,[2 1 3])])

maxes.u = 0; maxes.v = 0; maxes.www = 0; maxes.spd = 0;
dummy = -9999;

%load data
for ino = 1:ninst

    serialno = id_z_sn.sn(ino);
    
    disp('*************************************************************')
    disp(['Reading ' id_z_sn.inst{ino} ' - ' num2str(serialno)])
    disp('*************************************************************')
    
    pd = moor_inoutpaths(id_z_sn.instl{ino}, moor);
    infile = fullfile(pd.stage2path, sprintf(pd.stage2form, serialno));
    [dat(ino).data, maxes] = currents_load_calc(infile, id_z_sn.vars_cs{ino}, dummy, unfilt, maxes);
    
end


% ------------------------------
% Plotting section
% ------------------------------

% ------------------------------------------------
% Determine plot_interval if not input to function
if nargin == 1 || (exist('plot_interval','var') && plot_interval==0)
    plot_interval = zeros(2,4);
    plot_interval(1,:) = [s_d(1:2)' 1 0];
    plot_interval(2,:) = [e_d(1) e_d(2)+1 1 0];
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

% calculate limits of u and v axes using maxes.u and maxes.v
maxes.u=ceil(maxes.u/10)*10;
maxes.v=ceil(maxes.v/10)*10;
% calculate limits of spd axes using maxes.spd
maxes.spd=ceil(maxes.spd/10)*10;

%set figure size on screen for better viewing
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

%set print area of figure
pos1 = [1/8*scnsize(3),8*bdwidth,1/2*scnsize(3),(scnsize(4) - 30*bdwidth)];
plot_handle={'horspeed_plot','hordirection_plot','u_plot','v_plot','w_plot'};
ylabels={'speed (cm/s)','direction (deg)','velocity (cm/s)','velocity (cm/s)','velocity (cm/s)'};
yvalue={'spd', 'dir', 'u', 'v', 'w'};
suptitle_label={'hor. current speed','hor. current direction','u-velocity component','v-velocity component','w-velocity component'};
ylimits={[0 maxes.spd],[0 360],[-maxes.u maxes.u],[-maxes.v maxes.v],[-maxes.www maxes.www]};

for no = 1:length(plot_handle)
    ph.(plot_handle{no})=figure('Position',pos1);
    set(ph.(plot_handle{no}),'PaperUnits','centimeters');
    set(ph.(plot_handle{no}), 'PaperType', 'A4');
    set(ph.(plot_handle{no}), 'PaperOrientation',layout);
    if no==1
        papersize = get(ph.(plot_handle{1}),'PaperSize');
        width=17; height=26; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
        figuresize = [left, bottom, width, height];
    end
    set(ph.(plot_handle{no}), 'PaperPosition', figuresize);
end

for k=1:length(plot_handle)
    figure(ph.(plot_handle{k}))
    % setup subplot number of panels
    subplot(ninst,1,1);

    % create axes
    clear ax
    for ino=1:ninst
        % create axes
        ax(ino) = subplot(ninst,1,ino);
        
        % plot data
        
        plot(dat(ino).data.jd-jd1, dat(ino).data.(yvalue{k}));
        ylabel(ylabels{k});
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
        s = ['S/N ' num2str(id_z_sn.sn(ino)) ', ' num2str(id_z_sn.z(ino)) ' m'];
        text((X_limits(2)-X_limits(1))*1.005+X_limits(1), (Y_limits(2)-Y_limits(1))*0.5+Y_limits(1),s,'FontSize',10,'Rotation',90);

    end

    % Display year labels on bottom graph
    axes(ax(ino))
    a=(Y_limits(1)-Y_limits(2))*1.5+Y_limits(2);
    text(X_limits(1),a,num2str(xticks(1,1)),'FontSize',10);
    for yno=1:length(year_indexes)
        text((jd2-jd1)*(year_indexes(yno)-1)/(length(xticklabels)-1),a,num2str(xticks(year_indexes(yno),1)),'FontSize',10);
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
    
    figname = [moor '_currents_stacked_' plot_handle{k} '_proclvl_' num2str(proclvl)];
    print('-dpng',figname)
    savefig(figname)

end


function [datt, maxes] = currents_load_calc(infile, rodb_varstr, dummy, unfilt, maxes)

datc = rodbload(infile,rodb_varstr);
datt = rodb2masked_table(datc, rodb_varstr, dummy);

if unfilt==0
    %replace u and v with filtered versions***this will average
    %across gaps?***
    sampling_rate = 1/median(diff(datt.jd));
    if ~isempty(datt.u)
        ii = find(~isnan(datt.u));
        datt.u(ii) = auto_filt(datt.u(ii),sampling_rate,1/2,'low',4);
        ii = find(~isnan(datt.v));
        datt.v(ii) = auto_filt(datt.v(ii),sampling_rate,1/2,'low',4);
    end
end
% determine current_speed limits for setting y-axis after plotting
maxes.u = max(maxes.u,max(abs(datt.u)));
maxes.v = max(maxes.v,max(abs(datt.v)));

% calculate speed and direction
datt.spd = sqrt(datt.u.^2 + datt.v.^2);
maxes.spd=max(maxes.spd,max(datt.spd));
if ismember('w',datt.Properties.VariableNames)
    datt.www = datt.w;
    maxes.www = max(maxes.www,max(abs(datt.www)));
end
datt.dir = atan(datt.u./datt.v)*180/pi;
m = datt.v<0; datt.dir(m) = datt.dir(m)+180;
m = datt.v>=0 & datt.u<0; datt.dir(m) = datt.dir(m)+360;

