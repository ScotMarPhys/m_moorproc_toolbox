
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

global MOORPROC_G

if nargin <1
    help currents_overlay
    return
end

%defaults
layout = 'portrait';
procpath = fullfile(MOORPROC_G.moordatadir,'proc');
outpath = fullfile(MOORPROC_G.reportdir,'figs');
plot_interval = 0;
width=26; height=17;
num_to_plot = 2;
unfilt = 0;
%and optional inputs overwrite them
n = 1;
while n<=nargin-1
    if ischar(varargin{n})
        if strcmp('unfiltered',varargin{n})
            unfilt = 1;
            n = n+1;
        else
            eval([varargin{n} ' = varargin{n+1};'])
            n = n+2;
        end
    end
end

if ~exist('proclvl','var')
    proclvl = 2;
elseif ischar(proclvl)
    proclvl = str2double(proclvl);
end
proclvlstr0 = num2str(proclvl);
if unfilt
    proclvlstr = [proclvlstr0 '_unfilt'];
else
    proclvlstr = [proclvlstr0 '_lpfilt'];
end

infofile = fullfile(procpath,moor,[moor 'info.dat']);

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

MAX_spd=0;
MAX_www=0;


%get information about all the instruments including prefixes and filenames/paths, in a table we can loop through
id_z_sn = all_inst_table(id, z, sn);
id_z_sn.data_loaded = false(length(id),1);

for iid=1:length(id_z_sn.id)
    if ~isempty(id_z_sn.dirs{iid}) && contains(id_z_sn.vars{iid},'u')
        disp('*************************************************************')
        disp(['Reading ' id_z_sn.inst{iid} ' - ',num2str(id_z_sn.sn(iid))])
        disp('*************************************************************')
        iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
        infile = fullfile(procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%0.4d%s.use',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
        if ~exist(infile,'file')
            infile = fullfile(procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%3.3d%s.use',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
        end
        %read data into structure array
        fileopen=fopen(infile,'r');
        if fileopen>0
                        varstr = ['yy:mm:dd:hh:' id_z_sn.vars{iid}];
            clear a; a(1).a = [];
            [yy,mm,dd,hh,a(1).a,a(2).a,a(3).a,a(4).a,a(5).a,a(6).a,a(7).a,a(8).a,a(9).a] = rodbload(infile,varstr);
            vars = split(id_z_sn.vars{iid},':');
            for no = 1:length(vars)
                data.(vars{no}) = a(no).a;
                data.(vars{no})(data.(vars{no})==-9999) = NaN;
            end
            data.jd=julian(yy,mm,dd,hh);

        else
            disp('File does not exist!')
            disp(['infile = ' infile])
            data = [];
        end

        if ~isempty(data)
            if unfilt==0
                sampling_rate = 1/median(diff(data.jd));
                ii = find(~isnan(data.u));
                data.u(ii) = auto_filt(data.u(ii),sampling_rate,1/2,'low',4);
                ii = find(~isnan(data.v));
                data.v(ii) = auto_filt(data.v(ii),sampling_rate,1/2,'low',4);
                if isfield(data,'w')
                    ii = find(~isnan(data.w));
                    data.w(ii) = auto_filt(data.w(ii),sampling_rate,1/2,'low',4);
                end
            end
            data.spd = sqrt(data.u.^2+data.v.^2);
            data.dir = atan2(data.v,data.u)*180/pi;
            data.dir(data.dir<0) = data.dir(data.dir<0)+360;
            %current_speed for setting y-axis after plotting
            MAX_spd=max(max(data.spd),MAX_spd);
            if isfield(data,'w')
                data.www = data.w;
                MAX_www = max(max(data.www),MAX_www);
            end

            alldata.(iname) = data;
            id_z_sn.data_loaded(iid) = true;
        end
    end

end
id_z_sn = id_z_sn(id_z_sn.data_loaded,:);

% ------------------------------
% Plotting section
% ------------------------------
colours = 'brgcmybrgcmybrgcmy';
for iid = 1:length(id_z_sn.id)
    if iid<7 
        line_style = '-';
    elseif iid>12
        line_style = '.';
    else
        line_style = 'o';
    end
    a=1;%sampling_rate/num_to_plot; %value to decimate data by to ease plot visualisation doesn't actually average it, just decimates
    % if sampling rate equals num_to_plot entered in function then no decimation occurs.
    % Assumes same sampling rate for all current meters on the mooring.
    %a = 5; 
    iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
    b = length(alldata.(iname).jd);
    figure(direction_plot); hold on
    plot(alldata.(iname).jd(1:a:b)-jd1,alldata.(iname).dir(1:a:b),colours(iid))
    figure(speed_plot); hold on
    plot(alldata.(iname).jd(1:a:b)-jd1,alldata.(iname).spd(1:a:b),colours(iid))
    figure(www_plot); hold on
    if isfield(alldata.(iname),'www')
        plot(alldata.(iname).jd(1:a:b)-jd1,alldata.(iname).www(1:a:b),colours(iid))
    end
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
if unfilt==0
    title(['Current directions of low pass filtered currents at mooring ' s '.'])
else
    title(['Current directions of unfiltered currents at mooring ' s '.'])
end
for iid = 1:length(id_z_sn.id)
    iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
    legend_text(iid)={[num2str(id_z_sn.z(iid)) ' m']};
end
legend(legend_text);
text(1,-20,num2str(xticks(1,1)),'FontSize',10);

if plot_x_labels>0
    for i=1:length(year_indexes)
        text((jd2-jd1)*(year_indexes(i)-1)/(length(xticklabels)-1),-20,num2str(xticks(year_indexes(i),1)),'FontSize',10);
    end
end
outfile = fullfile(outpath,[moor '_horcurrents_overlay_dir_proclvl_' proclvlstr]);
print('-dpng',outfile)
savefig(outfile)

figure(speed_plot);
ylabel('current speed (cm/s)');
xlim([0 jd2-jd1]);
MAX_spd=ceil(MAX_spd/10)*10;
ylim([0 MAX_spd]);
set(gca,'YMinorTick','on');
set(gca,'xTickLabel',xticklabels);
set(gca,'XTick',jdxticks-jd1);
s = regexprep(moor,'_','\\_');
if unfilt==0
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

outfile = fullfile(outpath,[moor '_horcurrents_overlay_speed_proclvl_' proclvlstr]);
print('-dpng',outfile)
savefig(outfile)

figure(www_plot);
ylabel('current speed (cm/s)');
xlim([0 jd2-jd1]);
MAX_www=ceil(MAX_www/10)*10;
ylim([0 MAX_www]);
set(gca,'YMinorTick','on');
set(gca,'xTickLabel',xticklabels);
set(gca,'XTick',jdxticks-jd1);
s = regexprep(moor,'_','\\_');
if unfilt==0
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

outfile = fullfile(outpath,[moor '_vertcurrents_overlay_speed_proclvl_' proclvlstr]);
print('-dpng',outfile)
savefig(outfile)
