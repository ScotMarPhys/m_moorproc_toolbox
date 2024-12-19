function overlay_plots(moor,datatype,varargin)
%
% Function for plotting various types of data from a mooring, with data
% from different instruments overlaid on the same axis
%
% function overlay_plots('moor','currents','proclvl','layout','plot_interval','procpath')
%    plots 3 subplots: speed, direction, vertical speed
% function overlay_plots('moor','properties','proclvl','layout','plot_interval','procpath')
%    plots 3-4 subplots: pressure, temperature, salinity, (oxygen)
%
% required inputs:-
%   moor: complete mooring name as string. e.g. 'wb1_1_200420'
%   datatype: 'currents' or 'properties'
%
% optional inputs (otherwise set by plot_options.m):-
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
close all

if nargin <1
    help currents_overlay
    return
end
if nargin>1 && isstruct(varargin{1})
    plotpar = varargin{1};
else
    inargs = varargin; plot_options
end
if ~exist(plotpar.outpath,'dir')
    warning('making directory %s',plotpar.outpath)
    mkdir(plotpar.outpath)
end

% Load vectors of mooring information
% id instrument id, sn serial number, z nominal depth of each instrument
% s_t, e_t, s_d, e_d start and end times and dates
% lat lon mooring position, wd corrected water depth (m)
% mr mooring name
[id,sn,z,s_t,s_d,e_t,e_d]  =  rodbload(plotpar.infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date');
disp(['z : instrument id : serial number'])
for i = 1:length(id)
    disp([z(i),id(i),sn(i)])
end
plotpar.times = [s_d(1:2)'; e_d(1:2)'];

%initialise plots
plots = set_plots(datatype,moor,plotpar);

%get information about all the instruments including prefixes and filenames/paths, in a table we can loop through
id_z_sn = all_inst_table(id, z, sn);

%load data
[alldata, id_z_sn, plots] = load_alldata(moor, id_z_sn, plots, plotpar, datatype);
id_z_sn = id_z_sn(id_z_sn.data_loaded,:);

%make plots
plots = make_plots(plots, id_z_sn, alldata, plotpar);


% ------------------------------
% Setup section
% ------------------------------
function plots = set_plots(datatype,moor,plotpar)

s = regexprep(moor,'_','\\_');

switch datatype
    case 'currents'
        plots.direction.h = [];
        plots.direction.var = 'dir';
        plots.direction.ylabel = 'current direction (deg M)';
        plots.direction.title = ['Current directions of ' plotpar.proclvltext ' currents at mooring ' s];
        plots.direction.outfile = fullfile(plotpar.outpath,[moor '_horcurrents_overlay_dir_proclvl_' plotpar.proclvlstr]);
        plots.direction.ytick = [0 90 180 270 360];

        plots.speed.h = [];
        plots.speed.var = 'spd';
        plots.speed.ylabel = 'current speed (cm/s)';
        plots.speed.title = ['Hor. current speed of ' plotpar.proclvltext ' currents at mooring ' s];
        plots.speed.outfile = fullfile(plotpar.outpath,[moor '_horcurrents_overlay_speed_proclvl_' plotpar.proclvlstr]);
        plots.speed.ylim = [0 0];

        plots.www.h = [];
        plots.www.var = 'w';
        plots.www.ylabel = 'current speed (cm/s)';
        plots.www.title = ['Vert. current speed of ' plotpar.proclvltext ' currents at mooring ' s];
        plots.www.outfile = fullfile(plotpar.outpath,[moor '_vertcurrents_overlay_speed_proclvl_' plotpar.proclvlstr]);
        plots.www.ylim = [0 0];

    case 'properties'

        plots.pressure.h = [];
        plots.pressure.var = 't';
        plots.pressure.ylabel = 'pressure (dbar)';
        plots.pressure.title = [upper(plotpar.proclvltext(1)) plotpar.proclvltext(2:end) ' pressure from mooring: ' s];
        plots.pressure.outfile = fullfile(plotpar.outpath,[moor '_pressure_overlay_proclvl_' plotpar.proclvlstr]);

        plots.temperature.h = [];
        plots.temperature.var = 't';
        plots.temperature.ylabel = 'temperature (degC)';
        plots.temperature.title = [upper(plotpar.proclvltext(1)) plotpar.proclvltext(2:end) ' temperature from mooring: ' s];
        plots.temperature.outfile = fullfile(plotpar.outpath,[moor '_temperature_overlay_proclvl_' plotpar.proclvlstr]);

        plots.salinity.h = [];
        plots.salinity.var = 's';
        plots.salinity.ylabel = 'salinity (psu)';
        plots.salinity.title = [upper(plotpar.proclvltext(1)) plotpar.proclvltext(2:end) ' salinity from mooring: ' s];
        plots.salinity.outfile = fullfile(plotpar.outpath,[moor '_salinity_overlay_proclvl_' plotpar.proclvlstr]);

        plots.conductivity.h = [];
        plots.conductivity.var = 'c';
        plots.conductivity.ylabel = 'conductivity (mS/cm)';
        plots.conductivity.title = [upper(plotpar.proclvltext(1)) plotpar.proclvltext(2:end) ' conductivity from mooring: ' s];
        plots.conductivity.outfile = fullfile(plotpar.outpath,[moor '_conductivity_overlay_proclvl_' plotpar.proclvlstr]);

        plots.oxygen.h = [];
        plots.oxygen.var = 'o2';
        plots.oxygen.ylabel = 'oxygen (umol/kg)';
        plots.oxygen.title = [upper(plotpar.proclvltext(1)) plotpar.proclvltext(2:end) ' oxygen from mooring: ' s];
        plots.oxygen.outfile = fullfile(plotpar.outpath,[moor '_oxygen_overlay_proclvl_' plotpar.proclvlstr]);


end

% ------------------------------
% Loading section
% ------------------------------
function [alldata, id_z_sn, plots] = load_alldata(moor, id_z_sn, plots, plotpar, datatype)

id_z_sn.data_loaded = false(length(id_z_sn.id),1);

switch datatype
    case 'currents'
        diagntype = 'u';
    case 'properties'
        diagntype = 't';
end

for iid=1:length(id_z_sn.id)
    if ~isempty(id_z_sn.dirs{iid}) && contains(id_z_sn.vars{iid},diagntype)
        if ~plotpar.non_verbose
            disp('*************************************************************')
            disp(['Reading ' id_z_sn.inst{iid} ' - ',num2str(id_z_sn.sn(iid))])
            disp('*************************************************************')
        end
        infile = fullfile(plotpar.procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%0.4d%s.use',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
        if ~exist(infile,'file')
            infile = fullfile(plotpar.procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%3.3d%s.use',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
        end
        if sum(strcmp({'MC' 'ODOMC'},id_z_sn.inst{iid})) && plotpar.proclvl==3
            infile1 = fullfile(plotpar.procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%0.4d%s.microcat',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
            if exist(infile1,'file')
                infile = infile1;
            end
        end            
        iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));

        %read data into structure array
        fileopen=fopen(infile,'r');
        if fileopen>0
            varstr = ['yy:mm:dd:hh:' id_z_sn.vars{iid}];
            d = rodbload(infile,varstr);
            clear data
            data.jd=julian(d{1},d{2},d{3},d{4});
            vars = split(id_z_sn.vars{iid},':');
            for no = 1:length(vars)
                data.(vars{no}) = d{no+4};
                data.(vars{no})(data.(vars{no})==-9999) = NaN;
            end
            if ~plotpar.unfilt
                sampling_rate = 1/median(diff(data.jd));
                for no = 1:length(vars)
                    ii = find(~isnan(data.(vars{no})));
                    data.(vars{no})(ii) = auto_filt(data.(vars{no})(ii),sampling_rate,1/2,'low',4);
                end
            end
            switch datatype
                case 'currents'
                    data.spd = sqrt(data.u.^2+data.v.^2);
                    data.dir = atan2(data.v,data.u)*180/pi;
                    data.dir(data.dir<0) = data.dir(data.dir<0)+360;
                    %current_speed for setting y-axis after plotting
                    plots.speed.ylim(2)=max(max(data.spd),plots.speed.ylim(2));
                    if isfield(data,'w')
                        data.www = data.w;
                        plots.www.ylim(2) = max(max(data.www),plots.www.ylim(2));
                    end
                case 'properties'
            end

            alldata.(iname) = data;
            id_z_sn.data_loaded(iid) = true;
        else
            disp('File does not exist!')
            disp(['infile = ' infile])
            alldata.(iname) = [];
        end

    end
end


% ------------------------------
% Plotting section
% ------------------------------
function plots = make_plots(plots, id_z_sn, alldata, plotpar)

plot_types = fieldnames(plots);

if plotpar.plot_interval==0
    plotpar.plot_interval = [plotpar.times [1 0; 1 0]];
    plotpar.plot_interval(:,2) = plotpar.plot_interval(:,2)+[-1;1];
    if plotpar.plot_interval(1,2)==0
        plotpar.plot_interval(1,2)=12; plotpar.plot_interval(1,1)=plotpar.plot_interval(1,1)-1;
    end
    if plotpar.plot_interval(2,2)==13
        plotpar.plot_interval(2,2)=1; plotpar.plot_interval(2,1)=plotpar.plot_interval(2,1)+1;
    end
end

jd1 = julian(plotpar.plot_interval(1,:));
jd2 = julian(plotpar.plot_interval(2,:)); 

% create xtick spacings based on start of months     
check=0;
i=2;
xticks(1,:)=plotpar.plot_interval(1,:);
while check~=1
    xticks(i,:)=xticks(i-1,:);
    if xticks(i,2)<12
        xticks(i,2)=xticks(i-1,2)+1;
    else
        xticks(i,2)=1;
        xticks(i,1)=xticks(i-1,1)+1;
    end
    if xticks(i,:)>=plotpar.plot_interval(2,:)
        check = 1;
    end
    i=i+1;
end

if i<4
   jdxticks =julian(plotpar.plot_interval(1,:)):(julian(plotpar.plot_interval(2,:))-julian(plotpar.plot_interval(1,:)))/5:julian(plotpar.plot_interval(2,:));
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

%set print area of figure
pos1  = [1/8*plotpar.scnsize(3),8*plotpar.bdwidth,1/2*plotpar.scnsize(3),(plotpar.scnsize(4) - 30*plotpar.bdwidth)];
for no = 1:length(plot_types)
    plots.(plot_types{no}).h = figure('Position',pos1);
    set(plots.(plot_types{no}).h,'PaperUnits','centimeters','PaperType', 'A4','PaperOrientation',plotpar.layout);
    papersize = get(plots.(plot_types{no}).h,'PaperSize');
    left = (papersize(1)- plotpar.width)/2; bottom = (papersize(2)- plotpar.height)/2;
    figuresize = [left, bottom, plotpar.width, plotpar.height];
    set(plots.(plot_types{no}).h, 'PaperPosition', figuresize);
end

%actually plot
for iid = 1:length(id_z_sn.id)
    iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
    if strcmp(plots.(plot_types{no}).var,'o2') && ~strncmp('ODO',iname,3)
        continue
    end
    b = length(alldata.(iname).jd); a = plotpar.subsamp;
    for no = 1:length(plot_types)
        if isfield(alldata.(iname),plots.(plot_types{no}).var)
            figure(plots.(plot_types{no}).h); hold on
            if iid<=length(plotpar.colours)/2
                ls = '-';
            else
                ls = '--';
            end
            plot(alldata.(iname).jd(1:a:b)-jd1,alldata.(iname).(plots.(plot_types{no}).var)(1:a:b),ls,'color',plotpar.colours(iid))
        end
    end
end
legend_text = [num2str(id_z_sn.z(:)) repmat(' m',length(id_z_sn.z),1)];
for no = 1:length(plot_types)
    figure(plots.(plot_types{no}).h)
    ylabel(plots.(plot_types{no}).ylabel)
    xlim([0 jd2-jd1]);
    if isfield(plots.(plot_types{no}),'ylim')
        plots.(plot_types{no}).ylim = ceil(plots.(plot_types{no}).ylim/10)*10;
        ylim(plots.(plot_types{no}).ylim)
    end
    set(gca,'YMinorTick','on','xTickLabel',xticklabels,'XTick',jdxticks-jd1);
    set(gca,'fontsize',14);
    if isfield(plots.(plot_types{no}),'ytick')
        set(gca,'ytick',plots.(plot_types{no}).ytick)
    end
    title(plots.(plot_types{no}).title)
    legend(legend_text)
    text(1,-20,num2str(xticks(1,1)),'FontSize',12);
    if plotpar.plot_x_labels>0
        for i=1:length(year_indexes) %***
            text((jd2-jd1)*(year_indexes(i)-1)/(length(xticklabels)-1),-20,num2str(xticks(year_indexes(i),1)),'FontSize',10);
        end
    end
    print('-dpng',plots.(plot_types{no}).outfile)
    savefig(plots.(plot_types{no}).outfile)
end
