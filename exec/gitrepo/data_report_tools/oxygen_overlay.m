% function oxygen_overlay(moor,'procpath','proclvl','layout','plot_interval','unfiltered')
%
% Function for plotting oxygens of a mooring overlayed on the same axes
%
% required inputs:-
%   moor: complete mooring name as string. e.g. 'wb1_1_200420'
%
% optional inputs:-
%   layout: orientation of figure portrait/lanscape (default = portrait)
%           input of 'landscape' or 'portrait'
%           e.g. oxygen_overlay('wb1_1_200420','layout','landscape')
%   plot_interval: matrix of start and end dates for plot
%           e.g. oxygen_overlay('wb1_1_200420','plot_interval',[2004 02 01 00; 2005 06 01 00])
%           dates are:- yyyy mm dd hh. Default is calculated automatically
%   procpath: can specify exact procpath if not using standard data paths. 
%           e.g. oxygen_overlay('wb1_1_200420','inpath','/Volumes/noc/mpoc/hydro/rpdmoc/rapid/data/moor/proc/')
%   proclvl: can specify level of processing of the data to plot. 
%           e.g. 'proclvl','2': will plot the .use file ; 'proclvl','3' will plot the .microcat and .edt files
%   unfiltered: plot data in unfiltered format - input 'unfiltered'
%
% functions called:-
%   rodbload, julian, auto_filt
%   from .../exec/moor/tools and .../exec/moor/rodb paths
% 
% Routine written by Darren Rayner July 2006.
% adapted from oxygen_overlay.m

function oxygen_overlay(moor,varargin)

global MOORPROC_G

if nargin <1
    help oxygen_overlay
    return
end

%defaults
layout = 'portrait';
procpath = fullfile(MOORPROC_G.moordatadir,'proc');
outpath = fullfile(MOORPROC_G.reportdir,'figs');
non_verbose = 0;
proclvl = 2;
plot_interval = 0;
unfilt = 0;
%and optional inputs overwrite them
n = 1;
while n<=nargin-1
    if ischar(varargin{n}) 
        if strcmp(varargin{n},'non-verbose')
            non_verbose = 1;
            n = n+1;
        elseif strcmp(varargin{n},'unfiltered')
unfilt = 1;
n = n+1;
        else
            eval([varargin{n} ' = varargin{n+1};'])
            n = n+2;
        end
    end
end

infofile = fullfile(procpath,moor,[moor 'info.dat']);

proclvlstr0 = num2str(proclvl);
if unfilt
    proclvlstr = [proclvlstr0 '_unfilt'];
else
    proclvlstr = [proclvlstr0 '_lpfilt'];
end


% Load vectors of mooring information
% id instrument id, sn serial number, z nominal depth of each instrument
% s_t, e_t, s_d, e_d start and end times and dates
% not used: lat lon mooring position, wd corrected water depth (m)
% mr mooring name
[id,sn,z,s_t,s_d,e_t,e_d]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date');

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
% use year_indexes later for plotting on bottom graph

jd1 = julian(plot_interval(1,:));
jd2 = julian(plot_interval(2,:)); 

num_to_plot=2; % num_to_plot is the number of samples per day to plot and can be adjusted accordingly

id_z_sn = all_inst_table(id, z, sn);
id_z_sn.data_loaded = false(length(id),1);

%set figure size on screen for better viewing
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

%set print area of figure
pos1  = [1/8*scnsize(3),8*bdwidth,1/2*scnsize(3),(scnsize(4) - 30*bdwidth)];
plot1=figure('Position',pos1);
set(plot1,'PaperUnits','centimeters');
set(plot1, 'PaperType', 'A4');
set(plot1, 'PaperOrientation',layout);
papersize = get(plot1,'PaperSize');
width=17; height=26; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
figuresize = [left, bottom, width, height];
set(plot1, 'PaperPosition', figuresize);

plot_string={};

% -----------------------------------
% START OF READING IN INSTRUMENT DATA
% -----------------------------------
for iid=1:length(id_z_sn.id)
    if ~isempty(id_z_sn.dirs{iid}) && contains(id_z_sn.vars{iid},'o')
        if ~non_verbose
            disp('*************************************************************')
            disp(['Reading ' id_z_sn.inst{iid} ' - ',num2str(id_z_sn.sn(iid))])
            disp('*************************************************************')
        end
        iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
        infile = fullfile(procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%0.4d%s.use',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
        if strcmp(id_z_sn.inst{iid},'SG')
            infile = fullfile(procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%3.3d%s.use',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
        end
        if sum(strcmp({'MC' 'ODOMC'},id_z_sn.inst{iid})) && proclvl==3
            infile1 = fullfile(procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%0.4d%s.microcat',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
            if exist(infile1,'file')
                infile = infile1;
            end
        end
        iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
        %read data into structure array
        fileopen=fopen(infile,'r');
        if fileopen>0
            [yy,mm,dd,hh,data.o] = rodbload(infile,'yy:mm:dd:hh:o2');
            data.jd=julian(yy,mm,dd,hh);
        else
            disp('File does not exist!')
            disp(['infile = ' infile])
            data = [];
        end

        if ~isempty(data)
            data.o(data.o<-999) = NaN;
            if unfilt==0
                sampling_rate = 1/median(diff(data.jd));
                ii = find(~isnan(data.o));
                data.o(ii) = auto_filt(data.o(ii),sampling_rate,1/2,'low',4);
            end
            alldata.(iname) = data;
            id_z_sn.data_loaded(iid) = true;
        end
    end

end
%only keep the rows where we loaded data
id_z_sn = id_z_sn(id_z_sn.data_loaded,:);


% ------------------------------
% Plotting section
% ------------------------------

for iid = 1:length(id_z_sn.id)
    figure(plot1); hold on
    iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
    col = 'black'; if isodd(iid); col = 'blue'; end
    plot(alldata.(iname).jd-jd1,alldata.(iname).o,'color',col);
end

figure(plot1);
ylabel('oxygen (degC)');
xlim([0 jd2-jd1]);
set(gca,'YMinorTick','on');
set(gca,'xTickLabel',xticklabels);
set(gca,'XTick',jdxticks-jd1);
s = regexprep(moor,'_','\\_');
if unfilt==0
    title(['Low-pass filtered oxygen from mooring: ' s])
else
    title(['Unfiltered oxygen from mooring: ' s])
end

% Display year labels on bottom graph
Y_limits=ylim;   X_limits=xlim;
a=(Y_limits(1)-Y_limits(2))*1.05+Y_limits(2);
text(X_limits(1),a,num2str(xticks(1,1)),'FontSize',10);
for i=1:length(year_indexes)
    text((jd2-jd1)*(year_indexes(i)-1)/(length(xticklabels)-1),a,num2str(xticks(year_indexes(i),1)),'FontSize',10);
end


% label data with nominal instrument depth
label_x_positions=(X_limits(2)-X_limits(1))*1.005+X_limits(1);
for iid = 1:length(id_z_sn.id)
    iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
    label_y_positions = nanmedian(alldata.(iname).o);
    col = 'black'; if isodd(iid); col = 'blue'; end
    text(label_x_positions,label_y_positions,sprintf('%d m',id_z_sn.z(iid)),'FontSize',8,'color',col);
end

%keyboard
if ~exist(outpath,'dir')
    mkdir(outpath)
end
outfile = fullfile(outpath,[moor '_oxygen_overlay_proclvl_' proclvlstr]);
print('-dpng',outfile)
savefig(outfile)
