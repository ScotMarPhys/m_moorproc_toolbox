% function t_s_diagram(moor,'procpath','layout','time_period')
%
% Function for plotting T-S diagram of a whole mooring
%
% required inputs:-
%   moor: complete mooring name as string. e.g. 'wb1_1_200420'
%
% optional inputs:-
%   layout: orientation of figure portrait/lanscape (default = portrait)
%           input of 'landscape' or 'portrait'
%           e.g. t_s_diagram('wb1_1_200420','layout','landscape')
%   procpath: can specify exact procpath if not using standard data paths. 
%           e.g. t_s_diagram('wb1_1_200420','inpath','/Volumes/noc/mpoc/hydro/rpdmoc/rapid/data/moor/proc/')
%   time_period: if want to isolate a certain section of the timeseries
%           e.g. t_s_diagram('wb1_1_200420','time_period',[2004,08,23,12,00,00],[2004,12,17,12,00,00])
%                                                          y    m  d  H  M  S     y    m  d  H  M  S
%
% functions called:-
%   rodbload, gsw_SA_from_SP, gsw_pt_from_t, gsw_SP_from_C
%   from .../exec/working/tools and .../exec/working/rodb and .../exec/working/mfiles/gsw paths
% 
% Routine written by Darren Rayner July 2021.

function t_s_diagram(moor,varargin)

global MOORPROC_G

%defaults
layout = 'portrait';
procpath = fullfile(MOORPROC_G.moordatadir,'proc');
outpath = fullfile(MOORPROC_G.reportdir,'figs');
time_period = NaN; %use whole series
%and optional inputs overwrite them
n = 2; 
while n<=nargin
    if ischar(varargin{n}) 
        if strcmp(varargin{n},'time_period')
            time_period = [datenum(varargin{n+1}) datenum(varargin{n+2})];
            n = n+3;
        else
            eval([varargin{n} ' = varargin{n+1};'])
            n = n+2;
        end
    end
end

infofile = fullfile(procpath,moor,[moor 'info.dat']);

% Load vectors of mooring information
% id instrument id, sn serial number, z nominal depth of each instrument
% s_t, e_t, s_d, e_d start and end times and dates
% lat lon mooring position, wd corrected water depth (m)
% mr mooring name
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,...
    'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');


disp('z : instrument id : serial number')
for i = 1:length(id)
    disp([z(i),id(i),sn(i)])
end

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
    if ~isempty(id_z_sn.dirs{iid}) && ((contains(id_z_sn.vars{iid},'c') && contains(id_z_sn.vars{iid},'p')) || contains(id_z_sn.vars{iid},'s'))
        disp('*************************************************************')
        disp(['Reading ' id_z_sn.inst{iid} ' - ',num2str(id_z_sn.sn(iid))])
        disp('*************************************************************')
        iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
        infile = fullfile(procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%0.4d%s.use',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
        if strcmp(id_z_sn.inst{iid},'SG')
            infile = fullfile(procpath,moor,id_z_sn.dirs{iid},sprintf('%s_%3.3d%s.use',moor,id_z_sn.sn(iid),id_z_sn.suf{iid}));
        end
        iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
        %read data into structure array
        fileopen=fopen(infile,'r');
        if fileopen>0
            if contains(id_z_sn.vars{iid},'c')
                [yy,mm,dd,hh,data.t,data.c,data.p] = rodbload(infile,'yy:mm:dd:hh:t:c:p');
                data.c(data.c<-999) = NaN;
                data.s = gsw_SP_from_C(data.c,data.t,data.p);
            else
                [yy,mm,dd,hh,data.t,data.s] = rodbload(infile,'yy:mm:dd:hh:t:s');
            end
            data.jd=julian(yy,mm,dd,hh);
        else
            disp('File does not exist!')
            disp(['infile = ' infile])
            data = [];
        end

        if ~isempty(data)
            data.p(data.p<0) = NaN;
            data.s(data.s<30) = NaN;
            data.t(data.t<-2) = NaN;
            data.absS = gsw_SA_from_SP(data.s,data.p,lon,lat);
            data.potemp = gsw_pt_from_t(data.absS,data.t,data.p);
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
colours = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];

for iid = 1:length(id_z_sn.id)
    figure(plot1); hold on
    iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
    iic = iid-floor(iid/7)*7+1;
    plot(alldata.(iname).s,alldata.(iname).potemp,'.','color',colours(iic,:));
end

figure(plot1);
xlabel('salinity');
ylabel('potemp (deg C)');
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
s = regexprep(moor,'_','\\_');
if ~isnan(time_period)
    title(['T-S diagram for mooring: ' s ' for period ' datestr(time_period(1),'yyyy-mmm-dd HH:MM:SS') ' to ' datestr(time_period(2),'yyyy-mmm-dd HH:MM:SS')])
else
    title(['T-S diagram for mooring: ' s])
end


X_limits=get(gca,'xlim');
% label data with nominal instrument depth
label_x_positions=(X_limits(2)-X_limits(1))*1.005+X_limits(1);

for iid = 1:length(id_z_sn.id)
    iname = sprintf('%s_%d', id_z_sn.inst{iid}, id_z_sn.sn(iid));
    label_y_positions = nanmedian(alldata.(iname).potemp);
    iic = iid-floor(iid/7)*7+1;
    text(label_x_positions,label_y_positions,sprintf('%d m',id_z_sn.z(iid)),'FontSize',8,'color',colours(iic,:));
end

%keyboard
