function batch_data_report(moor,varargin)

% optional inputs = 'unfiltered' or 'procpath','.....' where the entry after
% procpath gives the path to the proc directory if not using standard paths
% relative to data_report_tools_directory

% check for optional arguments
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
    proclvl=char(varargin(a+1));
else
    proclvl='2';
end


a=strmatch('unfiltered',varargin,'exact');

b=strmatch('plot_interval',varargin,'exact');
if b>0
    intervalstr=varargin{b+1};
else
    intervalstr=0;
end

if a>0 & isempty(b)
   plot_stacked(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
    pressure_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
    conductivity_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
    currents_overlay(moor,'procpath',procpath,'proclvl',proclvl,'layout','landscape','unfiltered')
   currents_stacked(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
  stick_plot(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
   progressive_vector(moor,'procpath',procpath,'proclvl',proclvl)
    temperature_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
    salinity_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
    pden1000_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')      
%%    stats_table(moor,'procpath',procpath)
elseif a>0 & ~isempty(b)
     pressure_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered','plot_interval',intervalstr)
     conductivity_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered','plot_interval',intervalstr)
     currents_overlay(moor,'procpath',procpath,'proclvl',proclvl,'layout','landscape','unfiltered','plot_interval',intervalstr)
    currents_stacked(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered','plot_interval',intervalstr)
    stick_plot(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered','plot_interval',intervalstr)
    progressive_vector(moor,'procpath',procpath,'proclvl',proclvl,'plot_interval',intervalstr)
     temperature_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered','plot_interval',intervalstr)
     salinity_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered','plot_interval',intervalstr)
     pden1000_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered','plot_interval',intervalstr)      
%%%    stats_table(moor,'procpath',procpath)
elseif isempty(a) & ~isempty(b)
     pressure_overlay(moor,'procpath',procpath,'proclvl',proclvl,'plot_interval',intervalstr)
     conductivity_overlay(moor,'procpath',procpath,'proclvl',proclvl,'plot_interval',intervalstr)
     currents_overlay(moor,'procpath',procpath,'proclvl',proclvl,'layout','landscape','plot_interval',intervalstr)
    currents_stacked(moor,'procpath',procpath,'proclvl',proclvl,'plot_interval',intervalstr)
    stick_plot(moor,'procpath',procpath,'proclvl',proclvl,'plot_interval',intervalstr)
    progressive_vector(moor,'procpath',procpath,'proclvl',proclvl,'plot_interval',intervalstr)
     temperature_overlay(moor,'procpath',procpath,'proclvl',proclvl,'plot_interval',intervalstr)
     salinity_overlay(moor,'procpath',procpath,'proclvl',proclvl,'plot_interval',intervalstr)
     pden1000_overlay(moor,'procpath',procpath,'proclvl',proclvl,'plot_interval',intervalstr)   
elseif isempty(a) & isempty(b)
     pressure_overlay(moor,'procpath',procpath,'proclvl',proclvl)
     conductivity_overlay(moor,'procpath',procpath,'proclvl',proclvl)
     currents_overlay(moor,'procpath',procpath,'proclvl',proclvl,'layout','landscape')
     currents_stacked(moor,'procpath',procpath,'proclvl',proclvl)
     stick_plot(moor,'procpath',procpath,'proclvl',proclvl)
     progressive_vector(moor,'procpath',procpath,'proclvl',proclvl)
     temperature_overlay(moor,'procpath',procpath,'proclvl',proclvl)
     salinity_overlay(moor,'procpath',procpath,'proclvl',proclvl)
    pden1000_overlay(moor,'procpath',procpath,'proclvl',proclvl)    
%%% stats_table(moor,'procpath',procpath)
% see how to set up the plot correlation :
% currentdir = pwd;
% cd vel_vs_temp_figures
% vel_vs_temp_corr
% v_and_t_vs_press_corr
% Add a correlation function U,V with W
end

