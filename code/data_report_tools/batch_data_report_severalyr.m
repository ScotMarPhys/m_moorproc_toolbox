function batch_data_report_severalyr(moorlist,varargin)

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


if a>0 
% %     pressure_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
% %     conductivity_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
% %     currents_overlay(moor,'procpath',procpath,'proclvl',proclvl,'layout','landscape','unfiltered')
% %    currents_stacked(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
% %    stick_plot(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
% %    progressive_vector(moor,'procpath',procpath,'proclvl',proclvl)
% %     temperature_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
% %     salinity_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')
% %     pden1000_overlay(moor,'procpath',procpath,'proclvl',proclvl,'unfiltered')      
%%%    stats_table(moor,'procpath',procpath)
else
     overlay_plot_severalyr(moorlist,'procpath',procpath,'proclvl',proclvl) % for T-S data
     currents_stacked_severalyr(moorlist,'procpath',procpath,'proclvl',proclvl)

end

