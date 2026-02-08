function batch_data_report(moor,varargin)

% optional inputs = 'unfiltered' or 'procpath','.....' where the entry after
% procpath gives the path to the proc directory if not using standard paths
% relative to data_report_tools_directory

%parse defaults and optional inputs
if nargin>1
    plotpar = plot_options(moor, varargin{1});
else
    plotpar = plot_options(moor);
end

overlay_plots(moor, 'currents', plotpar)
overlay_plots(moor, 'properties', plotpar)

%below don't actually use plotpar (yet?)***
plot_stacked(moor, plotpar)
currents_stacked(moor, plotpar)
%fix these (use moor_inoutpaths, use loading code from currents_stacked,
%etc.) later***
%stick_plot(moor, plotpar)
%progressive_vector(moor, plotpar)

