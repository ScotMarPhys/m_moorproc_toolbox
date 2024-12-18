function batch_data_report(moor,varargin)

% optional inputs = 'unfiltered' or 'procpath','.....' where the entry after
% procpath gives the path to the proc directory if not using standard paths
% relative to data_report_tools_directory

%parse defaults and optional inputs
inargs = varargin;
plot_options

overlay_plots(moor, 'currents', plotpar)
overlay_plots(moor, 'properties', plotpar)

plot_stacked(moor, plotpar)
currents_stacked(moor, plotpar)
stick_plot(moor, plotpar)
progressive_vector(moor,plotpar)

