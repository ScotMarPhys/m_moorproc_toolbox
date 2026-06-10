function plotpar = plot_options(moor, varargin)
% sets parameters to be used by various plotting scripts
% to return defaults:
%   plotpar = plot_options(moor)
% 
% defaults can be overwritten by passing them in varargin as structure,
% e.g.
%   p.subsamp = 2; p.non_verbose = 1; plotpar = plot_options(moor, p)
% or as parameter-value pairs, e.g.
%   plotpar = plot_options(moor, 'unfiltered', 1, 'plot_x_labels', 1);
%

global MOORPROC_G

%defaults
pd = moor_inoutpaths('reports', moor);
plotpar.procpath = fullfile(MOORPROC_G.moordatadir,'proc');
plotpar.infofile = fullfile(plotpar.procpath,moor,[moor 'info.dat']);
plotpar.outpath = fullfile(pd.figsdir,moor);
plotpar.unfilt = 0;
plotpar.non_verbose = 0;
plotpar.proclvl = 2;
plotpar.plot_interval = 0;
plotpar.layout = 'portrait';

plotpar.subsamp = 1; %set higher than 1 to subsample data to ease visualisation

plotpar.width=26; 
plotpar.height=17;
plotpar.num_to_plot = 2;
plotpar.bdwidth = 5;
plotpar.topbdwidth = 30;
set(0,'Units','pixels') 
plotpar.scnsize = get(0,'ScreenSize');
plotpar.plot_x_labels = 0;
plotpar.colours = 'kbrgcmykbrgcmy'; %or just kbkbkb for properties?

%and optional inputs overwrite them
if nargin>1
    if iscell(varargin{1})
        varargin = varargin{1};
    end
    n = 1;
    while n<=length(varargin)
        if isstruct(varargin{n})
            fn = fieldnames(varargin{n});
            for no = 1:length(fn)
                plotpar.(fn{no}) = varargin{n}.(fn{no});
            end
            n = n+1;
        elseif strcmp(varargin{n},'unfiltered')
            plotpar.unfilt = 1;
            n = n+1;
        elseif strcmp(varargin{n},'non-verbose')
            plotpar.non_verbose = 1;
            n = n+1;
        else
            plotpar.(varargin{n}) = varargin{n+1};
            n = n+2;
        end
    end
end

if ischar(plotpar.proclvl); plotpar.proclvl = str2double(plotpar.proclvl); end
if plotpar.unfilt
    plotpar.proclvlstr = [num2str(plotpar.proclvl) '_unfilt'];
    plotpar.proclvltext = 'unfiltered';
else
    plotpar.proclvlstr = [num2str(plotpar.proclvl) '_lpfilt'];
    plotpar.proclvltext = 'low pass filtered';
end
