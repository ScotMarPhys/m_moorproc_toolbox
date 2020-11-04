% function x=big_fig(varargin)
% function to produce a large figure window for use with plotting in other routines
% varargin options = number (figure number)
%                  = 'portrait' (portrait layout)
%                  = 'landscape' (landscape layout, DEFAULT)
% e.g. x=big_fig(3,'portrait')
% or can call without any inputs whihc will just use the next available
% figure number and use landscape layout
% e.g. x=big_fig

function x=big_fig(varargin)

%get screen size
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 

% 2/8/12 - new line to address when have multiple screens
monitor_positions=get(0,'MonitorPositions');
%scnsize = get(0,'ScreenSize');

scnsize=monitor_positions(1,:);

if nargin ==1
    if isnumeric(varargin{1})
        fig_number=varargin{1};
        way_round='landscape';
    else
        fig_number=0;
        way_round=varargin{1};
    end
elseif nargin>1
    fig_number=varargin{1};
    way_round=varargin{2};
else
    fig_number=0;
    way_round='landscape';
end


%set print area of figure
if strncmpi(way_round,'portrait',8)
    pos = [1/8*scnsize(3),8*bdwidth,1/2*scnsize(3),(scnsize(4) - 30*bdwidth)];
else
    pos = [1/10*scnsize(3),8*bdwidth,8/10*scnsize(3),(scnsize(4) - 30*bdwidth)];
end

if fig_number==0;
    x=figure('Position',pos);
else
    x=figure(fig_number);
    set(x,'Position',pos)
end

set(x,'PaperUnits','centimeters');
set(x, 'PaperType', 'A4');
set(x, 'PaperOrientation',way_round);
papersize = get(x,'PaperSize');
if strncmpi(way_round,'portrait',8)
    width=17; height=26; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
else
    width=26; height=17; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
end
figuresize = [left, bottom, width, height];
set(x, 'PaperPosition', figuresize);

end
