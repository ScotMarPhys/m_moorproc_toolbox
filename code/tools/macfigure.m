function macfigure(vargin)
% gdm, rb1201
% modified figure to:
%   - set better default lines
%   - enlarge to fullscreen (macbook pro 13")

if nargin == 0
    figure;
elseif nargin == 1
    nfig = vargin;
    figure(nfig);
end;
% builtin('figure')

% routine by BAK to set better default line width and font size

set(gcf,'defaultaxeslinewidth',2)
set(gcf,'defaultlinelinewidth',2)
set(gcf,'defaultaxesfontsize',16)
set(gcf,'defaulttextfontsize',16)
set(gcf,'position',[1 108 893 702])