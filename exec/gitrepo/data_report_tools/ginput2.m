function  [X,Y,BUTTON,SCALEMAT] = ginput2(varargin)
%GINPUT2   Graphical input from mouse with zoom, plot and georeference.
%   
%   SYNTAX:
%                       XY  = ginput2;
%                     [X,Y] = ginput2;
%              [X,Y,BUTTON] = ginput2;
%     [X,Y,BUTTON,SCALEMAT] = ginput2;
%                     [...] = ginput2(SCALEOPT);
%                     [...] = ginput2(PLOTOPT);
%                     [...] = ginput2(SCALEOPT,PLOTOPT);
%                     [...] = ginput2(...,'keepzoom');
%                     [...] = ginput2(N,...);
%
%   INPUT:
%     SCALEOPT   - Single logical specifying whether the IMAGE should be
%                  interactively scaled (georeferenced), or it can be the
%                  2x4 SCALEMAT matrix for automatically scaling.
%                  DEFAULT: false (do not scales/georeferences)
%     PLOTOPT    - String and/or parameter/value pairs specifying the drawn
%                  points properties (see PLOT for details). 
%                  DEFAULT: 'none' (do not plot any point)
%     'keepzoom' - When finishing selection by default the zoom is
%                  restored. By using this option this feature is ignored.
%                  DEFAULT: 'unzoom' (restores original axis limits)
%     N          - Number of points to be selected. One of 0,1,2,...,Inf
%                  DEFAULT: Inf (selects until ENTER or ESCAPE is pressed)
%
%   OUTPUT:
%     XY        - X- and Y-coordinate(s).
%     X         - X-coordinate(s).
%     Y         - Y-coordinate(s).
%     BUTTON    - Pressed button.
%     SCALEMAT  - 2x4 matrix specifying the coordinates of two different
%                 points (1) and (2) in the Actual image coordinates and 
%                 the Real one: [Ax1 Ay1 Ax2 Ay1; Rx1 Ry1 Rx2 Ry2] to
%                 be use for scaling/georeferencing.
%
%   DESCRIPTION:
%     This program uses MATLAB's GINPUT function to get the coordinates
%     of a mouse-selected point in the current figure (see GINPUT for
%     details), but with four major improvements:
%                  1. ZOOMING
%                  2. DELETING
%                  3. PLOTING
%                  4. SCALING or GEOREFERENCE IMAGES.
%     The differences are:
%      a) Obviously, the SCALEOPT and PLOTOPT optional arguments.
%      b) When click is made outside the axes, it is ignored.
%      c) When LEFT-click, ZOOM-IN is performed right into the selected
%         point (PANNING).
%      d) When RIGHT-click, the point is selected (normal).
%      e) When DOUBLE-click, ZOOM-OUT is done.
%      f) When MIDDLE-click, ZOOM-RESET is done (see ZOOM for details).
%      g) When pressed any KEY follows the next rules: 
%          A) If ENTER is pressed, the selection is terminated. If no point
%             was already selected, the outputs are empty's.
%          B) If BACKSPACE key is pressed, the last selected point is
%             deleted and the selection continues.
%          C) If SPACEBAR the mouse current position or NANs coordinates
%             are saved, depending whether the mouse was inside or outside
%             any of the current figure axes, respectively. In this latter
%             case, the selection is NOT counted as one of the N points.
%             Besides, when drawing the color is changed. Then, the outputs
%             may not be of length N.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * The function can be used for interactively digitalize/vectorize
%       RASTER images.
%     * The function can be used only as a georeference function with 
%       >> ginput2(0,true)
%     * The scale/georeference only works when the current axes has an
%       IMAGE type children (see Image for details). 
%     * The x and y data from axes and image are changed when scale/
%       georeference is used.
%     * The drawn points are deleted from the graphics once the selection
%       is finished. 
%     * The priority of the inputs are: N, then SCALEOPT and finally
%       PLOTOPT. If the first (integer) is missing, the next is taken into
%       account (logical or 2x4 matrix) and so on.
%
%   EXAMPLE:
%     % Selects until ENTER pressed.
%       xy = ginput2;
%     % Selects 5 points.
%       [x,y] = ginput2(5);
%     % Gets pressed button.
%       [x,y,button] = ginput2(1);
%
%     imagesc(peaks(40))
%     [x,y,button,scalemat] = ginput2(4,true,'k*','keepzoom');
%     hold on, plot(x,y), hold off
%
%   SEE ALSO:
%     GINPUT, PLOT.
%
%
%   ---
%   MFILE:   ginput2.m
%   VERSION: 2.1 (Aug 20, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jul 09, 2008)
%   2.0      Changed default YESERROR value and fixed a bug there. Changed
%            behavior when N==1. Fixed bug with zoom out. Changed default
%            selection keys. Changed default selection click mouse: from
%            left one to the right one (Jun 08, 2009).
%   2.1      Fixed bugs related with points delection. Added new 'keepzoom'
%            feature (Aug 20, 2009). 

%   DISCLAIMER:
%   ginput.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2008-2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% PARAMETERS
% Defaults:
X        = [];
Y        = [];
BUTTON   = [];
SCALEMAT = [];
N         = Inf;
SCALEOPT  = false;
PLOTOPT   = {'none'};
UNZOOM    = 'unzoom';
% Constants KEYs (on my personal keyboard):
DOUBLECLICK    =   0;
LEFTCLICK      =   1;
MIDDLECLICK    =   2;
RIGHTCLICK     =   3;
BACKSPACE      =   8;
ESCAPE         =  27;
LEFTARROW      =  28;
RIGHTARROW     =  29;
UPARROW        =  30;
DOWNARROW      =  31;
SPACEBAR       =  32;
DELETE         = 127;
ASCII          = [ ...
                    33:64  ...  UP-KEYS
                    65:90  ...  UP-LETTERS
                    91:96  ... LOW-KEYS
                    97:122 ... LOW-LETTERS
                   123:126 ... LOW-KEY
                   161:255 ...     FOREING
                   ];
% Functionality:
% NOTE: I left all this KEYs because the user may use this special case for
% other special purposes outside this function.
% % First version default:
% % SELECTS   = [LEFTCLICK ASCII ESCAPE LEFTARROW RIGHTARROW ...
% %              UPARROW DOWNARROW SPACEBAR DELETE];
% % ZOOMIN    = RIGHTCLICK; 
SELECTS   = [RIGHTCLICK SPACEBAR]; % Selection  buttons
DELETES   = BACKSPACE;             % Deletion   buttons
FINISHES  = [];                    % Finishes   buttons
ZOOMIN    = LEFTCLICK;             % ZOOM(2)    buttons
ZOOMRESET = MIDDLECLICK;           % ZOOM RESET buttons
ZOOMOUT   = DOUBLECLICK;           % ZOOM OUT   buttons
% Other parameters
secpause  = 0.3;    % Seconds to wait for double-click response.
YESERROR  = false;  % If there is an error with GINPUT, it tells to display 
                    % an ERROR or a WARNING message.
                    
% Checks number of inputs:
if nargin>4
 error('CVARGAS:ginput2:tooManyInputs',...
  'At most 4 inputs are allowed.')
elseif nargout>4
 error('CVARGAS:ginput2:tooManyOutputs',...
  'At most 4 outputs are allowed.')
end

% Checks input:
Nargin = nargin;
% First checks N:
if Nargin && ~isempty(varargin{1}) && isfloat(varargin{1}) && ...
  numel(varargin{1})==1 && round(varargin{1})==varargin{1} && ...
  varargin{1}>=0
 N           = varargin{1};
 varargin(1) = [];
 Nargin      = Nargin-1;
end
% Second checks SCALEOPT:
if Nargin && ~isempty(varargin{1}) && (((length(varargin{1})==1) && ...
  islogical(varargin{1})) || all(size(varargin{1})==[2 4]))
 SCALEOPT    = varargin{1};
 varargin(1) = [];
 Nargin      = Nargin-1;
end
% Now checks PLOTOPT and UNZOOM:
if Nargin && ~isempty(varargin{1}) && ischar(varargin{1})
 if ischar(varargin{end}) && (strcmpi(varargin{end},'unzoom') || ...
    strcmpi(varargin{end},'keepzoom'))
  UNZOOM        = lower(varargin{end});
  varargin(end) = [];
  Nargin        = Nargin-1;
 end
 PLOTOPT = varargin;
end
clear varargin Nargin

% SCALES/GEOREFERENCE?:
ha       = [];
hf       = [];
if (islogical(SCALEOPT) && SCALEOPT) || (numel(SCALEOPT)==8)
 method = 'linear';
 extrap = 'extrap';
 ha     = gca;
 hi     = findobj(get(ha,'Children'),'Type','image');
 axes(ha)
 if ~isempty(hi)
  hi    = hi(1);
  xlim  = get(ha,'XLim');
  ylim  = get(ha,'YLim');
  zlim  = get(ha,'ZLim');
  z     = repmat(max(zlim),1,5);
  xdata = get(hi,'XData');
  ydata = get(hi,'YData');
  if islogical(SCALEOPT) % interactivelly
   Ax1 = round(min(xdata)); Ax2 = round(max(xdata));
   Ay1 = round(min(ydata)); Ay2 = round(max(ydata));
   % Default (equal):
   Rx1 = Ax1; Rx2 = Ax2;
   Ry1 = Ay1; Ry2 = Ay2;
   hgeo     = [];
   dlgTitle = 'Georeference image';
   lineNo   = 1;
  
   while true
    % Selects first corner:
    theans = ...
          questdlg('Select the first corner (1 of 2):',dlgTitle,'OK','OK');
    if ~strcmp(theans,'OK'), return, end
    pause(secpause)
    
    [Ax1,Ay1] = ginput2(1,false,'none','unzoom');
    Ax1       = round(Ax1);
    Ay1       = round(Ay1);
    if ~ishandle(ha), return, end
    if ha==gca && ~isempty(Ax1) && ~isnan(Ax1)
     axis(ha,[xlim ylim])
     hgeo(1) = line([xlim NaN Ax1 Ax1],[Ay1 Ay1 NaN ylim],z,'color','m');
     prompt  = {'X-coordinate at 1st corner:',...
                'Y-coordinate at 1st corner:'};
     def     = {int2str(Ax1),int2str(Ay1)};
     answer  = inputdlg(prompt,dlgTitle,lineNo,def);
     answer  = str2num(char(answer{:}));
     break
    end
   end
   axes(ha)
   
   % Checks inputs:
   if ~isempty(answer) && isfloat(answer) && length(answer)==2 && ...
                                                      all(isfinite(answer))
    Rx1 = answer(1); Ry1 = answer(2);
    secondcorner = true;
   else
    secondcorner = false;
    warning('CVARGAS:ginput2:incorrectGeoreference',...
            'Ignored incorrect georeference corners.')
   end
  
   while secondcorner
    % Selects second corner:
    theans = ...
         questdlg('Select the second corner (2 of 2):',dlgTitle,'OK','OK');
    if ~strcmp(theans,'OK'), return, end
    pause(secpause)
   
    [Ax2,Ay2] = ginput2(1,false,'none','unzoom');
    Ax2       = round(Ax2);
    Ay2       = round(Ay2);
    if ~ishandle(ha), return, end
    if ha==gca && ~isempty(Ax2) && ~isnan(Ax2) && Ax2~=Ax1 && Ay2~=Ay1
     axis(ha,[xlim ylim])
     hgeo(2) = line([xlim NaN Ax2 Ax2],[Ay2 Ay2 NaN ylim],z,'color','c');
     prompt  = {'X-coordinate at 2nd corner:',...
                'Y-coordinate at 2nd corner:'};
     def     = {int2str(Ax2),int2str(Ay2)};
     answer  = inputdlg(prompt,dlgTitle,lineNo,def);
     answer  = str2num(char(answer{:}));
     break
    end
   end
   axes(ha)
   
   % Checks inputs:
   if secondcorner && ~isempty(answer) && isfloat(answer) && ...
                         length(answer)==2 && all(isfinite(answer))
    Rx2 = answer(1); Ry2 = answer(2);
   else
    warning('CVARGAS:ginput2:incorrectGeoreference',...
            'Ignored incorrect georeference corners.')
   end
  
   % Deletes corner's lines:
   if any(ishandle(hgeo))
    delete(hgeo(ishandle(hgeo)))
   end
   
   % Scale matrix:
    SCALEMAT = [Ax1 Ay1 Ax2 Ay2; Rx1 Ry1 Rx2 Ry2];
  else
   SCALEMAT = SCALEOPT;
  end
 else
  warning('CVARGAS:ginput2:noImageFound',...
   'No image found in the current axes to georeference.')
 end
 
 % OK, set the scaling then:
 if ~isempty(SCALEMAT)
  xdata = interp1(SCALEMAT(1,[1 3]),SCALEMAT(2,[1 3]),xdata,method,extrap);
  ydata = interp1(SCALEMAT(1,[2 4]),SCALEMAT(2,[2 4]),ydata,method,extrap);
  xlim2 = interp1(SCALEMAT(1,[1 3]),SCALEMAT(2,[1 3]),xlim ,method,extrap);
  ylim2 = interp1(SCALEMAT(1,[2 4]),SCALEMAT(2,[2 4]),ylim ,method,extrap);
  set(hi,'XData',xdata);
  set(hi,'YData',ydata);  
  set(ha,'XLim' ,sort(xlim2,'ascend'));
  set(ha,'YLim' ,sort(ylim2,'ascend'));
  % Reverses axis directions:
  if diff(xlim)*diff(xlim2)<1
   if strcmp(get(ha,'XDir'),'normal')
    set(ha,'XDir','reverse')
   else
    set(ha,'XDir','normal')
   end
  end
  if diff(ylim)*diff(ylim2)<1
   if strcmp(get(ha,'YDir'),'normal')
    set(ha,'YDir','reverse')
   else
    set(ha,'YDir','normal')
   end
  end
 end
 axis(ha,'normal')
 
end

% DRAWS?:
if strcmpi(PLOTOPT{1},'none')
 yesdraw = false;
else
 yesdraw = true;
end
% Optional parameters:
if yesdraw
 hpoints  = [];
 % Check for linestyle color:
 yescolor     = true;
 Nplotopt     = length(PLOTOPT);
 yeslinestyle = rem(Nplotopt,2);
 if yeslinestyle % Given LineStyle
  for k=1:length(PLOTOPT{1})
   switch lower(PLOTOPT{1}(k))
    case 'y', yescolor = false; break
    case 'm', yescolor = false; break
    case 'c', yescolor = false; break
    case 'r', yescolor = false; break
    case 'g', yescolor = false; break
    case 'b', yescolor = false; break
    case 'w', yescolor = false; break
    case 'k', yescolor = false; break
    otherwise, % no color specified
   end
  end
 end
 if ~yescolor && (Nplotopt*yeslinestyle~=1)
  for k = yeslinestyle+1:2:Nplotopt    % Given 'Color'
   if strncmpi(PLOTOPT{k},'co',2), yescolor = false; break, end
  end
 end
 if yescolor
  contnan  = 1;
  colors   = get(gca,'ColorOrder');
  ncolors  = size(colors,1);
  color    = colors(1,:);
 end
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

cont     = 0;
alim.ha  = [];
alim.la  = {};

while cont<N     % Principal loop
 
 % GINPUT:
 try
  [x,y,button] = ginput(1);
 catch % Changed for compatibility.
  % GINPUT error:
  if YESERROR
     error('CVARGAS:ginput2:executionError',lasterr)
  else
   warning('CVARGAS:ginput2:executionError',lasterr)
   if nargout<=1  % Fixed BUG 29 SEP, 2008
    X = [X Y]; 
   end
   return
  end
 end
 
 % Axes clicked:
 ha = gca;
 
 % Gets limits:
 if ~any(alim.ha==ha)
  alim.ha(end+1) = ha;
  alim.la{end+1} = axis;
 end
 
 % Sets zoom:
 zlim = getappdata(ha,'zoom_zoomOrigAxesLimits');
 if isempty(zlim) % Fixed BUG, SEP 2008
  zoom reset
  zlim = getappdata(ha,'zoom_zoomOrigAxesLimits');
 end
 
 % Checks if DOUBLE clicks:
 pause(secpause) % Gives time for response
 if strcmp(get(gcf,'selectiontype'),'open')
  button = DOUBLECLICK;
 end

 % Checks if ENTER or FINISHES button:
 if isempty(button) || ismember(button,FINISHES)
  % Finishes selection:
  break
 end
 
 % Checks if DELETION button:
 if ismember(button,DELETES)
  if ~isempty(X)
   inan = isnan(X(end));
   if yesdraw
    if ~inan
     % Deletes last drawn point:
     if ~isempty(hpoints) && ishandle(hpoints(end)) % Fixed bug Aug 2009
      delete(hpoints(end)), hpoints(end) = [];
     end
    elseif yescolor
     % Change color as necesary:
     contnan = contnan-1;
     color   = colors(mod(contnan-1,ncolors)+1,:);
    end
   end
   % Deletes the last selected point:
   X(end)      = [];
   Y(end)      = [];
   BUTTON(end) = [];
   % Checks if the last point was NaN:
   if ~inan
    cont = cont-1;
   end
  end
  continue
 end
 
 % Checks if ZOOM OUT button:
 if ismember(button,ZOOMOUT)
  zoom out
  continue
 end
 
 % Checks if ZOOM RESET button:
 if ismember(button,ZOOMRESET)
  zoom reset
  continue
 end
 
 
 % Checks if the mouse was inside an axes of the current figure;
 lim     = axis;
 outside = x<lim(1) || x>lim(2) || y<lim(3) || y>lim(4);
 
 % Checks if ZOOM IN with PAN:
 if ismember(button,ZOOMIN) && ~outside
  % First centers (do the PANNING):
  lim = [x+diff(lim(1:2))/2*[-1 1] y+diff(lim(3:4))/2*[-1 1]];
  axis(lim)
  % Then do the ZOOMING:
  zoom(2)
  continue
 end
 
 % Checks if SELECTS button:
 if ismember(button,SELECTS)
  
  % Sets NaNs if outside the axes:
  if outside 
   if ~isnumeric(N) 
    % Change color:
    if yesdraw && yescolor 
     contnan = contnan+1;
     color   = colors(mod(contnan-1,ncolors)+1,:);
    end
    % Adds NaNs but the point counters do not take it into account:
    X      = [X;      NaN];
    Y      = [Y;      NaN];
    BUTTON = [BUTTON; button];
   else
    % Ignores the point
   end
  else
   % Draws the result:
   if yesdraw
    % Search for last point:
    x0 = []; y0 = []; z0 = [];
    inan = isnan([NaN; X; NaN]);
    if ~inan(end-1)
     inan        = find(inan);
     nlastpoints = inan(end)-inan(end-1)-1;
     npoints     = length(hpoints);
     range       = npoints-nlastpoints+1:npoints;
     hlastaxes   = get(hpoints(range),'Parent');
     if iscell(hlastaxes), hlastaxes  = cell2mat(hlastaxes); end
     [loc,loc] = ismember(ha,hlastaxes);
     if loc
      x0 = get(hpoints(range(loc)),'XData');
      y0 = get(hpoints(range(loc)),'YData');
      z0 = get(hpoints(range(loc)),'ZData');
     end
    end
    holdon = ishold;
    if ~holdon, hold on, end 
     h = plot([x0 x],[y0 y],PLOTOPT{:});
     % Elevates the value:
     z = get(ha,'ZLim'); z = z(2);
     set(h,'Zdata',[z0 z])
     % Sets the color:
     if yescolor
      set(h,'Color',color)
     end
     hpoints = [hpoints; h];
    if ~holdon, hold off, end 
    % Restores limits:
    axis(lim)
   end
   
   % Centers the selected point if ZOOM-IN: 29 SEP,2008
   if all((lim~=zlim))
    lim = [x+diff(lim(1:2))/2*[-1 1] y+diff(lim(3:4))/2*[-1 1]];
    axis(lim)
   end
   
   % Saves the result:
   X      = [X;      x];
   Y      = [Y;      y];
   BUTTON = [BUTTON; button];
   cont   = cont+1;
  end
  continue
 end

 % Checks if any other button pressed inside the axes:
 if ~outside
  if 1
   X      = [X;      x];
   Y      = [Y;      y];
   BUTTON = [BUTTON; button];
   cont   = cont+1;
  else
   break
  end
 else
  % ignores the selection
 end
 
end

% Deletes drawn points if still exist:
if yesdraw && any(ishandle(hpoints))
 delete(hpoints(ishandle(hpoints)))
end

% Returns original limits.
if ~strcmp(UNZOOM,'keepzoom') && ~isempty(alim.ha)
 alim.ha(~ishandle(alim.ha)) = [];
 for k = 1:length(alim.ha)
  temp = axis(alim.ha(k));
  if ~all(temp(1:4)==alim.la{k}(1:4))
   axis(alim.ha(k),alim.la{k})
  end
 end
end


% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

if nargout<=1 
 X = [X Y]; 
end


% [EOF]   ginput2.m