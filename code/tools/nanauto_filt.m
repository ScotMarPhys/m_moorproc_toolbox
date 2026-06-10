function yf = nanauto_filt(y,sr,co, jd)

% digital filter routine using a butterworth filter design  
%
%  function yf = nanauto_filt(y,sr,co,[type],[fo])
%    input :  
%       y    :  data array
%       sr   :  sampling rate (frequency)
%       co   :  cutoff rate   (frequency), 2-element vector for type =
%       'stop'
%       jd   :  time base

% Paul Wright Feb 2010 - added the nan bit!  (based on auto_filt.m)



% remove the NaNs
ix = find(~isnan(y));

yf1 = auto_filt(y(ix), sr, co);

yf = interp1(jd(ix), yf1, jd);


