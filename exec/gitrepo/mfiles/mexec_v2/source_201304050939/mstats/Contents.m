% mstats : mexec substitute routines for functions in the Matlab stats toolbox
%
% Files
%   m_nancumsum  - function y = m_nancumsum(x,[dim]); % version of cumsum to handle NaNs
%   m_nanmean    - function y = m_nanmean(x,[dim]); % version of nanmean to avoid stats toolbox
%   m_nanmedian  - function y = m_nanmedian(x,[dim]); % version of nanmedian to avoid stats toolbox
%   m_nanprctile - function y = m_nanprctile(x,pc,[dim]); % version of prctile to avoid stats toolbox
%   m_nanstd     - function y = m_nanstd(x,[dim]); % version of nanstd to avoid stats toolbox
%   m_nansum     - function y = m_nansum(x,[dim]); % version of nansum to avoid stats toolbox
%   m_nanvar     - function y = m_nanvar(x,[dim]); % version of nanvar to avoid stats toolbox
