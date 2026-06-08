function [arr]=onesnan(s1,s2);
% function [arr]=onesnan(s1,s2)
%
% onesnan is like 'ones' but accepts NaN and then gives []

% Gerd Krahmann, IfM Kiel, Apr 1993

if (isnan(s1)+isnan(s2))
  arr=[];
else
  arr=ones(s1,s2);
end
