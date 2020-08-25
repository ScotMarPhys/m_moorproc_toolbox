
function [n,nidx]=nearest(x0,x)

% NEAREST look up the nearest value in a series. 
% 
%         [n,idx] = nearest(x0,x)
%         
%         inputs  x0 : search value
%                 x  : series to look up
% 
%         outputs n   : nearest value
%                 idx : index of the nearest value 

% by Ulf Garternicht, Dec 95, IfM Kiel.
% accetpting nan's d.kieke,13.03.1997

nidx=find(abs(x-x0)==minnan(abs(x-x0)));
n=x(nidx);