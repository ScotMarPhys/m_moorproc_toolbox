% set values of integration domain for a given step size between two pressures
% or depths or whatever
% function inc = spacing(p1,p2,step)
% input:   p1   = lower integration boundary
%          p2   = uppet integration boundary
%          step = step size of integration between p1 and p2 (positive value)
%
% T.Kanzow

function inc = spacing(p1,p2,step) 

if p1 > p2
  step = -step;
end

fac = (p2 - p1) / step;

inc = [p1 : step : [p1+floor(fac)*step] p2];
%keyboard 
