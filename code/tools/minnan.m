function [y,i] = minmiss(x)
%MINMISS Column minimum with missing data.
%  Y = MINMISS(X) returns the smallest element in each column of X.
%  Missing data values must be encoded as NaNs. For vectors, MINMISS(X)
%  returns the smallest value of the elements in X.
%  [Y,I] = MINMISS(X) stores the indices of the minimum values in vector I.
%
%  See also MAXMISS, MEANMISS.

%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/03/08 14:27:31 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995
% added help function	G.Krahmann, IfM Kiel, Jun 1996

if nargin==0
  help minnan
  return
end

if ~isempty(x)
  [m,n] = size(x);
  y = Inf * ones(m,n);
  valid = ~isnan(x);
  y(valid) = x(valid);
  [y,i] = min(y);
else
  y=nan;
end
