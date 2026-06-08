function [y,i] = maxmiss(x)
%MAXMISS Column maximum with missing data.
%  Y = MAXMISS(X) returns the largest element in each column of X.
%  Missing data values must be encoded as NaNs. For vectors, MAXMISS(X)
%  returns the largest value of the elements in X.
%  [Y,I] = MAXMISS(Y) stores the indices of the maximum values in vector I.
%
%  See also MEANMISS, MINMISS.

%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/03/08 14:27:31 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995
% added help function	G.Krahmann, IfM Kiel, Jun 1996

if nargin==0
  help maxnan
  return
end

if ~isempty(x)
  [m,n] = size(x);
  y = -Inf * ones(m,n);
  valid = ~isnan(x);
  y(valid) = x(valid);
  [y,i] = max(y);
else
  y=nan;
end
