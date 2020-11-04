function y = cumsummiss(x)
%MINMISS Cumulative sum with missing data.
%  Y = CUMSUMMISS(X) returns the cumulative sums over each column of X.
%  Missing data values must be encoded as NaNs. For vectors, CUMSUMMISS(X)
%  returns the cumulative sum of elements in X.
%

%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/03/09 14:07:36 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995

if ~isempty(x)
  [m,n] = size(x);
  y = zeros(m,n);
  valid = ~isnan(x);
  y(valid) = x(valid);
  y = cumsum(y);
else
  y=nan;
end
