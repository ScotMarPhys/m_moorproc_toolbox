function [y,mi] = meanmiss(x,dim)
%MEANMISS Column mean with missing data.
%  Y = MEANMISS(X) returns the mean of each column of X as a row vector
%  where missing data values are encoded as NaNs. For vectors, MEANMISS(X)
%  returns the mean value of the elements in X.
%  Y = MEANMISS(X,DIM) averages over dimension DIM (only version 5 or higher)

%  This code has been suggested by Douglas M. Schwarz (schwarz@kodak.com) 
%  in the news group comp.soft-sys.matlab.
%
%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/03/08 14:27:31 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995
% added compatibility for version 5.	G.Krahmann, LODYC Paris, Jul 1997

% check for version
v=version;
if strcmp(v(1),'5') | strcmp(v(1),'6') | strcmp(v(1),'7')
  if nargin==1 
    dim = 1;
    if size(x,1)==1
      dim = 2;
    end
  end
else
  dim=0;
end

if ~isempty(x)
  if dim==0
    [m,n] = size(x);
    y = zeros(m,n);
    valid = ~isnan(x);
    y(valid) = x(valid);
    y = sum(y) ./ sum(valid);
    mi = sum(valid);
  else
    y = zeros(size(x));
    valid = ~isnan(x);
    y(valid) = x(valid);
    y = sum(y,dim) ./ sum(valid,dim);
    mi = sum(valid,dim);
  end
else
  y=nan;
end
