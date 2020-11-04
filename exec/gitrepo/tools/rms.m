function y = rms(x)
%RMS    root-mean square. For vectors, RMS(x) returns the standard
%       deviation.  For matrices, RMS(X) is a row vector containing
%       the root-mean-square of each column. The difference to STD is
%       that here the mean is NOT removed.  
%
%	See also STD,COV.

%       Uwe Send, IfM Kiel, Apr 1992
% added NaN handling   Gerd Krahmann, IfM Kiel, Oct 1993, Jun 1994
% removed bug in NaN handling   G.Krahmann, Aug 1994

[m,n] = size(x);

bad=find(isnan(x));
if ~isempty(bad)
  if (m == 1) + (n == 1)
    x=x(find(~isnan(x)));
    [m,n] = size(x);
    if isempty(x)
      y=nan;
      return
    end
  end
end    
  

if (m == 1) + (n == 1)
	m = max(m,n);
        y = norm(x);
	y = y / sqrt(m);
else
	y = zeros(1,n);
	for i=1:n
                y(i) = rms(x(:,i));
	end
end

