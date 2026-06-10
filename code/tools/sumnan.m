% function sn = sumnan(s,dim)
%
% sum of matrix s along dimension dim
% NaNs are discarded

function sn = sumnan(s,dim);


ii = find(isnan(s));
s(ii) = 0;

if nargin > 1
 sn = sum(s,dim);
else 
 sn = sum(s);
end

return
