% replaces 'dum' by NaN's in the input matrix a
% function b = dum2nan(a,dum);
%
% T.Kanzow

function b = dum2nan(a,dum)


ii = find(a == dum);
a(ii) = NaN;
b = a;




