function [m,i,j]=mmin(a)
%MMIN Matrix Minimum Value. (MM)
% MMIN(A) returns the minimum value in the matrix A.
% [M,I] = MMIN(A) or [M,R,C] = MMAX(A) in addition returns the
% indices of ALL minimum values in the rows of I = [rows cols]
% or in R and C respectively.

% D.C. Hanselman, University of Maine, Orono ME 04469
% 1/4/95, revised 4/20/96
% Mastering MATLAB, Prentice Hall, ISBN 0-13-191594-0

m=min(min(a));
if nargout>=2	%return indices
	[i,j]=find(m==a);
	if nargout==2,i=[i,j];end
end
