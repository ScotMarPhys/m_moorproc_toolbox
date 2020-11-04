function [m,i,j]=mmax(a)
%MMAX Matrix Maximum Value. (MM)
% MMAX(A) returns the maximum value in the matrix A.
%
% [M,I] = MMAX(A) or [M,R,C] = MMAX(A) in addition returns the
% indices of ALL maximum values in the rows of I = [rows cols]
% or in R and C respectively.

% D.C. Hanselman, University of Maine, Orono ME 04469
% 1/4/95, revised 5/20/96
% Mastering MATLAB, Prentice Hall, ISBN 0-13-191594-0

m=max(max(a));
if nargout>=2	%return indices
	[i,j]=find(m==a);
	if nargout==2,i=[i,j];end
end
