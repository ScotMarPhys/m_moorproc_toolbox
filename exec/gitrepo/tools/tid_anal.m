function [co,am,ph,fit,l,V,U,mest]=tid_anal(time,px,sigma,cutoff)

% function [co,am,ph,fit,l,V,U,mest]=tid_anal(time,px,sigma,cutoff)
%
%
% harmonic tidal analysis of px (given at time) 
%
% time in days (column vector)
% px column vector
% sigma: tidal frequencies to use (degree/hour)
% cutoff: relativ value for cutting off eigenvalues (SingularValueDecomposition)
%
% output:
% mest = V*L*U  SVD-matrix
% fit = co(1) + am*cos((sigma*time*24-ph)*pi/180)  fitted curve
% am,ph  amplitude,phase
% co  coefficients (offset;cos(sigma1),sin(sigma1), ...) 


time=time.*24;  
freq=(2*pi/360).*sigma;
A=ones(length(time),2*length(freq)+1);
A1=cos(time*freq'); A2=sin(time*freq');
A(:,2:2:2*length(freq)+1)=A1;
A(:,3:2:2*length(freq)+1)=A2;

%disp(' start svd analysis')
 [U,L,V]=svd(A,0);
%disp(' ready svd analysis')

si=size(L);
l=diag(L(1:si(1),:));
i=find(l./max(l)>cutoff);
l_i=max(i);
if l_i < length(l)
 disp([int2str(length(l)-l_i) ' von ' int2str(length(l)) '  Eigenwerten nicht verwendet'])
end

mest=V(:,1:l_i)*L(1:l_i,1:l_i)^(-1)*U(:,1:l_i)';
%keyboard
co=mest*px;
am=sqrt((co(2:2:2*length(freq)+1,:).^2)+(co(3:2:2*length(freq)+1,:).^2));
ph=atan2(co(3:2:2*length(freq)+1,:),co(2:2:2*length(freq)+1,:))*180/pi;
for i=1:length(ph)
 if ph(i)<0
 ph(i)=ph(i)+360;
 end
end

if nargout > 3
 fit=A*co;                            
else
 fit=[];
end
if nargout < 5 
 mest=[];
end
if nargout < 6
 V=[];
end
if nargout < 7
 L=[];
end
if nargout < 8
 U=[];
end






