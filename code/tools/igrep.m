% function ovec = igrep(ivec)
%
% convert input vector ivec = [a1 a2 a3 a4 ...]
% into output vector ovec = [a1:a2 a3:a4 ...]
%
% kanzow
% kanzow, 24.4.01 error corrected 

function ovec = igrep(ivec)

l = length(ivec);

if isodd(l) 

     disp('length of ivec must be even')
     return
end


b = reshape(ivec,2,l/2)';

b = num2str(b);
sb = size(b);

ii = findstr(b(1,:),' ');
jj = find(diff(ii) >1);
if isempty(jj) 
  jj = 0;
end
b(:,ii(jj+1)) = ':';
b(:,sb(2)+1) = ' ';

sb = size(b);

ovec = ['[',reshape(b',1,sb(1)*sb(2)),']'];
ovec = str2num(ovec);
return

