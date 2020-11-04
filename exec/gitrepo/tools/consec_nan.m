function [a,b] = consec_nan(vec)

lvec = length(vec);

inan = find(isnan(vec));

a =[];
b = [];

dinan = diff(inan);
%tmp = vec;
a(1) = inan(1);

cnt = 1;

cnt2 = 0;

i    = 1;


for i = 1 : length(dinan)
 
 if(dinan(i)) == 1
     cnt = cnt + 1;
 else
     cnt2      = cnt2 +1;
     b(cnt2)   = cnt;
     cnt       = 1;
     a(cnt2+1) = inan(i+1);   
 end    
end


b(cnt2+1) = cnt;
