% function num_legend(num,[pos],[digit],[handle])
%
%  generates plot legend from vector  
%
%     num    = vector of numbers in legend
%     pos    = position of legend within plot (see legend for details) 
%     digit  = number of digits for each legend entry (filled with '0's)
%     handle = plot handle
%
% kanzow, 23.10.01
% 03.05.02  -- plot handles can be handed over
% sent by email to scu on knorr on 18/5/05
% 19.10.15  -- modified by L. Houpert: pos can be a string of character
% indicating the legend location ('north', 'southwestoutside',...)

function num_legend(num,loc,digit,handle)

if nargin <4 
  handle = [];
end


if nargin <3 
  digit = 3;
end
if nargin <2 
  loc = 1;
end

if isempty(loc)
  loc = 1;
end
if isempty(digit)
  digit = 3;
end


fmt =['%',num2str(digit),'.',num2str(digit),'d'];
str  = num2str(num,fmt);

nnum = length(num); 
nstr = []; 

for i = 1:nnum,
 nstr = [nstr,'str(',num2str((i-1)*digit+1),':',num2str(i*digit),'),']; 
end
%keyboard
if ~isstr(loc)
  if isempty(handle)
    eval(['legend(',nstr,num2str(loc),')'])
  else
    eval(['legend(handle,',nstr,num2str(loc),')'])
  end
else
   if isempty(handle)
    eval(['legend(' nstr '''location'',''' loc ''')'])
  else
    eval(['legend(handle,' nstr '''location'',''' loc ''')'])
  end   
end
    
