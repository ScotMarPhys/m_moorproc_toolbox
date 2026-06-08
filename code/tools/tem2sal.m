% for given continuous temperature profiles  a  t/s climatology is used get continuous 
% salinity profiles 
%
% function s = tem2sal(t,TSclim)
%
% input:  t   -- temp matrix,dummies should be marked as 'NaN'
%         TSclim -- name of climatology
%        
% output: s       -- salinity matrix

%
% T.Kanzow 5.4.00

function [si] = tem2sal(ti,TSclim)


% ------- load ts relationship -------------------------

eval(['load ',TSclim])

% ---- get corresponding salinities ------
fprintf(1,'\n get salinity\n')
[m,n] = size(ti);


% check input data

i1  = find(isnan(ti)); % look for dummies
i2  = find(ti > max(tgrid)); % look for temp. above climatology temp. range
i3  = find(ti < min(tgrid)); % look for temp. beneath climatology temp. range

ti(i1) = -9999; % set NaN to out of range temp,
               % which will produce NaN in s

for i = 1 : n,
 si(1:m,i) = interp1(tgrid,s,ti(:,i),'*linear'); % table lookup
end



[xxx,i] = min(tgrid);
[xxx,j] = max(tgrid);

si(i2) = s(j);    % set s from temp. above climatology to S(Tmax)
si(i3) = s(i);    % set s from temp. beneath climatology to S(Tmin)




