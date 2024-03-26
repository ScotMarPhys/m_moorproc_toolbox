% guess temperatures between vertical grid points using a dtdp climatology
%
% function [innerT,innerP] = t_int0(t,s,p,time,step,cTg,cSg,cT,tTg,tT,tt)
%
% input:  t   -- discrete temperature Profile from sensors at one point of time
%                (must be column vector)
%         p   -- corresponding pressures of sensors [dbar]
%         time-- point of time: decimal month
%         step-- integation step size for dtdp method
%         cTg -- constant gradient [C/dbar] for : 2.3 < T < 18 
%         cT  -- corresponding temperature
%         tTg -- time dependend  temp. gradient [C/dbar] for: 18<T<26
%         tT  -- corresponding temperature
%         tt  -- corresponding time (decimal months)    
%
% output: innerT-- 
%         innerP-- 
%
% uses: spacing.m
%  
% T.Kanzow 4.4.00
%          6.4.00   nearest.m replaced to speed up code 

function [innerT,innerS,innerP] = t_int0(t,s,p,time,step,cTg,cSg,cT,tTg,tT,tt)

innerT = [];
innerP = [];
innerS = [];
%%innert1=[];
%%innert2=[];

ni = length(t) -1;  % number of seperate integrations needed

% --- dtdp method -------------------------------------
% now, between each two grid points there  an upward and a downward 
% integration of dT/dP (T) and dS/dT (T) is performed and the 
% temperature / salinity guess for values inbetween is calculated
% from a weighted average of both integrations 
%
 
%disp('creating continous temperature / salinity profile')

for i = 1 : ni,       % put-parts-together loop
  tgrad1 = []; % temperature gradient for upward integration 
  tgrad2 = []; % temperature gradient for downward integration
  sgrad1 = []; % sal gradient for upward integration 
  sgrad2 = []; % sal gradient for downward integration
 
  t1     = []; % temperature guess for upward integration
  t2     = []; % temperature guess for downward integration
  s1     = []; % salinity guess for upward integration
  s2     = []; % salinity guess for downward integration


  pb = p(i:i+1);      % pressure boundaries for i-th integration
  tb = t(i:i+1);      % corresponding temperatures
  sb = s(i:i+1);
  
  inc    = spacing(pb(1),pb(2),step); % values of integration increment  

  dinc1  = diff(inc);
  dinc2  = diff(fliplr(inc));

  w1     = 1 - abs(inc-pb(1)) / (pb(2)-pb(1));  % weight for upward integration
  w2     = 1 - abs(inc-pb(2)) / (pb(2)-pb(1));  % weight for downward integration

  t1(1)     = tb(1);  % upper boundary temp
  t2(1)     = tb(2);  % lower boundary temp
  s1(1)     = sb(1);  % upper boundary temp
  s2(1)     = sb(2);  % lower boundary temp
   
 %-- create continuous temp. profile between two successive sensors --
  
  for j = 1 : length(dinc1)+1,  % integration loop
    
    if t1(j) <= 9999 % --------- beneath seasonal thermocline ----------------

       
      if  t1(j)>=min(cT) & t1(j)<=max(cT)   
        tgrad1(j)= interp1(cT,cTg,t1(j),'*linear');
        sgrad1(j)= interp1(cT,cSg,t1(j),'*linear');
      elseif t1(j) < min(cT)
           [xx,II] = min(cT);
           tgrad1(j) = cTg(II);
           sgrad1(j) = cSg(II);
      elseif t1(j) > max(cT)
           [xx,II]   = max(cT);
           tgrad1(j) = cTg(II);
           sgrad1(j) = cSg(II);
      end     
       if j < length(dinc1)+1
         t1(j+1) = t1(j) + tgrad1(j) * dinc1(j); % new temperature at p+dp 
         s1(j+1) = s1(j) + sgrad1(j) * dinc1(j); % new salinity at p+dp 
       end

    else %-------- above seasonal thermocline ---------------------------

       [xYY,ii] = min(abs(tT-t1(j)));
        
       [xYY,jj] = min(abs(tt(ii)-time));
       tgrad1(j)= tTg(ii(jj));

       if j < length(dinc1)+1
          t1(j+1) = t1(j) + tgrad1(j) * dinc1(j); % new temperature at p+dp
       end
 
    end % end thermocline if
  
   if t2(j) <= 9999  %--------- beneath seasonal thermocline ---------------------

       %%tgrad2(j)= interp1(cT,cTg,t2(j),'*linear');
       %%sgrad2(j)= interp1(cT,cSg,t2(j),'*linear');

      if  t2(j)>=min(cT) & t2(j)<=max(cT)   
        tgrad2(j)= interp1(cT,cTg,t2(j),'*linear');
        sgrad2(j)= interp1(cT,cSg,t2(j),'*linear');
      elseif t2(j) < min(cT)
           [xx,II] = min(cT);
           tgrad2(j) = cTg(II);
           sgrad2(j) = cSg(II);
      elseif t2(j) > max(cT)
           [xx,II]   = max(cT);
           tgrad2(j) = cTg(II);
           sgrad2(j) = cSg(II);
      end     
       if j < length(dinc1)+1
         t2(j+1)  = t2(j) + tgrad2(j) * dinc2(j); % new temperature at p+dp
         s2(j+1)  = s2(j) + sgrad2(j) * dinc2(j); % new salinity at p+dp
       end

    else %--------- above seasonal thermocline -------------------------------


         [xYY,ii] = min(abs(tT-t2(j)));


       [xYY,jj] = min(abs(tt(ii)-time));

       tgrad2(j)= tTg(ii(jj));   % 

       if j < length(dinc1)+1
          t2(j+1) = t2(j) + tgrad2(j) * dinc2(j); % new temperature at p+dp
       end

    end % end thermocline  if
  
   
  end    % end of integration loop 
  
  T1 = w1.* (tb(1) + cumtrapz(inc,tgrad1));  % weighted upward integration
 
  T2 = w2.* fliplr(tb(2) + cumtrapz(fliplr(inc),tgrad2)); % weighted downward integration

  S1 = w1.* (sb(1) + cumtrapz(inc,sgrad1));  % weighted upward integration
 
  S2 = w2.* fliplr(sb(2) + cumtrapz(fliplr(inc),sgrad2)); % weighted downward integration
   
  
  innerT = [innerT;[T1+T2]'];
  innerS = [innerS;[S1+S2]'];
  innerP = [innerP;inc'];
  
  %%innert1 = [innert1;t1'];
  %%innert2 = [innert2;flipud(t2')];
  

end    % end of put-parts-togehter loop   
 
% ---- eliminate equal values of successive data segments -----

ii = find(diff(innerP)~=0);
ii = [ii; [ii(length(ii))+1]]; 
innerT = innerT(ii); % 
innerP = innerP(ii);
innerS = innerS(ii);
