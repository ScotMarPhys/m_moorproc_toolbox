% Guesses temperature / salinity profile between two pressures with by integrating
% a dtdp climatology with temperature given at the starting point. 
% Temperature is not allowed to take higher values than the mean 
% sst from the corresponding month 
%
%function [bdT,bdP] = t_bound0(T,S,P,time,int_step,p_bound,cTg,cSg,cT,tTg,tT,tt);
%
% input:   T/S      --- temperature/salinity of boundary
%          P      --- pressure of boundary 
%          time   --- corresponding time (decimal month)       
%          int_step-- integration step size 
%          p_bound -- pressure, up to where to perform the dtdp integration
%          cTg/cSg -- time constant temperature/salinity gradient   
%          cT      -- corresponding temperature
%          tTg     -- time dependend temperature gradient 18 <t<26
%          tT      -- corresponding temperature
%          tt      -- corresponding time (decimal month)
%
% output:  bdT     -- temperature profile between P and b_bound
%          bdP     -- corresponding pressure
%
% uses: spacing.m, sst_check.m
%
% T.Kanzow 4.4.00
%          6.4.00 nearest.m replaced to speed up code
 
function [bdT,bdS,bdP] = t_bound0(T,S,P,time,int_step,p_bound,cTg,cSg,cT,tTg,tT,tt);

fprintf(1,' oneway dtdp integration')
inc   = spacing(P,p_bound,int_step); 
%keyboard
dinc1  = diff(inc);
t1(1) = T;
s1(1) = S;

for j = 1:length(dinc1)+1,

  if t1 <= 9999              % beneath seasonal thermocline
      %%   tgrad1(j)= interp1(cT,cTg,t1(j)); 
      %%   sgrad1(j)= interp1(cT,cSg,t1(j)); 

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
            s1(j+1) = s1(j) + sgrad1(j) * dinc1(j); % new temperature at p+dp
  
      end
 
  else % above seasonal thermocline 
%%        [I,ii]   = nearest(t1(j),tT);    %look up nearest temp.
%%%        ii       = find(abs(tT-t1(j))==min(abs(tT-t1(j)))); % find nearest temp. in lookup table
        [xYY,ii] = min(abs(tT-t1(j)));
%%        [I,jj]   = nearest(time,tt(ii)); %look up nearest point of time
%%%        jj       = find(abs(tt(ii)-time)==min(abs(tt(ii)-time))); % find nearest time in lookup table
%%%        jj = jj(1);                      
        [xYY,jj] = min(abs(tt(ii)-time));
        tgrad1(j)= tTg(ii(jj));           

        if j < length(dinc1)+1
           t1(j+1) = t1(j) + tgrad1(j) * dinc1(j); % new temperature at p+dp
         end

  end % end if
end

bdT = (t1(1) + cumtrapz(inc,tgrad1))';
bdS = (s1(1) + cumtrapz(inc,sgrad1))';
bdP = inc'; 

if (P-p_bound) > 0   % sort temperature corresponding to increasing pressures
  bdT = flipud(bdT);
  bdS = flipud(bdS);
  bdP = flipud(bdP);
end


