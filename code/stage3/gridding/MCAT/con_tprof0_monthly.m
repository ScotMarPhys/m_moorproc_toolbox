%% GRIDDING CODE
%
% $Rev: 255 $
% $Author: bim $
% $Date: 2019-06-27 16:50:34 +0100 (Thu, 27 Jun 2019) $
%
% make from a vertival dicrete temperature and salinity profiles  quasi 
% continous ones by guessing the course of the profiles between the vertical 
% grid points using a geographically local dt/dp climatology 
% output data is then interpolated onto a given vertical grid
% the dt/dp climatology method is dicribed in:
%   William E. Johns et al:
%   The Kuroshio east of Taiwan: Moored transport observations from the WOCE
%   PCM-1 Array
%
%   March 2012
%   con_tprof0_v3_bim has been adapted to use a monthly climatology
%
% 
%       bim: Notes from bim (March 2012).
%	  
%	1) for a given discrete temperature, estimate the next measurement at pressure+increment
% 	   by finding the temperature gradient ** at that temperature (from climatology) 
%          and projecting it forwards a distance incrememnt.
%          Repeat for each pressure level between the discrete measurements.
%          Do it starting at the shallower measurement (downwards)
%          and then again starting at the deeper measurement (upwards).
%	2) for a given discrete salinity estimate the next measurement at pressure+increment
%          by finding the salinity gradient ** at that salinity (from climatology) 
%          and projecting it forwards a distance incrememnt.
%          Repeat for each pressure level between the discrete measurements.
%          Do it starting at the shallower measurement (downwards) 
%          and then again starting at the deeper measurement (upwards).
%	3) calculate weights based on the distance from the shallower measurement (downwards)
%          and distance from the deeper measurement (upwards).
%       4) Temperature profile between discrete measurements is estimated from
%          solving dT/Dp (interpolated from climatology, see ** above) over the required pressure
%          levels. Add in the start temperature and apply the weighting.
%          This is done in the downwards and upwards directions.
%          The temperature profile between the discrete measurements is the sum of the
%          upwards and downwards profiles.
%        5) Salinity is treated in the same way (see 4).
%
% function [t_con,s_con] = con_tprof0(tg,sg,pg,PG,time,int_step,TSclim,TS_CLIMATOLOGY,TS_CLIMATOLOGY_NAME)
%    
%   input:
%      tg       = mooring temperature time series (stored as rows in tg)
%                 rows are the microcats at different depths
%      sg       = mooring salinity time series (stored as rows in sg)  
%                 rows are the microcats at different depths
%      pg       = corresponding pressure time series, (stored as rows in pg)
%                 dummies in tg and pg
%                 have to be marked as NaNs! ! !    
%      PG       = vertical pressure grid onto which quasi continous profiles
%                 will be interpolated -- must be column vector 
%                 e.g. pg = 0:20:4820; % depths in 20dbar bins
%      time     = month -- must be row vector
%                 correspnding to tg
%                 note: all values passed to the function were zero in previous versions
%      int_step = integration step size [dbar] between grid points for 
%                 dp/dt method, if empty will be set to 20 dbar
%
%      TSclim   = path to the file containing the dT/dP and dS/dP climatology
%
%	TS_CLIMATOLOGY -- climatology, e.g. deep, slope or full region.
%       TS_CLIMATOLOGY_NAME -- climatology name, e.g. argo or hbase
%
%   output:
%      t_con and  s_con   = output temperature interpolated onto the pressure grid (PG)
%
%   uses t_int0.m (that one  needs further user defined functions),
%        t_bound0.m, sst_check.m 
%
% T.Kanzow 3.4.00
% 
% CODE HISTORY
% ------------
% v3:  25 May 2010 added t_int0.m and t_bound0.m as sub-functions, with
% preverse variable added.  Previous problems with multiple copies of
% different edits.  Pressure grid variable changed to stop integration at
% the upper MicroCAT at each time step.
% 1 July 2010 changed to require a minimum of 4 MCs rather than 2.
% seasonal: 27th March 2012 Ben Moat

function [t_con,s_con] = con_tprof0_monthly(tg,sg,pg,PG,time,int_step,TSclim,TS_CLIMATOLOGY,TS_CLIMATOLOGY_NAME)  

% bim: why bring the gridded pressure read in, then cleared and reset ?
clear PG
PG = [0: 20: 4820]';
p_gridsize = 20;

% ---------------------------------------------
% select the minimum pressure at each time step
pmin = ceil(min(pg)/p_gridsize) * p_gridsize;
% select the maximum pressure at each time step
pmax = floor(max(pg)/p_gridsize) * p_gridsize; % remove???????

% find the maximum depths used
%depths=unique(pmax);
%depths=depths(find(~isnan(depths)))
%         for kk=1:length(depths),depths_sum(kk)=sum(pmax==depths(kk));end
%depths_sum
%	index=find(pmax==800)
%        datestr(time(index(1)))
%        datestr(time(index(2)))
%pmax = 4820;     % floor(max(pg)/p_gridsize) * p_gridsize;  changed yet
%again!!
% ---------------------------------------------

% set the vertical pressure grid
p_grid = pmin: p_gridsize: pmax;

t_con = []; 
s_con = [];
interp.error = [];

% --- load lookup tables for dt/dp climatology --------
disp(['loading ' TSclim])
eval(['load ',TSclim])
display([TSclim])

% added sept 2013
	if strmatch(TS_CLIMATOLOGY_NAME,'hbase')
%	change array name to allow us to use gerrards new hydrobase climatology
%	TS_CLIMATOLOGY will be slope, mareast, marwest or wall.
        ['dsdp_' TS_CLIMATOLOGY '_hbase = dsdp_' TS_CLIMATOLOGY '_on_t_hbase;clear dsdp_' TS_CLIMATOLOGY '_on_t_hbase;']
        ['dtdp_' TS_CLIMATOLOGY '_hbase = dtdp_' TS_CLIMATOLOGY '_on_t_hbase;clear dtdp_' TS_CLIMATOLOGY '_on_t_hbase;']
        eval( ['dsdp_' TS_CLIMATOLOGY '_hbase = dsdp_' TS_CLIMATOLOGY '_on_t_hbase;clear dsdp_' TS_CLIMATOLOGY '_on_t_hbase;'] );
        eval( ['dtdp_' TS_CLIMATOLOGY '_hbase = dtdp_' TS_CLIMATOLOGY '_on_t_hbase;clear dtdp_' TS_CLIMATOLOGY '_on_t_hbase;'] );
	end

tTg = 0;
tT = 0;
tt = 0;

% ------------- check input --------------------------------
if isempty(int_step)   % default integration step size [dbar] for dtdp-method
    int_step = 20;
    disp(['setting default integration size to ' num2str(int_step) ' dbar']);
end

% why only check the first column ?
if diff(pg([1:2],1)) < 0 % low pressure must be on top of matrix
  disp('lowest pressure is not on the top of the matrix!')
  tg = flipud(tg);
  pg = flipud(pg);
  sg = flipud(sg);
end

if diff(PG(1:2)) < 0 
  PG = flipud(PG);
end

% ----- call t_int.m successively for every point of time to get a time series
% ----- of vertical quasi continuous temperature profiles from the discrete 
% ----- mooring data -------------------------------------------------------

for ti = 1 : size(tg,2)  % for each time step (time loop)
	clear cTg cSg cT

%set climatology based on the month of the measurement
%choose a column of the climatology which corresponds to the month

	eval(['cTg = dtdp_' TS_CLIMATOLOGY '_' TS_CLIMATOLOGY_NAME '(:,time(ti));'])
	cTg = cTg'; % reorder 

	eval(['cSg = dsdp_' TS_CLIMATOLOGY '_' TS_CLIMATOLOGY_NAME '(:,time(ti));'])
	cSg = cSg'; % reorder

	eval(['cT  = tgrid_' TS_CLIMATOLOGY_NAME ';']);

%	cT  = tgrid_argo; % (1x475)
  
	  fprintf(1,['\r time step: ',num2str(ti),' of ',num2str(size(tg,2)),'                 '])
% index of measurements in the time series (not nan's) 
% where there are measurements of temperature and presure
	  ii = find(~isnan(tg(:,ti)) & ~isnan(pg(:,ti)));

	  if length(ii) < 5  % not-enough-values-if
%   create two columns with length of PG -- presure grid
	    innerTi = ones(length(PG),1)*NaN;
	    innerSi = ones(length(PG),1)*NaN;
	  else
%   set temperature (T), Pressure (P) and Salinity (S)
%   as the non nan discrete measurements from the moorings at each time setp	
%   
	    T  = tg(ii,ti);
	    P  = pg(ii,ti);
	    S  = sg(ii,ti);

% estimte the temperature between gridpoints 
% pass in T,S,P from the mooring at a time step ti

%****************************
    [innerT,innerS,innerP] = t_int0(T,S,P,time(ti),int_step,cTg,cSg,cT,tTg,tT,tt, 4000);
%****************************
% innerT -- temperature between gridpoints
% innerS -- Salinity between gridpoints
% innerP -- Pressure between gridpoints

%%%%%%%%%%%%%%%%%%%%%%%%%%
 % oneway integration for temp. guess above uppermost sensor
    if(P(1) > pmin(ti)) %  if the shallowest pressure measured is deeper than the minimum for that depth
                        %  interpolate to that depth using climatology.
                        %  This could happen if there are Nan's for the last few pressure values
         display('nans above shallow sensor')
%      [upT,upS,upP] = ...
%      t_bound0(T(1),S(1),P(1),time(ti),int_step,pmin,cTg,cSg,cT,tTg,tT,tt);% temperature above

      [upT,upS,upP] = ...
      t_bound0(T(1),S(1),P(1),time(ti),int_step,pmin(ti),cTg,cSg,cT,tTg,tT,tt);% temperature above
                                                        % uppermost sensor
      %-- don't allow temperatures > monthly mean SST, if exist, ------
      %-- replace by monthly mean SST ---------------------------------
          
      %%upT =  sst_check(upT,floor(time(ti)));  

      % --- put profile segments together ---------- 
      innerT = [upT;innerT];
      innerP = [upP;innerP];      
      innerS = [upS;innerS];  
    end    

%%%%%%%%%%%%%%%%%%%%%%%%%
% P(length(P)) %  -- lowest pressure measurement

%    if P(length(P)) < pmax %oneway integr. for temp. guess beneath downmost sensor
	if P(length(P)) < pmax(ti) %  if the deepest pressure is less than the maximum for that depth
				   %  interpolate to that depth using climatology.
				   %  This could happen if there are Nan's for the last few pressure values
	 display('nans below the deepest sensor')
%      [lowT,lowS,lowP] = ...
%      t_bound0(T(length(P)),S(length(P)),P(length(P)),time(ti),int_step,pmax,cTg,cSg,cT,tg,tT,tt);
      [lowT,lowS,lowP] = ...
      t_bound0(T(length(P)),S(length(P)),P(length(P)),time(ti),int_step,pmax(ti),cTg,cSg,cT,tg,tT,tt);

                                                        % temperature beneath
                                                        % deepest sensor
      % --- put profile segments together ---------- 
      innerT = [innerT;lowT];
      innerP = [innerP;lowP];
      innerS = [innerS;lowS];
    end

%%%%%%%%%%%%%%%%%
   
    % ---- eliminate equal values of successive data segments -----

    ii = find(diff(innerP)~=0);
    ii = [ii; [ii(length(ii))+1]];
   
    innerP = innerP(ii);
    innerT = innerT(ii);
    innerS = innerS(ii);
    
    %---- interpolate temperature output onto given grid ---------
% grid temperature and salinity onto the pressure grid in 20dbar levels
    innerTi= interp1(innerP,innerT,PG);
    innerSi= interp1(innerP,innerS,PG);
    
    %---- attach successive profile to one another as column vectors

  end % end not-enough-values-if
  
  % --- function output: time series of temperature profile
  % ---   
  t_con = [t_con, innerTi];
  s_con = [s_con, innerSi];
%  interp_error(ti).downward=uncertainty(:,1);
%  interp_error(ti).upward=uncertainty(:,2);
%  interp_error(ti).pressure=P;
  
end %  end of time loop

fprintf(1,'%\n')
end

%--------------------------------------------------------------------

%% SUB-function t_int0.m
% guess temperatures between vertical grid points using a dtdp climatology
%
% function [innerT,innerS,innerP,uncert] = t_int0(t,s,p,time,step,cTg,cSg,cT,tTg,tT,tt,preverse)
%
% input:  t   -- discrete temperature Profile from sensors at one point of time
%                (must be column vector)
%         s   -- discrete salinity Profile from sensors at one point of time
%                (must be column vector)
%         p   -- corresponding pressures of sensors [dbar]
%         time-- month
%         step-- integation step size for dtdp method (int_step (e.g. 10) passed to con_tprof0_v3.m)
%         cTg --  2.3 < T < 18 (dT/dP) from climatology
%         cT  -- corresponding temperature (tgrid)
%         cSg -- (dS/dP) from climatology
%         tTg -- time dependend  temp. gradient [C/dbar] for: 18<T<26  (passing in 0 )
%         tT  -- corresponding temperature (tgrid)                     (passing in 0 )
%         tt  -- corresponding time (decimal months)                   (passing in 0 )
%         preverse -- If a temperature inversion (i.e. increasing temperature with depth) is 
%                     found below preverse [dbar] instead of dT/dP and dS/dP a simple linear interpolation
%                     scheme is used to avoid ambiguities in the climatologies.%
%                   
%
% output: innerT-- Temperature profile between discreate data for a given time step 
%                  gridded onto the size 'step' 
%         innerS-- Salinity  profile between discreate data for a given time step 
%                  gridded onto the size 'step'
%         innerP-- Pressure profile between discreate data for a given time step 
%                  gridded onto the size 'step'
%
% uses: spacing.m
%  
% T.Kanzow 4.4.00
%          6.4.00   nearest.m replaced to speed up code 
%          29.06.07 input variable preverse added (see description avove) 

function [innerT,innerS,innerP,uncert] = t_int0(t,s,p,time,step,cTg,cSg,cT,tTg,tT,tt,preverse)

uncert=[];
innerT = [];
innerP = [];
innerS = [];

%%innert1=[];
%%innert2=[];

ni = length(t) -1;  % number of seperate integrations needed
                    % e.g. the number of gaps between the discrete temperature data

% --- dtdp method -------------------------------------
% now, between each two grid points there  an upward and a downward 
% integration of dT/dP (T) and dS/dT (T) is performed and the 
% temperature / salinity guess for values inbetween is calculated
% from a weighted average of both integrations 
%
 
%disp('creating continous temperature / salinity profile')

% step through each gap between measurements 
for i = 1 : ni,       % put-parts-together loop

  tgrad1 = []; % temperature gradient for upward integration 
  tgrad2 = []; % temperature gradient for downward integration
  sgrad1 = []; % sal gradient for upward integration 
  sgrad2 = []; % sal gradient for downward integration
 
  t1     = []; % temperature guess for upward integration
  t2     = []; % temperature guess for downward integration
  s1     = []; % salinity guess for upward integration
  s2     = []; % salinity guess for downward integration


% store two successive mooring pressure, temperature and salinity measurememnts
  pb = p(i:[i+1]);     % pressure boundaries for i-th integration
  tb = t(i:[i+1]);      % corresponding temperatures
  sb = s(i:[i+1]);      % corresponding salinities
  
%  Values of the pressures between pb(1) and pb(2) with a set spacing of step
% e.g. pb(1)=240, pb(2)=265, step = 10
% e.g. inc = [240 250 260 265]
  inc    = spacing(pb(1),pb(2),step); % values of integration increment  

  dinc1  = diff(inc); % difference between the inc (most will be step, but the end may not be)
  dinc2  = diff(fliplr(inc)); % put the last not step value first and make negative
    
% w1 decreases from 1 at pb(1) to zero at pb(2)
% e.g. pb(1)=50,pb(2)=100, inc = 50,60,70,80,90,100
% therefore, w1= [1,0.8,0.6,0.4,0.2,0]
  w1     = 1 - abs(inc-pb(1)) / (pb(2)-pb(1)); % weight for upward integration
                                               % pb(1) weight = 1
% w2 decreases from 1 at pb(2) to zero at pb(1)
  w2     = 1 - abs(inc-pb(2)) / (pb(2)-pb(1));  % weight for downward integration
					       % pb(2) weight = 1

% t1 -- the interpolated temperature profile between the discrete observations 
%       for depth increments of step (dbar) (e.g. 10)
%       Moves from the shallower discrete measurement to the deeper measurement
% t2 -- the interpolated temperature profile between the discrete observations
%       for depth increments of step (dbar) (e.g. 10)
%       Moves from the deeper discrete measurement to the shallower measurement
% s1 -- the interpolated salinity profile between the discrete observations
%       for depth increments of step (dbar) (e.g. 10)
%       Moves from the shallower discrete measurement to the deeper measurement
% s2 -- the interpolated salinity profile between the discrete observations
%       for depth increments of step (dbar) (e.g. 10)
%
% set the upper and lower discrete measurements
  t1(1)     = tb(1);   % upper boundary temp
  t2(1)     = tb(2);   % lower boundary temp
  s1(1)     = sb(1);  % upper boundary temp
  s2(1)     = sb(2);  % lower boundary temp
   
 %-- create continuous temp. profile between two successive sensors --
  
% for deep ocean cases where temperature gradient reverses
% (i.e., increases with depth),
% a linear interpolation should be used because the
% climatology is ambigous
   if pb(1) > preverse & pb(2) > preverse & tb(2) > tb(1)  

     	  T_linint = interp1(pb,tb, inc); % interpolated temperatures between tb(1) and tb(2)
      	  S_linint = interp1(pb,sb, inc); % interpolated salinities between sb(1) and sb(2)
          uncert_inversion=zeros(length(T_linint(:)) , 2); % interpolation so no error

      innerT = [innerT;T_linint(:)]; % build up profile
      innerS = [innerS;S_linint(:)]; % build up profile
      innerP = [innerP;inc']; 
      uncert = [uncert;uncert_inversion];
      
      fprintf(1,'\r t_int0.m: deep T inversion: p1 = %3.3d ; p2 = %3.3d ; t1 = %3.2f ; t2 =%3.2f',...
                                                               round(pb(1)),round(pb(2)),t1,t2) 
   else            % for all other cases the dT/dP (T) interpolation should be used     
 
     for j = 1 : length(dinc1)+1  % integration loop between pb(1) and pb(2)
% UPPER BOUNDARY TEMPERATURE  
       if t1(j) <= 9999 % --------- beneath seasonal thermocline ----------------
%	the temperature will always be less than 9999 ??

%    for deep climatology  2.4<=cT<=24
%       cTg = dT/dP
         if  t1(j)>=min(cT) & t1(j)<=max(cT)   % cT is the climatology temperature range (tgrid)
% if the measured temperature is within the tgrid range
% then store the value of dT/Dp and dS/dT at that temperature

           tgrad1(j)= interp1(cT,cTg,t1(j),'*linear');
           sgrad1(j)= interp1(cT,cSg,t1(j),'*linear');
         elseif t1(j) < min(cT)
%		display(['measured temperature ' num2str(t1(j)) 'is below the tgrid range'])
% if the measured temperature is below the tgrid range
% store the the value of dT/Dp and dS/dT at the minimum temperature
           [xx,II] = min(cT);
           tgrad1(j) = cTg(II);
           sgrad1(j) = cSg(II);
         elseif t1(j) > max(cT)
            %   display(['measured temperature ' num2str(t1(j)) 'is above the tgrid range'])
% if the measured temperature is above the tgrid range
% store the the value of dT/Dp and dS/dT at the maximum temperature
           [xx,II]   = max(cT);
           tgrad1(j) = cTg(II);
           sgrad1(j) = cSg(II);
         end     

         if j < length(dinc1) + 1 
% give a temperature at depth j project it to the next point using the rate of change
% and the distance to the next depth
           t1(j+1) = t1(j) + tgrad1(j) * dinc1(j); % new temperature at p+dp 
           s1(j+1) = s1(j) + sgrad1(j) * dinc1(j); % new salinity at p+dp 
         end

       else %-------- above seasonal thermocline ---------------------------
	   display('t1: Above 9999 degC')
	
         [xYY,ii] = min(abs(tT-t1(j)));
         [xYY,jj] = min(abs(tt(ii)-time));
         tgrad1(j)= tTg(ii(jj));

         if j < length(dinc1)+1
           t1(j+1) = t1(j) + tgrad1(j) * dinc1(j); % new temperature at p+dp
         end
 
       end % end thermocline if for t1


%       use the difference between the last projected point and the mooring measurement as an
%	error or uncertainty measurement.
%[tb(1) tb(2) t1(1) t1(end) t2(1) t2(end)]

% LOWER BOUNDARY TEMPERATURE

       if t2(j) <= 9999  %--------- beneath seasonal thermocline for t2---------------------
                         %          t2 lower boundary temperature

       %%tgrad2(j)= interp1(cT,cTg,t2(j),'*linear'); % alreagy commented out. bim
       %%sgrad2(j)= interp1(cT,cSg,t2(j),'*linear'); % alreagy commented out. bim

         if  t2(j)>=min(cT) & t2(j)<=max(cT)   
% if the measured temperature is within the tgrid range
% then find the value of dT/Dp and dS/dT at that temperature
           tgrad2(j)= interp1(cT,cTg,t2(j),'*linear');
           sgrad2(j)= interp1(cT,cSg,t2(j),'*linear');
         elseif t2(j) < min(cT)
% if the measured temperature is below the tgrid range
% find the the value of dT/Dp and dS/dT at the minimum temperature
           [xx,II] = min(cT);
           tgrad2(j) = cTg(II);
           sgrad2(j) = cSg(II);
         elseif t2(j) > max(cT)
% if the measured temperature is above the tgrid range
% find the the value of dT/Dp and dS/dT at the maximum temperature
           [xx,II]   = max(cT);
           tgrad2(j) = cTg(II);
           sgrad2(j) = cSg(II);
         end     
         if j < length(dinc1)+1
% given a temperature at depth j, project it to the next point using the temperature gradient
% and the distance to the next depth.
% dinc2 is used so the temperature at the next depth is shallower
% i.e. moving towards the surface
           t2(j+1)  = t2(j) + tgrad2(j) * dinc2(j); % new temperature at p+dp
% given a salinity at depth j, project it to the next point using the salinity gradient
% and the distance to the next depth.
           s2(j+1)  = s2(j) + sgrad2(j) * dinc2(j); % new salinity at p+dp
         end

       else %--------- above seasonal thermocline for t2-------------------------------
         display('t2: Above 9999 degC')
	t2(j-1)
	t2(j)
	j
	tgrad2(j)
	dinc2(j)
	error('stop')

         [xYY,ii] = min(abs(tT-t2(j)));
         [xYY,jj] = min(abs(tt(ii)-time));

         tgrad2(j)= tTg(ii(jj));   % 

         if j < length(dinc1)+1
          t2(j+1) = t2(j) + tgrad2(j) * dinc2(j); % new temperature at p+dp
         end

       end % end thermocline  if for t2

    end  % end of integration loop 
%%%%%%%%%%%%%%%%%%%%%%%%
% solve dT/dp=f(p) to get the temperature profile between the upper and lower bounds 
% 
    T1 = w1.* (tb(1) + cumtrapz(inc,tgrad1));  % weighted upward integration
    T2 = w2.* fliplr(tb(2) + cumtrapz(fliplr(inc),tgrad2)); % weighted downward integration

    S1 = w1.* (sb(1) + cumtrapz(inc,sgrad1));  % weighted upward integration
    S2 = w2.* fliplr(sb(2) + cumtrapz(fliplr(inc),sgrad2)); % weighted downward integration
   
%----
% added as in some situations there can be nans to below preverse
% so innerT is empty. bim Jan 2013
% if pressure >= preverse and there is an index
        if isempty(innerT)
             preverse_temp = t(1);  % no innerT so use the last temperature value
        else
             preverse_temp = innerT(end);
        end

%    if pb(1) >= preverse & ~isempty(find(T1 + T2 > innerT(end)))
    if pb(1) >= preverse & ~isempty(find(T1 + T2 > preverse_temp ))
%----
%  if pressure is below the depth preserve and the interpolated temperature profile
%  between pb(1) and pb(2) [ T1 + T2 ] has a value greater than the whole prtofile
%  Deep Temperature Overshoot due to dT/dP

       T_linint = interp1(pb,tb, inc); % interpolated temperatures between tb(1) and tb(2)
       S_linint = interp1(pb,sb, inc); % interpolated temperatures between sb(1) and sb(2)
       
      innerT = [innerT;T_linint(:)];
      innerS = [innerS;S_linint(:)];
      innerP = [innerP;inc']; 
      fprintf(1,'\r t_int0.m:  Deep Temperature Overshoot due to dT/dP - using linear interpolation instead') 
    else    
% sum T1 and T2 and add to full profile
% sum S1 and S2 and add to full profile
      innerT = [innerT;[T1+T2]']; % build up the temperature profile
      innerS = [innerS;[S1+S2]']; % build up the temperature profile
      innerP = [innerP;inc'];     % build up the pressure profile
    
    end
  end  % end distiction between normal (dT/dP) and linear interpolation (for deep temperature inversion) 
  

end    % end of put-parts-togehter loop   
 
% ---- eliminate equal values of successive data segments -----
% if two pressures the same them remove the second pressure,temp and salinity
% index of non sucessive repeats is ii
	ii = find(diff(innerP)~=0);
% add in a last value which is ii(end)+1
	ii = [ii; [ii(length(ii))+1]]; 
% 
	innerT = innerT(ii); % temperature between the two measurements with no repeats
	innerP = innerP(ii); % Pressure between the two measurements with no repeats
	innerS = innerS(ii); % salinity between the two measurements with no repeats
end

%% SUB FUNCTION t_bound0

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
%          p_bound -- pressure, up to where to perform the dtdp integration (pmax)
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
%P
%p_bound
%int_step
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


end
