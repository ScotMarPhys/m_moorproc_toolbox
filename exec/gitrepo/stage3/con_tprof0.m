% make from a vertival dicrete temperature and salinity profiles  quasi 
% continous ones by guessing the course of the profiles between the vertical 
% grid points using a geographically local dt/dp climatology 
% output data is then interpolated onto a given vertical grid
%
% the dt/dp climatology method is dicribed in: 
%   William E. Johns et al: 
%   The Kuroshio east of Taiwan: Moored transport observations from the WOCE
%   PCM-1 Array 
%
% function [t_con,s_con] = con_tprof0(tg,sg,pg,p_grid,time,int_step,TSclim,preverse)
%    
%   input:
%      tg       = mooring temperature time series (stored as rows in tg)
%      sg       = mooring salinity time series (stored as rows n sg)  
%      pg       = corresponding pressure time series, dummies in tg and pg
%                 have to be marked as NaNs! ! !    
%      p_grid   = vertical pressure grid onto which quasi continous profiles
%                 will be interpolated -- must be column vector 
%      time     = decimal months -- must be row vector
%                 correspnding to tg
%      int_step = integration step size [dbar] between grid points for 
%                 dp/dt method, if empty will be set to 20 dbar
%
%      TSclim   = dT/dP and dS/dP climatology
%
%      preverse = If a temperature inversion (i.e. increasing temperature with depth) is found below preverse [dbar]
%                 instead of dT/dP and dS/dP a simple linear interpolation
%                 scheme is used to avoid ambiguities in the climatologies.
%                 Default = 4000 dbar
%                 
%                 
%                       
%
%   output:
%      t_con / s_con   = output temperature interpolated onto p_grid   
%
%   uses t_int0.m (that one  needs further user defined functions),
%        t_bound0.m, sst_check.m 
%
% T.Kanzow 3.4.00
% 29.06.07 input variable preverse added (see description avove) 

function [t_con,s_con] = con_tprof0(tg,sg,pg,p_grid,time,int_step,TSclim,preverse)  


if nargin == 7
    preverse = 4000;
    fprintf(1,'con_tprof0.m: Default value preverse = %d is used \n',preverse)
end    
t_con = []; 
s_con = [];

% --- load lookup tables for dt/dp climatology --------
disp('loading dtdp climatology')

eval(['load ',TSclim])
cTg = dtdp;  % dT/dP(T)  [C / dbar] 
cSg = dsdp;  % dS/dP(T)  [1 / dbar] 
cT  = tgrid;  % corresponding temperature [C]


tTg = 0;
tT = 0;
tt = 0;
disp('done.')

% ------------- check input --------------------------------
if isempty(int_step)   % default integration step size [dbar] for dtdp-method
  int_step = 20;
end
if diff(pg(1:2,1)) < 0 % low pressure must be on top of matrix
  tg = flipud(tg);
  pg = flipud(pg);
  sg = flipud(sg);
end
if diff(p_grid(1:2)) < 0 
  p_grid = flipud(p_grid);
end

% ----- call t_int.m successively for every point of time to get a time series
% ----- of vertical quasi continuous temperature profiles from the discrete 
% ----- mooring data -------------------------------------------------------
pmax = max(p_grid);
pmin = min(p_grid);

for ti = 1 : size(tg,2),  % time loop
  %if ti == 455
  %    keyboard
  %end
  
  fprintf(1,['\r time step: ',num2str(ti),' of ',num2str(size(tg,2)),'                 '])

  ii = find(~isnan(tg(:,ti)) & ~isnan(pg(:,ti)));

  if length(ii) <2  % not-enough-values-if
    innerTi = ones(length(p_grid),1)*NaN;
    innerSi = ones(length(p_grid),1)*NaN;
  else
    T  = tg(ii,ti);
    P  = pg(ii,ti);
    S  = sg(ii,ti);
   
    [innerT,innerS,innerP] = t_int0(T,S,P,time(ti),int_step,cTg,cSg,cT,tTg,tT,tt,preverse,ti);  %temperature between grid points  

    if(P(1) > pmin) % oneway integration for temp. guess above uppermost sensor
   
      [upT,upS,upP] = ...
      t_bound0(T(1),S(1),P(1),time(ti),int_step,pmin,cTg,cSg,cT,tTg,tT,tt);% temperature above
                                                        % uppermost sensor
      %-- don't allow temperatures > monthly mean SST, if exist, ------
      %-- replace by monthly mean SST ---------------------------------
          
      %%upT =  sst_check(upT,floor(time(ti)));  

      % --- put profile segments together ---------- 
      innerT = [upT;innerT];
      innerP = [upP;innerP];      
      innerS = [upS;innerS];  
    end    

    if P(length(P)) < pmax %oneway integr. for temp. guess beneath downmost sensor

      [lowT,lowS,lowP] = ...
      t_bound0(T(length(P)),S(length(P)),P(length(P)),time(ti),int_step,pmax,cTg,cSg,cT,tg,tT,tt);
                                                        % temperature beneath
                                                        % deepest sensor
      % --- put profile segments together ---------- 
      innerT = [innerT;lowT];
      innerP = [innerP;lowP];
      innerS = [innerS;lowS];
    end
   
    % ---- eliminate equal values of successive data segments -----

    ii = find(diff(innerP)~=0);
    ii = [ii; [ii(length(ii))+1]];
    innerP = innerP(ii);
    innerT = innerT(ii);
    innerS = innerS(ii);
    
    %---- interpolate temperature output onto given grid ---------

    innerTi= interp1(innerP,innerT,p_grid);
    innerSi= interp1(innerP,innerS,p_grid);
 
    %---- attach successive profile to one another as column vectors

  end % end not-enough-values-if
  
  % --- function output: time series of temperature profile
    
  t_con = [t_con, innerTi];
  s_con = [s_con, innerSi];

end %  end of time loop

fprintf(1,'%\n')