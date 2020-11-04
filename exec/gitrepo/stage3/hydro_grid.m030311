%function [SGfs,TGfs,p_grid] = hydro_grid(moor)
%
% convert C into S, lowpass filter and 
% interpolate onto regular grid (pressure x time)
%
% uses rodbload.m, sal78.m, ddspike.m, auto_filt.m, tem2sal.m, theta2sal.m,
%      con_tprof0.m, igrep.m
%

% T. Kanzow, Nov 2005

%function [SGfs,TGfs,p_grid] = hydro_grid(moor)
moor = 'mar2_4_200729'
%Global variables defined in config file

global co fss iss p_gridsize v_interp TSclim Theta_S_clim
global mc_id mc_cproblem mc_p0 mc_tproblem
global ac_id ac_cproblem ac_p0 ac_tproblem
global s4_id s4_cproblem s4_p0 s4_tproblem

warning off
%%%moor2 = 'eb1_2_200516_subsampled';


% paths

ini_path    = '/noc/ooc/rpdmoc/rapid/data/moor/proc/hydro_grid/';
rodbpath  =    '/noc/ooc/rpdmoc/rapid/data/moor'  % added an equals sign PW

info_path   = [rodbpath,'/proc/',moor,'/'];
%out_path    = '/noc/ooc/rpdmoc/users/tok/work/';
out_path    = '/noc/ooc/rpdmoc/rapid/data/moor/proc/hydro_grid/';  % PW
outname     = [moor,'_grid.mat'];

%external_ctd        = '/noc/ooc/rpdmoc/users/tok/rapid/data/cruise/d279/ctd/';
rodb_ctd            = '/noc/ooc/rpdmoc/rapid/data/'; % rodbpath for CTD PW

external_ctd        = '/noc/ooc/rpdmoc/rapid/data/cruise/d279/ctd/';
external_ctd_file   = 'd279_pos.mat';



% general settings 

 dum       = -9999;
 c1535     = 42.914;
 t90_68    = 1.00024;  % convert its90 to its68 for cond. to sal. conversion
 mcat      = [332:337];
 int_step  = 10;       % vertical interpolation step

 % repair / despike settings 
 
 gap_max    = 10;  % allow for a maximum of gap [days] in data to interpolate across  
 y_tol      = [-10 10];
 stddy_tol  = 4;
 [nloop]    = 5;
 graphics   = 'y';

 % -------- execute mooring ini file ----- 
 if exist(moor) ~= 2
    path(path,ini_path)
 end
 %eval([moor2])
 eval([moor])
 % ------ load data -----------------

info       = [info_path,moor,'info.dat']; 
 [sn_info,id_info,wd,lt,ln,start_date,start_time,end_date,end_time] = ...
     rodbload(info,'serialnumber:instrument:Waterdepth:Latitude:Longitude:Start_Date:Start_Time:End_Date:End_Time');
 
start = [start_date' start_time(1)+start_time(2)/60];
stop  = [end_date' end_time(1)+end_time(2)/60];
 
 if ~isempty(mc_ind)
   mc_path = ['proc:',moor,':microcat:[',num2str(mc_ind),']'];
   [yy_mc,mm,dd,hh,t,c,p,sn_mc,depth_mc] = ...
        rodbload(mc_path,'yy:mm:dd:hh:t:c:p:SerialNumber:Instrdepth');
   jd_mc     = julian(yy_mc,mm,dd,hh);
   t_mc      = dum2nan(t,dum);
   c_mc      = dum2nan(c,dum); 
   p_mc      = dum2nan(p,dum); 
 end    
 
 if ~isempty(ac_ind)
   ac_path = ['proc:',moor,':arg:[',num2str(ac_ind),']'];
   
   [yy_ac,mm,dd,hh,t,c,p,p_arg,sn_ac,depth_ac] = ...
       rodbload(ac_path,'yy:mm:dd:hh:tcat:c:pcat:p:SerialNumber:Instrdepth');
   jd_ac     = julian(yy_ac,mm,dd,hh);
   t_ac      = dum2nan(t,dum);
   c_ac      = dum2nan(c,dum); 
   p_ac      = dum2nan(p,dum); 

 end
 if ~isempty(s4_ind)
    disp('Baustelle')
 end

 mc_info = sn_info(find(id_info<=mcat(end) & id_info>=mcat(1)));
 
 % ----- interpolate T and C onto pressure time grid

 jd_grid = ceil(julian(start)):1/iss:floor(julian(stop));
 tres    = diff(jd_grid(1:2)); % zeitliche Aufloesung
 
 T  = []; C = []; P = [];
 disp(' interpolate T and C onto time grid ...')
 for inst = 1 : length(mc_ind)
     val  = find(yy_mc(:,inst) > 0 ); 
     T    = [T; interp1(jd_mc(val,inst),t_mc(val,inst),jd_grid)];
     C    = [C; interp1(jd_mc(val,inst),c_mc(val,inst),jd_grid)];
     P    = [P; interp1(jd_mc(val,inst),p_mc(val,inst),jd_grid)];
     instrdepth = [depth_mc];
     sn         = [sn_mc];
 end 
for inst = 1 : length(ac_ind)
     val  = find(yy_ac(:,inst) > 0 ); 
     T    = [T; interp1(jd_ac(val,inst),t_ac(val,inst),jd_grid)];
     C    = [C; interp1(jd_ac(val,inst),c_ac(val,inst),jd_grid)];
     P    = [P; interp1(jd_ac(val,inst),p_ac(val,inst),jd_grid)];
     instrdepth = [instrdepth depth_ac];
     sn         = [sn sn_ac];
end 

for inst = 1 : length(s4_ind)
     val  = find(yy_s4(:,inst) > 0 );  
     T    = [T; interp1(jd_s4(val,inst),t_s4(val,inst),jd_grid)];
     C    = [C; interp1(jd_s4(val,inst),c_s4(val,inst),jd_grid)];
     P    = [P; interp1(jd_s4(val,inst),p_s4(val,inst),jd_grid)];
     instrdepth = [instrdepth depth_s4];
     sn         = [sn sn_s4];
end 


 [m,n]  = size(P);
 
 % ------ repair P ------------------------------------

%%return 
 
P_nan   = isnan(P);   
cnt_nan = sum(P_nan');
P_std   = nanstd(P');



for prep = 1 : m  % repair Pressures where parts of timeseries are ok
  if cnt_nan(prep) > 0 & cnt_nan(prep) < n & find(cnt_nan == 0) 
    index_good = find(~isnan(P(prep,:)));
    index_bad  = find(isnan(P(prep,:)));
    comp       = nearest(prep,find(cnt_nan==0));
    comp       = comp(end); 
    pol        = polyfit(P(comp,index_good),P(comp,index_good)-P(prep,index_good),2);

    val               = polyval(pol,P(comp,index_bad));
    pcorr             = P(prep,:);
    P(prep,index_bad) = P(comp,index_bad)-val;
    
  end 
end

P_nan   = isnan(P);
cnt_nan = sum(P_nan');
P_std   = nanstd(P');

for prep = 1 : m % repair Pressure where complete time series is bad
  if cnt_nan(prep) == n & sum(cnt_nan == 0) > 0
     p0I    = find(sn(prep)== mc_p0(:,1));
     P0     = mc_p0(p0I,2);
     if isempty(P0)
       disp(['Starting pressure for sensor #',num2str(sn(prep)),' needs to be defined in ini file'])  
       return
     end  
     comp = find(cnt_nan == 0); 

    if  sum(cnt_nan == 0) > 1
      [XX,I] = sort(abs(instrdepth(comp)-instrdepth(prep)));
      I      = comp(I);
      
         
      P2f  = P(I(2),:) - mean(P(I(2),:));
      P1f  = P(I(1),:) - mean(P(I(1),:));
      fac  = diff(instrdepth([I(1) prep]))/diff(instrdepth(I([1 2])));

    elseif sum(cnt_nan == 0) == 1 

      P2f  = zeros(1,n);
      P1f  = P(comp,:) - mean(P(comp,:));
      fac  = diff([instrdepth([comp prep])])/diff([instrdepth(comp) wd]);
     end
    pol       = polyfit(P1f,P2f-P1f,1);
    P(prep,:) = polyval(pol,P1f)*fac + P0 + P1f;
    
  
      
  end
  if cnt_nan(prep) == n & sum(cnt_nan == 0) == 0  
     p0I    = find(sn(prep)== mc_p0(:,1));
     P0     = mc_p0(p0I,2);
     if isempty(P0)
       disp(['Starting pressure for sensor #',num2str(sn(prep)),' needs to be defined in ini file'])  
       return
     end  

    P(prep,:) = P0;
    
  end    
end  
     

% ------ converting C into S and despike--------------

 disp('converting C into S ...')

 S = sal78(P,T*t90_68,C,c1535,0);

 % ----- despike S -------------------------------------
 
  
 for i = 1 : m   % despike
    if ~isempty(find(~isnan(S(i,:)))) 
      [S(i,:),dx,ndx] = ddspike(S(i,:),y_tol,stddy_tol,[nloop],'y',NaN);
      T(i,dx)         = NaN;
    end 
 end


 % ------  temporal low pass filter ------

 jd   = ceil(julian(start)+2): 1/fss:floor(julian(stop)-2); 
 
 
 tnan_sum = sum(isnan(T'));
 snan_sum = sum(isnan(S'));
 pnan_sum = sum(isnan(P'));
 
 Tf = NaN * ones(m,n);
 Sf = NaN * ones(m,n);
 Pf = NaN * ones(m,n);
 
  for pg = 1 : m
    tnnan = find(~isnan(T(pg,:)));
    snnan = find(~isnan(S(pg,:)));
    pnnan = find(~isnan(P(pg,:)));
    
    if length(tnnan) < 30  % Anzahl muss mind. 3fache Filterordnung sein  
      tnnan = [];
    end
    if length(snnan) < 30 
      snnan = [];
    end
    if length(pnnan) < 30 
      pnnan = [];
    end
    if ~isempty(tnnan)
    
      Tf(pg,tnnan) = auto_filt(T(pg,tnnan),1/tres,co);
      Tf(pg,:)     = interp1(jd_grid(tnnan),Tf(pg,tnnan)',jd_grid)';
      if tnan_sum(pg)/iss > gap_max
        gapI   = gap_mark(T(pg,:),gap_max,iss);
        Tf(pg,gapI) = NaN;
      end     
    end
    
   if ~isempty(snnan)  
     Sf(pg,snnan)  = auto_filt(S(pg,snnan),1/tres,co);
     Sf(pg,:)      = interp1(jd_grid(snnan),Sf(pg,snnan)',jd_grid)';
     if snan_sum(pg)/iss > gap_max
        gapI        = gap_mark(S(pg,:),gap_max,iss);
        Sf(pg,gapI) = NaN;
      end     

   end
   if ~isempty(pnnan)  
     Pf(pg,pnnan)  = auto_filt(P(pg,pnnan),1/tres,co);
     Pf(pg,:)      = interp1(jd_grid(pnnan),Pf(pg,pnnan)',jd_grid)';
     if pnan_sum(pg)/iss > gap_max
        gapI        = gap_mark(P(pg,:),gap_max,iss);
        Pf(pg,gapI) = NaN;
      end     

     
   end
   
 end
 
  
 Tfs      = interp1(jd_grid,Tf',jd)';  
 Sfs      = interp1(jd_grid,Sf',jd)';
 Pfs      = interp1(jd_grid,Pf',jd)';

 pmin     = ceil(mmin(Pfs)/p_gridsize)*p_gridsize;
 pmax     = floor(mmax(Pfs)/p_gridsize)*p_gridsize;
 p_grid   = [pmin:p_gridsize:pmax]';  

% close "big" salinity gaps

for pg = 1 : m
    
    inan  = find(isnan(Sfs(pg,:)) & ~isnan(Tfs(pg,:)));
    if ~isempty(inan)
      innan = find(~isnan(Sfs(pg,:)) & ~isnan(Tfs(pg,:)));
      pol   = polyfit(Tfs(pg,innan),Sfs(pg,innan),1);
      Sfs(pg,inan) = polyval(pol,Tfs(pg,inan));
    end 
end  

% use T/S climatology to get salinities where the sal. record is completely useless

noS = find(sum(~isnan(Sfs'))==0 | sum(Sfs')==0);
 for i = noS
     fprintf(1,'No salinities in record %d - using Theta/S climatology',i)
     sguess       = tem2sal(Tfs(i,:),TSclim);  % first guess for salinities
     theta_guess  = theta(Pfs(i,:),Tfs(i,:),sguess,0); % guess for theta
     Sfs(i,:)     = theta2sal(theta_guess,Theta_S_clim); 
 end

 if  pg ==1 
   figure(1);   plot( Tfs)
   figure(2);   plot( Sfs)
   figure(3);   plot( Pfs)
   figure(4);   plot(Sfs,theta(Pfs,Tfs,Sfs,0),'r.'); hold on
   
 else
    figure(1);clf;contourf(Tfs)
    figure(2);clf;contourf(Sfs)
    figure(3);clf;plot(Pfs')
    figure(4);clf;%plot(Sfs,theta(Pfs,Tfs,Sfs,0),'k');hold on
   
  end
 
%% load CTD data
  if exist([external_ctd,external_ctd_file]) == 2
    rodbpath(rodb_ctd)  
    eval(['load ',external_ctd,external_ctd_file])
    for cnt = 1 :length(ctd_prof)
      dis(cnt) = dist2([ lt ctd_lat(cnt)],[ln ctd_lon(cnt)]);
    end

    near = find(dis<200e3);

    [ctd_p,ctd_t,ctd_s] = ...
        rodbload(['cruise:d279:ctd:[',num2str(near),']'],'p:t:s');
  end

if size(Tf,1) >1 
  
% ----- vertical interpolation ----------

  fprintf(1,'\n Carrying out vertical interpolation of mooring \n %s using climatology\n %s\n\n',moor,TSclim)

  [TGfs,SGfs] = con_tprof0(Tfs,Sfs,Pfs,p_grid,0*ones(1,size(Tfs,2)),int_step,TSclim);

else
    TGfs = [];
    SGfs = [];
    Pfs = Pfs';
    Sfs = Sfs';
    Tfs = Tfs';

end    
%%[scon,dscon,dsi] = ts2tscon(Tfs,Sfs,Pfs,TGfs,p_grid,TSclim);

 % ---save data ------

 eval(['save ',out_path,outname,' Tfs Sfs Pfs TGfs SGfs jd p_grid co T C S P jd_grid'])

% ------ add interpolated data to graphics -----------
    
   [m,n] = size(TGfs); 

   figure(4);plot(SGfs,theta(p_grid*ones(1,size(TGfs,2)),TGfs,SGfs,0),'m');hold on
    plot(Sf',theta(Pf,Tf,Sf,0)','.');

   hold on
    xli = get(gca,'Xlim');
    yli = get(gca,'Ylim');
    plot(ctd_s,theta(ctd_p,ctd_t,ctd_s,0),'g')
    xlim(xli);ylim(yli);
       grid on

   title('\theta-S')  
   eval(['print -depsc hydro_grid_',moor,'_theta_s.epsc'])     

   datum = gregorian(jd);
   monI  = find(datum(:,3)==1 & datum(:,4) == 0 & ~isodd(datum(:,2)));
   datum = datenum(datum);
   
   if size(Tf,1) >1 
       
    ta    = TGfs - meannan(TGfs,2)* ones(1,n);
    sa    = SGfs - meannan(SGfs,2)* ones(1,n);
   
   
   
    figure(5)
    clf
       subplot(2,1,1)
        contourf(datum,p_grid,ta,11)
        set(gca,'xtick',datum(monI))
        shading flat
        hold on
        plot(datum,Pfs','k','Linewidth',.5)
        datetick('x',12,'keepticks')
        colormap(jet(10))
        colorbar
        set(gca,'ydir','reverse','FontSize',12,'layer','top','tickdir','out')
        xlim([min(datum) max(datum)])
        title([moor,' T anomalies'])
        ylabel('Pressure [dbar]')
     
     
      subplot(2,1,2)
        contourf(datum,p_grid,sa,11)
        set(gca,'xtick',datum(monI))
        shading flat
        hold on
        plot(datum,Pfs','k','Linewidth',.5)
        datetick('x',12,'keepticks')
        colormap(jet(10))
        colorbar
        set(gca,'ydir','reverse','FontSize',12,'layer','top','tickdir','out')
        xlim([min(datum) max(datum)])
        title([moor,' S anomalies'])
        ylabel('Pressure [dbar]')
     
        orient landscape
     
        eval(['print -depsc hydro_grid_',moor,'_ta_sa.epsc'])     
   end 
  if 0      
% -----------------------------
% sub routines 
% ----------------------------------
%function gapI = gap_mark(vec,gap_max,iss)

 %  [a,b]    = consec_nan(vec); 
 %  gap      = b/iss;
 %  gapI     = find(gap>gap_max);
 %  if ~isempty(gapI)
 %    gapI        = igrep(sort([a(gapI) a(gapI)+b(gapI)-1]));
%   else
 %    gapI        = [];   
 %  end
   
end