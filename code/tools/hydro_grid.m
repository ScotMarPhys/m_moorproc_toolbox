%function [SGfs,TGfs,p_grid] = hydro_grid(moor)
%
% convert C into S, lowpass filter and 
% interpolate onto regular grid (pressure x time)
%
% uses rodbload.m, sal78.m, ddspike.m, auto_filt.m, tem2sal.m, theta2sal.m,
%      con_tprof0.m, igrep.m
%

% T. Kanzow, Nov 2005

function [SGfs,TGfs,p_grid] = hydro_grid(moor)

%Global variables defined in config file

global co fss iss p_gridsize v_interp TSclim Theta_S_clim
global mc_id mc_cproblem mc_p0 mc_tproblem
global ac_id ac_cproblem ac_p0 ac_tproblem
global s4_id s4_cproblem s4_p0 s4_tproblem
global sc_id sc_cproblem sc_p0 sc_tproblem
global sg_id sg_cproblem sg_p0 sg_tproblem sg_moor

warning off
%%%moor2 = 'eb1_2_200516_subsampled';


% paths

%ini_path    = '/platte/rapid/data/mooring/proc/hydro_grid/';
ini_path    = '/noc/ooc/rpdmoc/rapid/data/moor/proc/hydro_grid/';
%rodbpath      '/platte/rapid/data/mooring/'
rodbpath      '/noc/ooc/rpdmoc/rapid/data/moor/'
%rodbctdpath  =    '/platte/rapid/data/';
rodbctdpath  =    '/noc/ooc/rpdmoc/users/tok/rapid/data/';
info_path   = [rodbpath,'/proc/',moor,'/'];
%out_path    = '/platte/rapid/data/mooring/proc/hydro_grid/';
out_path    = '/noc/ooc/rpdmoc/users/jcol/work/';
outname     = [moor,'_grid.mat'];

%external_ctd        = '/platte/rapid/data/cruise/d279/ctd/';
%external_ctd        = '/noc/ooc/rpdmoc/users/tok/rapid/data/cruise/d279/ctd/';
%external_ctd_file   = 'd279_pos.mat';
external_ctd        = ['/noc/ooc/rpdmoc/d279/ctd/        ';
                       '/noc/ooc/rpdmoc/d344/ctd/rodb/   ';
                       '/noc/ooc/rpdmoc/d334/ctd/rodb/   ';
                       '/noc/ooc/rpdmoc/sj08/ctd/newctd/ ';
                       '/noc/ooc/rpdmoc/rb0901/ctd/rodb/ ';
                       '/noc/ooc/rpdmoc/oc459/ctd/rodb/  ';
                       '/noc/ooc/rpdmoc/d346/ctd/rodb/   '];
external_ctd_file   = ['d279_pos.mat  ';
                       'd344_pos.mat  ';
                       'd334_pos.mat  ';
                       'sj0804_pos.mat';
                       'rb0901_pos.mat';
                       'o459_pos.mat  ';
                       'd346_pos.mat  '];


% general settings 

 dum       = -9999.0000;
 c1535     = 42.914;
 t90_68    = 1.00024;  % convert its90 to its68 for cond. to sal. conversion
 mcat      = [332:337];
 int_step  = 10;       % vertical interpolation step
 preverse  = 4000; %4000 pressure level below whch deep temperature reversion may occur
 
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
   mc_path = ['../../../users/jcol:',moor,':microcat:[',num2str(mc_ind),']'];
   [yy_mc,mm,dd,hh,t,c,p,sn_mc,depth_mc] = ...
        rodbload(mc_path,'yy:mm:dd:hh:t:c:p:SerialNumber:Instrdepth');
   jd_mc     = julian(yy_mc,mm,dd,hh);
   t_mc      = dum2nan(t,dum);
   c_mc      = dum2nan(c,dum); 
   p_mc      = dum2nan(p,dum); 
 end   
 if exist('rb_ind')
 if ~isempty(rb_ind)
   rb_path = ['../../../users/jcol:',moor,':rbr:[',num2str(rb_ind),']'];
   [yy_rb,mm,dd,hh,t,c,p,sn_rb,depth_rb] = ...
       rodbload(rb_path,'yy:mm:dd:hh:t:c:p:SerialNumber:Instrdepth');
   jd_rb     = julian(yy_rb,mm,dd,hh);
   t_rb      = dum2nan(t,dum);
   c_rb      = dum2nan(c,dum); 
   p_rb      = dum2nan(p,dum);
 end
 end
 if exist('ac_ind')
 if ~isempty(ac_ind)
   ac_path = ['../../../users/jcol:',moor,':arg:[',num2str(ac_ind),']'];
   [yy_ac,mm,dd,hh,t,c,p,p_arg,sn_ac,depth_ac] = ...
       rodbload(ac_path,'yy:mm:dd:hh:tcat:c:pcat:p:SerialNumber:Instrdepth');
   jd_ac     = julian(yy_ac,mm,dd,hh);
   t_ac      = dum2nan(t,dum);
   c_ac      = dum2nan(c,dum); 
   p_ac      = dum2nan(p,dum);
  if isempty(find(~isnan(t_ac)))     % NO ARGOCAT BUT ARGONAUT
        [yy_ac,mm,dd,hh,t,c,p,sn_ac,depth_ac] = ...
       rodbload(ac_path,'yy:mm:dd:hh:t:c:p:SerialNumber:Instrdepth');
   
        jd_ac     = julian(yy_ac,mm,dd,hh);
        t_ac      = dum2nan(t,dum);
        c_ac      = NaN * ones(length(t_ac),1); 
        p_ac      = dum2nan(p,dum); 
   end
 
 end
 end
 if exist('sc_ind')
 if ~isempty(sc_ind)
   sc_path = ['../../../users/jcol:',moor,':seacat:[',num2str(sc_ind),']'];
   [yy_sc,mm,dd,hh,t,p,sn_sc,depth_sc] = ...
       rodbload(sc_path,'yy:mm:dd:hh:t:p:SerialNumber:Instrdepth');
   jd_sc     = julian(yy_sc,mm,dd,hh);
   t_sc      = dum2nan(t,dum);
   c_sc      = t_sc*NaN; 
   p_sc      = dum2nan(p,dum); 

 end
 end
 if exist('sg_ind')
 if ~isempty(sg_ind)
   if ~isempty(sg_moor)
       moor = sg_moor;
   end    
   sg_path = ['../../../users/jcol:',moor,':seagauge:[',num2str(sg_ind(1)),']'];
   [yy_sg,mm,dd,hh,t,p,sn_sg,depth_sg] = ...
       rodbload(sg_path,'yy:mm:dd:hh:t:p:SerialNumber:Instrdepth');
   jd_sg     = julian(yy_sg,mm,dd,hh);
   t_sg      = dum2nan(t,dum);
   c_sg      = t_sg*NaN; 
   p_sg      = dum2nan(p,dum); 

 end
 end
 if exist('s4_ind')
 if ~isempty(s4_ind)
   s4_path = ['../../../users/jcol:',moor,':s4:[',num2str(s4_ind),']'];
   [yy_s4,mm,dd,hh,t,p,sn_s4,depth_s4] = ...
       rodbload(s4_path,'yy:mm:dd:hh:t:p:SerialNumber:Instrdepth');
   jd_s4     = julian(yy_s4,mm,dd,hh);
   t_s4      = dum2nan(t,dum);
   c_s4      = t_s4*NaN; 
   p_s4      = dum2nan(p,dum); 
   
 end
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
  for inst = 1 : length(rb_ind)
     val  = find(yy_rb(:,inst) > 0 ); 
     T    = [T; interp1(jd_rb(val,inst),t_rb(val,inst),jd_grid)];
     C    = [C; interp1(jd_rb(val,inst),c_rb(val,inst),jd_grid)];
     P    = [P; interp1(jd_rb(val,inst),p_rb(val,inst),jd_grid)];
     instrdepth = [instrdepth depth_rb];
     sn         = [sn sn_rb];
 end
for inst = 1 : length(ac_ind)
     val  = find(yy_ac(:,inst) > 0 ) ;
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

for inst = 1 : length(sc_ind)
     val  = find(yy_sc(:,inst) > 0 );  
     T    = [T; interp1(jd_sc(val,inst),t_sc(val,inst),jd_grid)];
     C    = [C; interp1(jd_sc(val,inst),c_sc(val,inst),jd_grid)];
     P    = [P; interp1(jd_sc(val,inst),p_sc(val,inst),jd_grid) + sc_p0(inst,2)];
     instrdepth = [instrdepth depth_sc];
     sn         = [sn sn_sc];
end 

for inst = 1 : length(sg_ind)
     val  = find(yy_sg(:,inst) > 0 );  
     T    = [T; interp1(jd_sg(val,inst),t_sg(val,inst),jd_grid)];
     C    = [C; interp1(jd_sg(val,inst),c_sg(val,inst),jd_grid)];
     P    = [P; interp1(jd_sg(val,inst),p_sg(val,inst),jd_grid) + sg_p0(inst,2)];
     instrdepth = [instrdepth depth_sg];
     sn         = [sn sn_sg];
end 


 [m,n]  = size(P);
 
 % ------ repair P ------------------------------------

%%return 
 
P_nan   = isnan(P);   
cnt_nan = sum(P_nan');
P_std   = nanstd(P');

%if strcmp(moor,'mochab_2_367')
%    if mc_ind(9) > 100 
%      C(9,:) = NaN;
%      disp(['MicroCAT #',num2str(sn(9)),' : C = NaN gesetzt!!!'])
%    end
%end    

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
    if(~isempty(find(P(prep,index_bad) <= 0)))
      a = find(P(prep,index_bad) <= 0);
      P(prep,index_bad(a))=NaN;      
    end
    
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
  %   Tf(pg,tnnan) = auto_filt(Tf(pg,tnnan),1/tres,co);
      if tnan_sum(pg)/iss > gap_max
        gapI   = gap_mark(T(pg,:),gap_max,iss);
        Tf(pg,gapI) = NaN;
      end     
    end
    
   if ~isempty(snnan)  
     Sf(pg,snnan)  = auto_filt(S(pg,snnan),1/tres,co);
     Sf(pg,:)      = interp1(jd_grid(snnan),Sf(pg,snnan)',jd_grid)';
  %   Sf(pg,snnan)  = auto_filt(Sf(pg,snnan),1/tres,co);
     if snan_sum(pg)/iss > gap_max
        gapI        = gap_mark(S(pg,:),gap_max,iss);
        Sf(pg,gapI) = NaN;
      end     

   end
   if ~isempty(pnnan)  
     Pf(pg,pnnan)  = auto_filt(P(pg,pnnan),1/tres,co);
     Pf(pg,:)      = interp1(jd_grid(pnnan),Pf(pg,pnnan)',jd_grid)';
  %   Pf(pg,pnnan)  = auto_filt(Pf(pg,pnnan),1/tres,co);
     if pnan_sum(pg)/iss > gap_max
        gapI        = gap_mark(P(pg,:),gap_max,iss);
        Pf(pg,gapI) = NaN;
      end     

     
   end
   
 end
if strcmp(moor,'eb1_4_200619')
   
   a = input('Salinity shall be manipulated in an ad hoc fashion y/n','s');  
   if strcmp(a,'y') 
     Sf(19,:) = Sf(19,:) + .01;
     Sf(20,:) = Sf(20,:) + .01;
     Sf(21,:) = Sf(21,:) + .005;
     Sf(22,:) = Sf(22,:) + .02;
     Sf(23,:) = Sf(23,:) + .025;
     Sf(24,:) = Sf(24,:) + .025;
     disp('Manipulation carried out')
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
      if isempty(innan)
         disp(['Pressure level ',num2str(pg),': No salinities'])
          
      else   
         pol   = polyfit(Tfs(pg,innan),Sfs(pg,innan),1);
         Sfs(pg,inan) = polyval(pol,Tfs(pg,inan));
      end   
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

% if  ~isempty(noS) %pg ==1
 if  pg ==1
   figure(1);   plot( Tfs); title('Temp')
   figure(2);   plot( Sfs); title('Salinity')
   figure(3);   plot( Pfs)
   figure(4);   plot(Sfs,theta(Pfs,Tfs,Sfs,0),'r.'); hold on
   
 else
    figure(1);clf;contourf(Tfs); title('temp')
    figure(2);clf;contourf(Sfs); title('Salinity')
    figure(3);clf;plot(Pfs')
    figure(4);clf;%plot(Sfs,theta(Pfs,Tfs,Sfs,0),'k');hold on
    dens=sw_pden(Sfs,Tfs,Pfs,5000);
    figure(5);clf;contourf(dens);title('density')
   
  end
 
%% load CTD data
%  if exist([external_ctd,external_ctd_file]) == 2
%    eval(['load ',external_ctd,external_ctd_file])
%    for cnt = 1 :length(ctd_prof)
%      dis(cnt) = dist2([ lt ctd_lat(cnt)],[ln ctd_lon(cnt)]);
%    end
%
%    near = find(dis<200e3);
%    rodbpath(rodbpath,rodbctdpath) 
%    [ctd_p,ctd_t,ctd_s] = ...
%        rodbload(['cruise:d279:ctd:[',num2str(near),']'],'p:t:s');
%  end
ct_p=NaN*ones(6000,50);
ct_t=NaN*ones(6000,50);
ct_s=NaN*ones(6000,50);
ct_pp=NaN*ones(6000,50);
ct_tt=NaN*ones(6000,50);
ct_ss=NaN*ones(6000,50);
ct_ppp=NaN*ones(6000,50);
ct_ttt=NaN*ones(6000,50);
ct_sss=NaN*ones(6000,50);
  ne=0;

for zz=1:size(external_ctd_file,1)
    if exist([external_ctd(zz,1:max(findstr(external_ctd(zz,:),'/'))),external_ctd_file(zz,(~isspace(external_ctd_file(zz,:))))]) == 2
       eval(['load ',external_ctd(zz,1:max(findstr(external_ctd(zz,:),'/'))),external_ctd_file(zz,:)])
       clear cnt dis
       rodbpath(rodbpath,rodbctdpath) 
       for cnt = 1 :length(ctd_prof)
         dis(cnt) = dist2([ lt ctd_lat(cnt)],[ln ctd_lon(cnt)]);
       end
       clear near
       near = find(dis<200e3)
       distance = dis(near)/1000
   for qq=1:length(near)
      if ~isempty(strfind(external_ctd(zz,:),'d279'))
        ctd_file = [external_ctd(zz,1:max(findstr(external_ctd(zz,:),'/'))),external_ctd_file(zz,1:4),'_',sprintf('%3.3d',ctd_prof(near(qq))),'.ctd'];
        [ctd_p,ctd_t,ctd_s] = rodbload(ctd_file,'p:t:s');
 %       [ctd_p,ctd_t,ctd_s] = rodbload(['cruise:d279:ctd:[',num2str(ctd_prof(near)),']'],'p:t:s');
        hold on
     ct_pp(1:size(ctd_p,1),qq+ne)=ctd_p;
     ct_tt(1:size(ctd_p,1),qq+ne)=ctd_t;
     ct_ss(1:size(ctd_p,1),qq+ne)=ctd_s;        
        xli = get(gca,'Xlim');
        yli = get(gca,'Ylim');
%        if ~isempty(find(pn>dum))
          %%plot(ctd_s,theta(ctd_p,ctd_t,ctd_s,median(pn(val))),'g')
  %        plot(ctd_s,sw_ptmp(ctd_s,ctd_t,ctd_p,median(pn(val))),'g')
%        else     
 %         plot(ctd_s,ctd_t,'g')
%        end 
      elseif ~isempty(strfind(external_ctd(zz,:),'d346'))
     
    %  [ctd_p,ctd_t,ctd_s] = rodbload(['cruise:d279:ctd:[',num2str(ctd_prof(near)),']'],'p:t:s');
      ctd_file = [external_ctd(zz,1:max(findstr(external_ctd(zz,:),'/'))),external_ctd_file(zz,1:4),'_',sprintf('%3.3d',ctd_prof(near(qq))),'.ctd'];
       [ctd_pp,ctd_tt,ctd_ss] = rodbload(ctd_file,'p:t:s');
          % [ctd_p,ctd_t,ctd_s] = rodbload(['..:..:..:..:d279:ctd:[',num2str(ctd_prof(near)),']'],'p:t:s');
      hold on
      xli = get(gca,'Xlim');
      yli = get(gca,'Ylim');
     ct_ppp(1:size(ctd_pp,1),qq+ne)=ctd_pp;
     ct_ttt(1:size(ctd_pp,1),qq+ne)=ctd_tt;
     ct_sss(1:size(ctd_pp,1),qq+ne)=ctd_ss; 
      else
                  
      % for qq=1:length(near)
        ctd_file = [external_ctd(zz,1:max(findstr(external_ctd(zz,:),'/'))),external_ctd_file(zz,1:4),'_',sprintf('%3.3d',ctd_prof(near(qq))),'.ctd'];
        [c_p,c_t,c_s] = rodbload(ctd_file,'p:t:s');
        hold on
        xli = get(gca,'Xlim');
        yli = get(gca,'Ylim');
 %       if ~isempty(find(pn>dum))
          %%plot(ctd_s,theta(ctd_p,ctd_t,ctd_s,median(pn(val))),'g')
  %        plot(ctd_s,sw_ptmp(ctd_s,ctd_t,ctd_p,median(pn(val))),'b')
   %     else     
 %         plot(ctd_s,ctd_t,'b')
    %    end 
    
     ct_p(1:size(c_p,1),qq+ne)=c_p;
     ct_t(1:size(c_p,1),qq+ne)=c_t;
     ct_s(1:size(c_p,1),qq+ne)=c_s;
     %ne=ne+length(near);
       end 
    ne=ne+length(near);
      end
  end
end
  
  
  
  
  
if size(Tf,1) >1 
  
% ----- vertical interpolation ----------

% Deal with collapsed moorings on a mooring by mooring basis
% Grid all until time that pressure of instruments goes out of sequence.
% After collapse, grid only those instruments that remain in pressure order.

  if strcmp(moor,'wb2_2_200528')
     len  = length(jd);
     Pfs1 = Pfs(:,1:351);
     Pfs2 = Pfs(7:11,352:end);
     Tfs1 = Tfs(:,1:351);
     Tfs2 = Tfs(7:11,352:end);
     Sfs1 = Sfs(:,1:351);
     Sfs2 = Sfs(7:11,352:end);
     pmin2     = ceil(mmin(Pfs2)/p_gridsize)*p_gridsize;
     pmin1     = ceil(mmin(Pfs)/p_gridsize)*p_gridsize;
     pmax      = floor(mmax(Pfs)/p_gridsize)*p_gridsize;
     p_grid    = [pmin1:p_gridsize:pmax]';
     p_grid2   = [pmin2:p_gridsize:pmax]';
     nlev      = length(p_grid);
     nlev2      = length(p_grid2);
     [TGfs(:,1:351),SGfs(:,1:351)] = ...
         con_tprof0(Tfs1,Sfs1,Pfs1,p_grid,0*ones(1,size(Tfs1,2)),int_step,TSclim,preverse);
     [TGfs(nlev-nlev2+1:nlev,352:len),SGfs(nlev-nlev2+1:nlev,352:len)] = ...
         con_tprof0(Tfs2,Sfs2,Pfs2,p_grid2,0*ones(1,size(Tfs2,2)),int_step,TSclim,preverse);  
      TGfs = dum2nan(TGfs,0);
      SGfs = dum2nan(SGfs,0);
  elseif strcmp(moor,'wb4_4_200703')
     len  = length(jd);
     Pfs1 = Pfs(:,1:152);
     Pfs2 = Pfs(8:12,153:end);
     Tfs1 = Tfs(:,1:152);
     Tfs2 = Tfs(8:12,153:end);
     Sfs1 = Sfs(:,1:152);
     Sfs2 = Sfs(8:12,153:end);
     pmin2     = ceil(mmin(Pfs2)/p_gridsize)*p_gridsize;
     pmin1     = ceil(mmin(Pfs)/p_gridsize)*p_gridsize;
     pmax      = floor(mmax(Pfs)/p_gridsize)*p_gridsize;
     p_grid    = [pmin1:p_gridsize:pmax]';
     p_grid2   = [pmin2:p_gridsize:pmax]';
     nlev      = length(p_grid);
     nlev2      = length(p_grid2);
     [TGfs(:,1:152),SGfs(:,1:152)] = ...
         con_tprof0(Tfs1,Sfs1,Pfs1,p_grid,0*ones(1,size(Tfs1,2)),int_step,TSclim,preverse);
     [TGfs(nlev-nlev2+1:nlev,153:len),SGfs(nlev-nlev2+1:nlev,153:len)] = ...
         con_tprof0(Tfs2,Sfs2,Pfs2,p_grid2,0*ones(1,size(Tfs2,2)),int_step,TSclim,preverse);  
      TGfs = dum2nan(TGfs,0);
  elseif strcmp(moor,'mar1_6_200940')
     len  = length(jd);
     Pfs1 = Pfs(:,1:694);
     Pfs2 = Pfs(9:18,695:end);
     Tfs1 = Tfs(:,1:694);
     Tfs2 = Tfs(9:18,695:end);
     Sfs1 = Sfs(:,1:694);
     Sfs2 = Sfs(9:18,695:end);
     pmin2     = ceil(mmin(Pfs2)/p_gridsize)*p_gridsize;
     pmin1     = ceil(mmin(Pfs)/p_gridsize)*p_gridsize;
     pmax      = floor(mmax(Pfs)/p_gridsize)*p_gridsize;
     p_grid    = [pmin1:p_gridsize:pmax]';
     p_grid2   = [pmin2:p_gridsize:pmax]';
     nlev      = length(p_grid);
     nlev2      = length(p_grid2);
     [TGfs(:,1:694),SGfs(:,1:694)] = ...
         con_tprof0(Tfs1,Sfs1,Pfs1,p_grid,0*ones(1,size(Tfs1,2)),int_step,TSclim,preverse);
     [TGfs(nlev-nlev2+1:nlev,695:len),SGfs(nlev-nlev2+1:nlev,695:len)] = ...
         con_tprof0(Tfs2,Sfs2,Pfs2,p_grid2,0*ones(1,size(Tfs2,2)),int_step,TSclim,preverse);  
      TGfs = dum2nan(TGfs,0);
  else
    fprintf(1,'\n Carrying out vertical interpolation of mooring \n %s using climatology\n %s\n\n',moor,TSclim)

    [TGfs,SGfs] = con_tprof0(Tfs,Sfs,Pfs,p_grid,0*ones(1,size(Tfs,2)),int_step,TSclim,preverse);
  end
else
    TGfs = [];
    SGfs = [];
    Pfs = Pfs';
    Sfs = Sfs';
    Tfs = Tfs';

end    
%%[scon,dscon,dsi] = ts2tscon(Tfs,Sfs,Pfs,TGfs,p_grid,TSclim);

 % ---save data ------

 % eval(['save ',out_path,outname,' Tfs Sfs Pfs TGfs SGfs jd p_grid co T C S P jd_grid Pf Tf Sf ctd_s ctd_p ctd_t'])
   eval(['save ',out_path,outname,' Tfs Sfs Pfs TGfs SGfs jd p_grid co T C S P jd_grid Pf Tf Sf'])
% ------ add interpolated data to graphics -----------
    
   [m,n] = size(TGfs); 

   figure(4);%plot(ct_s,theta(ct_p,ct_t,ct_s,0),'r--');hold on
   plot(SGfs,theta(p_grid*ones(1,size(TGfs,2)),TGfs,SGfs,0),'m');hold on
    plot(Sf',theta(Pf,Tf,Sf,0)','.');
    if ~isempty(sg_ind)
         plot(Sfs(end,:)',theta(Pfs(end,:),Tfs(end,:),Sfs(end,:),0)','.')
    end
   hold on
    xli = get(gca,'Xlim');
    yli = get(gca,'Ylim');
       plot(ct_ss,theta(ct_pp,ct_tt,ct_ss,0),'g')
       plot(ct_sss,theta(ct_ppp,ct_ttt,ct_sss,0),'c--')
    plot(ct_s,theta(ct_p,ct_t,ct_s,0),'b--')
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
        
% -----------------------------
% sub routines 
% ----------------------------------
function gapI = gap_mark(vec,gap_max,iss)

   [a,b]    = consec_nan(vec); 
   gap      = b/iss;
   gapI     = find(gap>gap_max);
   if ~isempty(gapI)
     gapI        = igrep(sort([a(gapI) a(gapI)+b(gapI)-1]));
   else
     gapI        = [];   
   end
