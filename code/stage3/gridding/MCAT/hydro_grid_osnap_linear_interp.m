function [SGfs,TGfs,p_grid] = hydro_grid_osnap_linear_interp(p_hydrogrid)
%% HYDRO_GRID_OSNAP_LINEAR_INTERP *function [SGfs,TGfs,p_grid] = hydro_grid_osnap_linear_interp(p_hydrogrid)*
%% 
% * convert C into S and TEOS-10
% * lowpass filter
% * interpolate onto regular grid (pressure x time)
%% 
% Requires:  _*rodbload.m, sal78.m, ddspike.m, auto_filt.m, tem2sal.m, theta2sal.m, 
% con_tprof0.m, igrep.m*_
% 
% _*cmocean*_ <https://uk.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps 
% https://uk.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps>
% 
% _*gsw toolbox*_ <http://www.teos-10.org/software.htm http://www.teos-10.org/software.htm>
% 
% Author list: T. Kanzow, Nov 2005,  GDM, Sept 2011,  Loic Houpert October 2015, 
% Lewis Drysdale, 2020
warning off
% Load all the structure paths and paramter options
use p_hydrogrid
%  Load data
info       = [info_path,moor,'info.dat'];
[sn_info,id_info,wd,lt,ln,start_date,start_time,end_date,end_time] = ...
    rodbload(info,'serialnumber:instrument:Waterdepth:Latitude:Longitude:Start_Date:Start_Time:End_Date:End_Time');
start = [start_date' start_time(1)+start_time(2)/60];
stop  = [end_date' end_time(1)+end_time(2)/60];
% Read microcat data
% Loop option to bypass the calibration stage for shipboard processsing to look 
% at data
% if ~isempty(mc_int)
%     mc_path = [mooringpath,':',moor,':microcat:[',num2str(mc_int),']'];  
%     [yy_mc,mm,dd,hh,t,c,p,sn_mc,depth_mc] = ...
%         rodbload(mc_path,'yy:mm:dd:hh:t:c:p:SerialNumber:Instrdepth');
%     jd_mc     = julian(yy_mc,mm,dd,hh);
%     t_mc      = dum2nan(t,dum);
%     c_mc      = dum2nan(c,dum);
%     p_mc      = dum2nan(p,dum);
% else
    mc_path = [mooringpath,':',moor,':microcat:[',num2str(mc_ind),']'];  
    [yy_mc,mm,dd,hh,t,c,p,sn_mc,depth_mc] = ...
        rodbload(mc_path,'yy:mm:dd:hh:t:c:p:SerialNumber:Instrdepth');
    jd_mc     = julian(yy_mc,mm,dd,hh);
    t_mc      = dum2nan(t,dum);
    c_mc      = dum2nan(c,dum);
    p_mc      = dum2nan(p,dum);    
% end
% Interpolate T and C onto pressure time grid
jd_grid = ceil(julian(start)):1/iss:floor(julian(stop));
tres    = diff(jd_grid(1:2)); % get temporal resolution of new grid
T  = []; C = []; P = []; % set empty arrays
disp(' interpolate T and C onto time grid ...')
if exist('mc_ind','var')
    for inst = 1 : length(mc_ind)
        val  = find(yy_mc(:,inst) > 0 );
        T    = [T; interp1(jd_mc(val,inst),t_mc(val,inst),jd_grid)];
        C    = [C; interp1(jd_mc(val,inst),c_mc(val,inst),jd_grid)];
        P    = [P; interp1(jd_mc(val,inst),p_mc(val,inst),jd_grid)];
        instrdepth = [depth_mc];
        sn         = [sn_mc];
    end
end
% Repair pressures 
%% 
% # Where parts of timeseries are ok
[m,n]  = size(P);
P_nan   = isnan(P);
cnt_nan = sum(P_nan');
P_std   = nanstd(P');
for prep = 1 : m  
    if cnt_nan(prep) > 0 & cnt_nan(prep) < n & find(cnt_nan == 0)
        index_good = find(~isnan(P(prep,:)));
        index_bad  = find(isnan(P(prep,:)));
        comp       = nearest(prep,find(cnt_nan==0));
        comp       = comp(end);
        pol        = polyfit(P(comp,index_good),P(comp,index_good)-P(prep,index_good),2);
        val        = polyval(pol,P(comp,index_bad));
        pcorr      = P(prep,:);
        
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
%% 
% 2. Repair Pressure where complete time series is bad _(this loop will not 
% work if mc_P0 is not set in ini/*moor*.m script. Could add an error if mc_P0 
% is empty?)_
for prep = 1 : m 
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
% Converting C into S and convert TO TEOS 10
disp('converting C into S ...')
% ----------------  ------------------------------------
SP      = gsw_SP_from_C(C,T,P);         % Practical salinity from conductivity
[SA, ~]  = gsw_SA_from_SP(SP,P,ln,lt);   % Absolute salinity from practical salinity
T       = gsw_CT_from_t(SA,T,P);
% Despike
for i = 1 : m 
    if ~isempty(find(~isnan(SA(i,:)), 1))
        %   Despike salinity, replace bad values with NAN
        [SA(i,:),dx,~] = ddspike(SA(i,:),y_tol,stddy_tol,[nloop],'y',NaN);
        % Replce contemperaneous temperatures with NAN
        T(i,dx)         = NaN;
        % save de-spike plots
        ylabel('S_{A} (g kg^{-1})');
        title(['intrument' num2str(i)]);
        savename=[basedir '/Figures/' moor '/despike_instrument' num2str(i)];
        print(figure(989),'-dpng',savename);
    end
    close('all')
end
% Temporal low pass filter
% identify NaNs in the interpolated data
tnan_sum = sum(isnan(T'));
snan_sum = sum(isnan(SA'));
pnan_sum = sum(isnan(P'));
% make NaN matrix for filtered data
Tf = NaN * ones(m,n);
Sf = NaN * ones(m,n);
Pf = NaN * ones(m,n);
for pg = 1 : m
    tnnan = find(~isnan(T(pg,:)));
    snnan = find(~isnan(SA(pg,:)));
    pnnan = find(~isnan(P(pg,:)));
    
    if length(tnnan) < 30  % The number must be at least 3 times the filter order
        tnnan = [];
    end
    if length(snnan) < 30
        snnan = [];
    end
    if length(pnnan) < 30
        pnnan = [];
    end
    
    % Filter temperature as long as there are more than 30 non-nan records
    if ~isempty(tnnan)      
        % filter Temperature
        Tf(pg,tnnan) = auto_filt(T(pg,tnnan),1/tres,co);
        % interpolate on to original grid
        Tf(pg,:)     = interp1(jd_grid(tnnan),Tf(pg,tnnan)',jd_grid)';
        if tnan_sum(pg)/iss > gap_max
            gapI   = gap_mark(T(pg,:),gap_max,iss);
            Tf(pg,gapI) = NaN;
        end
    end
    
    % Filter salinity as long as there are more than 30 non-nan records
    if ~isempty(snnan)
        Sf(pg,snnan)  = auto_filt(SA(pg,snnan),1/tres,co);
        % interpolate on to original grid
        Sf(pg,:)      = interp1(jd_grid(snnan),Sf(pg,snnan)',jd_grid)';
        %   Sf(pg,snnan)  = auto_filt(Sf(pg,snnan),1/tres,co);
        
        if snan_sum(pg)/iss > gap_max
            gapI        = gap_mark(SA(pg,:),gap_max,iss);
            Sf(pg,gapI) = NaN;
        end        
    end
    
    % Filter pressue as long as there are more than 30 non-nan records
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
% Create new time grid
jd   = ceil(julian(start)+2): 1/fss:floor(julian(stop)-2);
% interpolate filtered data on to new grid
Tfs      = interp1(jd_grid,Tf',jd)';
Sfs      = interp1(jd_grid,Sf',jd)';
Pfs      = interp1(jd_grid,Pf',jd)';
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
%% 
% not sure about this stage?
%         % remove also bad data in the orginal grid files:
%         %figure;surf(t_RTWB,RTEB_merg.Pfs,RTEB_merg.Tfs);view([0 -90])
%         S_RTEBini(P_RTWB<100,ibadEB1topmc)=nan;
%         T_RTEBini(P_RTWB<100,ibadEB1topmc)=nan;
%         
%% 
% Get depth mean across time
%         WB(1).T = nanmean(RTWB_merg.Tfs(WBmc1,:));
%         WB(1).S = nanmean(RTWB_merg.Sfs(WBmc1,:));
%         WB(3).T = nanmean(RTWB_merg.Tfs(WBmc3,:));
%         WB(3).S = nanmean(RTWB_merg.Sfs(WBmc3,:));
%         
%         EB(1).T = nanmean(Tfs(EBmc1,:));
%         EB(1).S = nanmean(Sfs(EBmc1,:));
%         EB(3).T = nanmean(Tfs(EBmc3,:));
%         EB(3).S = nanmean(Sfs(EBmc3,:));
%     
%% 
% Low pass filtering of the temperature and salinity
%         daynber = 10;
%         Fc=1/daynber;
%         
%         % low pass WB depth mean temperature
%         abc=find(~isnan(WB(1).T ));
%         WB(1).Tlp =nan*WB(1).T ;
%         WB(1).Tlp(abc)=auto_filt(WB(1).T(abc),1/median(diff(jd)),Fc);
%         
%         % low pass EB depth mean temperature
%         abc=find(~isnan(EB(1).T ));
%         EB(1).Tlp =nan*EB(1).T ;
%         EB(1).Tlp(abc)=auto_filt(EB(1).T(abc),1/median(diff(t_RTWB)),Fc);
%         
%         % low pass WB depth mean salinity
%         abc=find(~isnan(WB(1).S));
%         WB(1).Slp =nan*WB(1).S;
%         WB(1).Slp(abc)=auto_filt(WB(1).S(abc),1/median(diff(t_RTWB)),Fc);
%         
%         % low pass WB depth mean salinity
%         abc=find(~isnan(EB(1).S));
%         EB(1).Slp =nan*EB(1).S;
%         EB(1).Slp(abc)=auto_filt(EB(1).S(abc),1/median(diff(t_RTWB)),Fc);
%         
%% 
% Plot temperature
%         
%         figure;
%         plot(EB(1).Tlp,WB(1).Tlp,'+');
%         
%         topmc = table(WB(1).Tlp',EB(1).Tlp','VariableNames',{'WB1','EB1'});
%         mdl = fitlm(topmc,'linear');
%         
%         %[b,bint,r,rint,stats] = regress(EB(1).T',[ones(size(EB(1).T)); WB(1).T]');
%         b= mdl.Coefficients.Estimate;
%         figure;plot(EB(1).T,WB(1).T,'+');
%         hold on; plot(EB(1).T,b(1) + b(2)*WB(1).T,'x')
%         
%         EB(1).Tnew = b(1) + b(2)*WB(1).T;
%         
%         figure;
%         plot(jd,EB(1).T);hold on;
%         plot(jd,EB(1).Tnew);
%         
%% 
% Plot Salinity
%         
%         figure;
%         plot(EB(1).Slp,WB(1).Slp,'+');
%         
%         topmc = table(WB(1).Slp',EB(1).Slp','VariableNames',{'WB1','EB1'});
%         mdl = fitlm(topmc,'linear');
%         
%         %[b,bint,r,rint,stats] = regress(EB(1).S',[ones(size(EB(1).S)); WB(1).S]');
%         b= mdl.Coefficients.Estimate;
%         figure;plot(EB(1).S,WB(1).S,'+');
%         hold on; plot(EB(1).S,b(1) + b(2)*WB(1).S,'x')
%         
%         EB(1).Snew = b(1) + b(2)*WB(1).S;
%         
%         figure;
%         plot(jd,EB(1).S);hold on;
%         plot(jd,EB(1).Snew );
%         
%         
%         Pfs(EBmc1(3),ibadEB1topmc)=RTWB_merg.Pfs(WBmc1(3),ibadEB1topmc);
%         Sfs(EBmc1(3),ibadEB1topmc)=EB(1).Snew(ibadEB1topmc);
%         Tfs(EBmc1(3),ibadEB1topmc)=EB(1).Tnew(ibadEB1topmc);
%     else
%            warning('2016 Western boundary mooring must be processed first!!')
%     end
% else
%     disp('No correction for missing data applied')
% end
% Create figures
if  pg ==1
    figure;   plot( Tfs); title('Temp')
    figure;   plot( Sfs); title('Salinity')
    figure;   plot( Pfs)
    figure;   plot(Sfs,theta(Pfs,Tfs,Sfs,0),'r.'); hold on
    
else
    figure;
    contourf(Tfs); title('{\Theta}');
    axis ij
    c=colorbar;
    ylabel('Instrument number');
    ylabel(c,'\celsius C')
    savename=[basedir '/Figures/' moor '/T'];
    print(gcf, '-dpng',savename);
    
    figure;
    contourf(Sfs); title('S_{a}');
    axis ij
    c=colorbar;
    ylabel('Instrument number');
    ylabel(c,'g kg^{-1}')
    savename=[basedir '/Figures/' moor '/SA'];
    print(gcf, '-dpng',savename);
    
    figure;    
    plot(Pfs'); title('Pressure');
    ylabel('db')
    axis ij;
    savename=[basedir '/Figures/' moor '/P'];
    print(gcf, '-dpng',savename);
    
    figure;
    dens=gsw_rho(Sfs,Tfs,Pfs);
    contourf(dens);title('potential density')
    axis ij
    c=colorbar;
    ylabel('Instrument number');
    ylabel(c,'kg m^{3}')      
    savename=[basedir '/Figures/' moor '/rho'];
    print(gcf, '-dpng',savename);
    
end
% Vertical linear interpolation
pmin     = ceil(mmin(Pfs)/p_gridsize)*p_gridsize;
pmax     = floor(mmax(Pfs)/p_gridsize)*p_gridsize;
p_grid   = [pmin:p_gridsize:pmax]';
TGfs = nan(length(p_grid),length(jd));
SGfs = nan(length(p_grid),length(jd)); 
for ijj=1:length(jd)
    TGfs(:,ijj) = interp1(Pfs(:,ijj),Tfs(:,ijj),p_grid) ; 
    SGfs(:,ijj) = interp1(Pfs(:,ijj),Sfs(:,ijj),p_grid) ;        
end
    
% Save data
save([out_path outname],'Tfs','Sfs', 'Pfs', 'TGfs', 'SGfs', 'jd',...
    'p_grid','co', 'T', 'C', 'SA', 'P', 'jd_grid', 'Pf', 'Tf', 'Sf');
% Add interpolated data to graphics 
[~,n] = size(TGfs);
figure;
pt_Gfs=gsw_pt_from_CT(SGfs,TGfs);
pt=gsw_pt_from_CT(Sf,Tf);
Y=[floor(min(min(pt))):0.1:ceil(max(max(pt)))];
X=[floor(min(min(Sf))):0.1:ceil(max(max(Sf)))];
[X,Y]=meshgrid(X,Y);
z=gsw_rho(X,Y,0);
[h,c]=contour(X,Y,z);
clabel(h,c);
hold on
plot(SGfs,pt_Gfs,'.k');
hold on
plot(Sf',pt','.');
hold on
xli = get(gca,'Xlim');
yli = get(gca,'Ylim');
xlim(xli);ylim(yli);
grid on
xlabel('S_{A} (g kg^{-1})');
ylabel('CT (\celsius C)')
title('\Theta-S_{A}')
print([out_path filesep 'hydro_grid_' moor '_theta_s'],'-dpng')
savename=[basedir '/Figures/' moor '/_theta_s'];
print(gcf, '-dpng',savename);
    
% Plot anomalies of non-gridded data
datum = gregorian(jd);
monI  = find(datum(:,3)==1 & datum(:,4) == 0 & ~isodd(datum(:,2)));
datum = datenum(datum);
if size(Tf,1) >1
    
    ta    = TGfs - meannan(TGfs',2)'* ones(1,n);
    sa    = SGfs - meannan(SGfs',2)'* ones(1,n);
    
    Pfwithnan1 = Pf;
    Pfwithnan1(isnan(Tf))=nan;
    Pfwithnan2 = Pf;
    Pfwithnan2(isnan(Sf))=nan;   
    datumgrid = datenum(gregorian(jd_grid));    
    
    figure;
    clf
    subplot(2,1,1)
    contourf(datum,p_grid,ta,11)
    set(gca,'xtick',datum(monI))
    shading flat
    hold on
    plot(datumgrid,Pfwithnan1','k','Linewidth',.5)
    datetick('x',12,'keepticks')
%     cmocean('balance','pivot',0)
    colorbar
    set(gca,'ydir','reverse','FontSize',12,'layer','top','tickdir','out')
    xlim([min(datum) max(datum)])
    title([moor,' T anomalies'], 'Interpreter', 'none')
    ylabel('Pressure [dbar]')
    
    
    subplot(2,1,2)
    contourf(datum,p_grid,sa,11)
    set(gca,'xtick',datum(monI))
    shading flat
    hold on
    plot(datumgrid,Pfwithnan2','k','Linewidth',.5)
    datetick('x',12,'keepticks')
%     cmocean('delta','pivot',0)
    colorbar
    set(gca,'ydir','reverse','FontSize',12,'layer','top','tickdir','out')
    xlim([min(datum) max(datum)])
    title([moor,' S anomalies'], 'Interpreter', 'none')
    ylabel('Pressure [dbar]')
    
    orient landscape
    
    print([out_path filesep 'hydro_grid_' moor '_ta_sa'],'-dpng')
        
    savename=[basedir '/Figures/' moor '/TS_ANOM'];
    print(gcf, '-dpng',savename);
end
%% *Sub routines*
function gapI=gap_mark(vec,gap_max,iss)
[a,b]    = consec_nan(vec);
gap      = b/iss;
gapI     = find(gap>gap_max);
if ~isempty(gapI)
    gapI        = igrep(sort([a(gapI) a(gapI)+b(gapI)-1]));
else
    gapI        = [];
end
function []=use(x)
%USE  Copies structure fields into named variables in workspace.
%
%   USE STRUCT, where STRUCT is a structure, copies all fields of the form
%   STRUCT.X into variables named X.
%  
%   This is useful for handling multiple datasets with the same variable
%   names.  The structures can be then kept in memory and 'mapped' into
%   variables as needed.
%
%   See also MAKE, MATSAVE.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2014 J.M. Lilly --- type 'help jlab_license' for details    
clear str
str{1}    =['if ~isempty(' x '),'];
str{end+1}=['  ZZFNAMES=fieldnames(' x ');' ];
str{end+1}='  for ZZi=1:length(ZZFNAMES),';
str{end+1}=[' 	  eval([ZZFNAMES{ZZi}, ''=getfield(' x ',ZZFNAMES{ZZi});'']);'];
str{end+1}='  end;';
str{end+1}='else;';
str{end+1}='  disp([''Contains no data.'']);'; 
str{end+1}='end;';
str{end+1}='clear ZZi ZZFNAMES';
str=strs2sray(str);
evalin('caller',str)
function [row]=strs2sray(x)
%STRS2SRAY  Converts a cell array of strings into a string array /w returns
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details      
  
if ~iscell(x)
  xtemp=x;
  clear x
  x{1}=xtemp;
end
M=length(x);
for i=1:M
    n(i)=length(x{i});
end
N=max(n);
row=[];
for i=1:M
   row=[row,x{i},char(10)]; 
end