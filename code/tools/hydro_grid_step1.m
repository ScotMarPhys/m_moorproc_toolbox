function [Sfs, Tfs, Pfs, jd] = hydro_grid_step1(moor_list,inst_list)
% 
% function [Sfs,Tfs,p_grid,jd] = hydro_grid_step1(moor_list,inst_list)
% 
% function [Sfs,Tfs,p_grid,jd] = hydro_grid_step1(['eb1_9_201113 ...],[4710 ...])
%  
% GDM - JC064 - Sept 2011 - quick and dirty access to the .use files in the
% proc directory for a quick look at the uncalibrated MOC

 basedir      = '/noc/users/pstar/';
%basedir = '/Volumes/';

proc_path     = [basedir 'rpdmoc/rapid/data/moor/proc/'];

% general settings

dum       = -9999.0000;
c1535     = 42.914;
t90_68    = 1.00024;      % convert its90 to its68 for cond. to sal. conversion
mcat      = [332:337];
int_step  = 10;           % vertical interpolation step
preverse  = 4000;         % 4000 pressure level below whch deep temperature reversion may occur

% repair / despike settings

gap_max    = 10;  % allow for a maximum of gap [days] in data to interpolate across
y_tol      = [-10 10];
stddy_tol  = 4;
[nloop]    = 5;
graphics   = 'y';

% usually defined in an initialisation file... which we're not going to
% bother with

iss = 12;
fss = 2;
co  = 0.5;


for k = 1: length(inst_list)
    mc_path = [proc_path moor_list{k} '/microcat/' moor_list{k} '_' num2str(inst_list(k)) '.use'];
    [yy,mm,dd,hh,t,c,p,sn,depth] = ...
        rodbload(mc_path,'yy:mm:dd:hh:t:c:p:SerialNumber:Instrdepth');
    
    eval(['mc' num2str(inst_list(k)) '.jd = julian(yy,mm,dd,hh);'])
    eval(['mc' num2str(inst_list(k)) '.t = t;'])
    eval(['mc' num2str(inst_list(k)) '.c = c;'])
    eval(['mc' num2str(inst_list(k)) '.p = p;'])

end;

% mc_info = sn_info(find(id_info<=mcat(end) & id_info>=mcat(1)));

% ----- interpolate T and C onto pressure time grid

jstart = inf; jend = 0;
for k = 1: length(inst_list)
    
    local_jd = eval(['mc' num2str(inst_list(k)) '.jd']);
    if min(local_jd)<jstart;
        jstart = min(local_jd);
    end;
    if max(local_jd)>jend
        jend = max(local_jd);
    end;
    
end;
    
% jd_grid = ceil(julian(start)):1/iss:floor(julian(stop));
jd_grid = jstart:1/iss:jend;

tres    = diff(jd_grid(1:2)); % 

T  = []; C = []; P = [];
disp(' interpolate T and C onto time grid ...')
for k = 1 : length(inst_list)
%     val  = find(yy_mc(:,inst) > 0 );
    
    local_jd = eval(['mc' num2str(inst_list(k)) '.jd']);
    local_t  = eval(['mc' num2str(inst_list(k)) '.t']);
    local_c  = eval(['mc' num2str(inst_list(k)) '.c']);
    local_p  = eval(['mc' num2str(inst_list(k)) '.p']);
    
%     T    = [T; interp1(jd_mc(val,inst),t_mc(val,inst),jd_grid)];
%     C    = [C; interp1(jd_mc(val,inst),c_mc(val,inst),jd_grid)];
%     P    = [P; interp1(jd_mc(val,inst),p_mc(val,inst),jd_grid)];

    T    = [T; interp1(local_jd,local_t,jd_grid)];
    C    = [C; interp1(local_jd,local_c,jd_grid)];
    P    = [P; interp1(local_jd,local_p,jd_grid)];

    
%     instrdepth = [depth_mc];
%     sn         = [sn_mc];
end

% ------  temporal low pass filter ------

[m,n]  = size(P);

% ------ repair P ------------------------------------

% P_nan   = isnan(P);
% cnt_nan = sum(P_nan');
% P_std   = nanstd(P');
% 
% for prep = 1 : m  % repair Pressures where parts of timeseries are ok
%     if cnt_nan(prep) > 0 & cnt_nan(prep) < n & find(cnt_nan == 0)
%         index_good = find(~isnan(P(prep,:)));
%         index_bad  = find(isnan(P(prep,:)));
%         comp       = nearest(prep,find(cnt_nan==0));
%         comp       = comp(end);
%         pol        = polyfit(P(comp,index_good),P(comp,index_good)-P(prep,index_good),2);
%         
%         val               = polyval(pol,P(comp,index_bad));
%         pcorr             = P(prep,:);
%         P(prep,index_bad) = P(comp,index_bad)-val;
%         if(~isempty(find(P(prep,index_bad) <= 0)))
%             a = find(P(prep,index_bad) <= 0);
%             P(prep,index_bad(a))=NaN;
%         end
%         
%     end
% end
% 
% P_nan   = isnan(P);
% cnt_nan = sum(P_nan');
% P_std   = nanstd(P');
% 
% for prep = 1 : m % repair Pressure where complete time series is bad
%     if cnt_nan(prep) == n & sum(cnt_nan == 0) > 0
%         p0I    = find(sn(prep)== mc_p0(:,1));
%         P0     = mc_p0(p0I,2);
%         if isempty(P0)
%             disp(['Starting pressure for sensor #',num2str(sn(prep)),' needs to be defined in ini file'])
%             return
%         end
%         comp = find(cnt_nan == 0);
%         
%         if  sum(cnt_nan == 0) > 1
%             [XX,I] = sort(abs(instrdepth(comp)-instrdepth(prep)));
%             I      = comp(I);
%             
%             
%             P2f  = P(I(2),:) - mean(P(I(2),:));
%             P1f  = P(I(1),:) - mean(P(I(1),:));
%             fac  = diff(instrdepth([I(1) prep]))/diff(instrdepth(I([1 2])));
%             
%         elseif sum(cnt_nan == 0) == 1
%             
%             P2f  = zeros(1,n);
%             P1f  = P(comp,:) - mean(P(comp,:));
%             fac  = diff([instrdepth([comp prep])])/diff([instrdepth(comp) wd]);
%         end
%         pol       = polyfit(P1f,P2f-P1f,1);
%         P(prep,:) = polyval(pol,P1f)*fac + P0 + P1f;
%         
%         
%         
%     end
%     if cnt_nan(prep) == n & sum(cnt_nan == 0) == 0
%         p0I    = find(sn(prep)== mc_p0(:,1));
%         P0     = mc_p0(p0I,2);
%         if isempty(P0)
%             disp(['Starting pressure for sensor #',num2str(sn(prep)),' needs to be defined in ini file'])
%             return
%         end
%         
%         P(prep,:) = P0;
%         
%     end
% end


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

% jd   = ceil(julian(start)+2): 1/fss:floor(julian(stop)-2);
jd   = jstart+2: 1/fss: jend-2;

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

Tfs      = interp1(jd_grid, Tf', jd)';
Sfs      = interp1(jd_grid, Sf', jd)';
Pfs      = interp1(jd_grid, Pf', jd)';

