function [Tfs,Sfs,Pfs,TGfs,SGfs] = extrap_missing_data(pathosnap,Tfs,Sfs,Pfs,TGfs,SGfs,JG,pgg)
%% EXTRAP_MISSING_DATA extrap_missing_data
% *Script performing additonal editing to the gridded mooring data particularly 
% in regards to the data loss at RTEB1 between Dec 2016 and May 2017*
% 
% _*Authors: May 2020 - L. Houpert confined from Durham, Lewis Drysdale, 2021*_
%% Load mooring data
%% Missing data from the upper microcat on RTEB1.
% reconstruction of EB1 top microcat for the Mar 17 - May 17 period from data 
% from the RTWB microcat 
    
% load western boundary data
wbfile=fullfile(pathosnap, '\data\moor\proc','hydro_grid_merged','RTWB_merg_linear_interp_2018.mat');
load(wbfile)
S_RTWB = RTWB_merg.SGfs;
T_RTWB = RTWB_merg.TGfs;
P_RTWB = RTWB_merg.PGfs;
t_RTWB = RTWB_merg.JG;
S_RTEBini = SGfs;
T_RTEBini = TGfs;
WBmc1 = [1 10 19];
EBmc1 = [1 8 16];
WBmc3 = [1 10 19]+2;
EBmc3 = [1 8 15]+2;
for ijk=1:length(EBmc1)
    EBmc1index(ijk,1)= min(find(~isnan(Tfs(EBmc3(ijk),:))));
    EBmc1index(ijk,2)= max(find(~isnan(Tfs(EBmc3(ijk),:))));    
end
% Set data to bad for 3rd deployment of upper microcat on EB1 once the top of the mooring line broke (28th March 2017) 
ibad=find(JG>datenum(2016,12,10) & JG<datenum(2017,5,10));
ibadEB1topmc = ibad(1): EBmc1index(3,2);
Sfs(EBmc1(3),ibadEB1topmc)=nan;
Tfs(EBmc1(3),ibadEB1topmc)=nan;
% Remove also bad data in the orginal grid files:
S_RTEBini(P_RTWB<100,ibadEB1topmc)=nan;
T_RTEBini(P_RTWB<100,ibadEB1topmc)=nan;
WB(1).T = nanmean(RTWB_merg.Tfs(WBmc1,:));
WB(1).S = nanmean(RTWB_merg.Sfs(WBmc1,:));
WB(3).T = nanmean(RTWB_merg.Tfs(WBmc3,:));
WB(3).S = nanmean(RTWB_merg.Sfs(WBmc3,:));
EB(1).T = nanmean(Tfs(EBmc1,:));
EB(1).S = nanmean(Sfs(EBmc1,:));
EB(3).T = nanmean(Tfs(EBmc3,:));
EB(3).S = nanmean(Sfs(EBmc3,:));
% Low pass filtering
daynber = 10;
Fc=1/daynber;
abc=find(~isnan(WB(1).T ));
WB(1).Tlp =nan*WB(1).T ;
WB(1).Tlp(abc)=auto_filt(WB(1).T(abc),1/median(diff(t_RTWB)),Fc);
abc=find(~isnan(EB(1).T ));
EB(1).Tlp =nan*EB(1).T ;
EB(1).Tlp(abc)=auto_filt(EB(1).T(abc),1/median(diff(t_RTWB)),Fc);
abc=find(~isnan(WB(1).S));
WB(1).Slp =nan*WB(1).S;
WB(1).Slp(abc)=auto_filt(WB(1).S(abc),1/median(diff(t_RTWB)),Fc);
abc=find(~isnan(EB(1).S));
EB(1).Slp =nan*EB(1).S;
EB(1).Slp(abc)=auto_filt(EB(1).S(abc),1/median(diff(t_RTWB)),Fc);
% Temperature
figure;
plot(EB(1).Tlp,WB(1).Tlp,'+');
topmc = table(WB(1).Tlp',EB(1).Tlp','VariableNames',{'WB1','EB1'});
mdl = fitlm(topmc,'linear');
b= mdl.Coefficients.Estimate;
figure;plot(EB(1).T,WB(1).T,'+');
hold on; plot(EB(1).T,b(1) + b(2)*WB(1).T,'x')
EB(1).Tnew = b(1) + b(2)*WB(1).T;
figure;
plot(t_RTWB,EB(1).T);hold on;
plot(t_RTWB,EB(1).Tnew);

% Salinity
figure;
plot(EB(1).Slp,WB(1).Slp,'+');
topmc = table(WB(1).Slp',EB(1).Slp','VariableNames',{'WB1','EB1'});
mdl = fitlm(topmc,'linear');
b= mdl.Coefficients.Estimate;

figure;plot(EB(1).S,WB(1).S,'+');
hold on; plot(EB(1).S,b(1) + b(2)*WB(1).S,'x')
EB(1).Snew = b(1) + b(2)*WB(1).S;
print(gcf,'-dpng',fullfile(pathosnap, '\data\moor\proc','hydro_grid_merged','otherfigure','TSplusLinear.png'));
figure;
plot(t_RTWB,EB(1).S);hold on;
plot(t_RTWB,EB(1).Snew );
print(gcf,'-dpng',fullfile(pathosnap, '\data\moor\proc','hydro_grid_merged','otherfigure','EBold_EBnew.png'));



Sfslm = Sfs;
Tfslm = Tfs;
Pfs(EBmc1(3),ibadEB1topmc)=RTWB_merg.Pfs(WBmc1(3),ibadEB1topmc);
Sfs(EBmc1(3),ibadEB1topmc)=EB(1).Snew(ibadEB1topmc);
Tfs(EBmc1(3),ibadEB1topmc)=EB(1).Tnew(ibadEB1topmc);
RTEB_merg2.JG   = JG;
RTEB_merg2.Pfs  = Pfs;
RTEB_merg2.Sfs  = Sfs;
RTEB_merg2.Tfs  = Tfs;
%% Regridding with the reconstructed top microcat data for 2017 (following original RAPID script method) 
% Order the matrices at every time step to avoid too many NaNs creeping in. 
% 2004 removed...
P_sort = NaN .* ones(size(Pfs)); T_sort = NaN .* ones(size(Tfs)); S_sort = NaN .* ones(size(Sfs));
j = 1;
for ii = 1: length(JG)
    [~, ix] = sort(Pfs(:, ii));
    P_sort(:,j) = Pfs(ix,ii);
    T_sort(:,j) = Tfs(ix,ii);
    S_sort(:,j) = Sfs(ix,ii);
    j = j + 1;
end
%% 
% Removing unused rows of the sorted matrices
Pfss = nan(size(Pfs));
Tfss = nan(size(Pfs));
Sfss = nan(size(Pfs));
i = 1; j = 1;
for i = 1: length(P_sort(:,1))
    ix = find(isnan(P_sort(i,:)));
    if length(ix) < length(JG)
        Pfss(j,:) = P_sort(i, :);
        Tfss(j,:) = T_sort(i, :);
        Sfss(j,:) = S_sort(i, :);
        j = j + 1;
    end
end
TGfs = nan(length(pgg),length(JG));
SGfs = nan(length(pgg),length(JG)); 
for ijj=1:length(JG)
    itok = find(~isnan(Tfss(:,ijj)));
    isok = find(~isnan(Sfss(:,ijj)));
    if length(itok)>1
        TGfs(:,ijj) = interp1(Pfss(itok,ijj),Tfss(itok,ijj),pgg) ; 
    end
    if length(isok)>1
        SGfs(:,ijj) = interp1(Pfss(isok,ijj),Sfss(isok,ijj),pgg) ;     
    end
end
% Regridding with the reconstructed top microcat data for the whole deployment (diagnostic)
Pfslm = Pfs;
% order the matrices at every time step to avoid too many NaNs creeping in
% 2004 removed....
P_sort = NaN .* ones(size(Pfslm)); T_sort = NaN .* ones(size(Tfslm)); S_sort = NaN .* ones(size(Sfslm));
j = 1;
for ii = 1: length(JG)
    [~, ix] = sort(Pfslm(:, ii));
    P_sort(:,j) = Pfslm(ix,ii);
    T_sort(:,j) = Tfslm(ix,ii);
    S_sort(:,j) = Sfslm(ix,ii);
    j = j + 1;
end
% removing unused rows of the sorted matrices
Pfslms = nan(size(Pfslm));
Tfslms = nan(size(Pfslm));
Sfslms = nan(size(Pfslm));
i = 1; j = 1;
for i = 1: length(P_sort(:,1))
    ix = find(isnan(P_sort(i,:)));
    if length(ix) < length(JG)
        Pfslms(j,:) = P_sort(i, :);
        Tfslms(j,:) = T_sort(i, :);
        Sfslms(j,:) = S_sort(i, :);
        j = j + 1;
    end
end
   
TGfslm = nan(length(pgg),length(JG));
SGfslm = nan(length(pgg),length(JG)); 
for ijj=1:length(JG)
    itok = find(~isnan(Tfslms(:,ijj)));
    isok = find(~isnan(Sfslms(:,ijj)));
    if length(itok)>1
        TGfslm(:,ijj) = interp1(Pfslms(itok,ijj),Tfslms(itok,ijj),pgg) ; 
    end
    if length(isok)>1
        SGfslm(:,ijj) = interp1(Pfslms(isok,ijj),Sfslms(isok,ijj),pgg) ;     
    end
end
%% Plots for visual check
c200 = [9:0.2:11];
xlimit = [datenum(2014,08,01) datenum(2018,08,01)];
figure('units','centimeters','paperposition',[0 0 40 25])
subplot(2,1,1)
hold on
[~,hf]=contourf(t_RTWB,-P_RTWB,T_RTEBini,c200);
hf.LineStyle ='none';
plot(t_RTWB,Pfs,'-k')
ylim([-600 -50])
colorbar
caxis([min(c200) max(c200)])
title('Original EB1 temperature field without bad data removed after March 2017')
xlim(xlimit)
datetick('keeplimits')
subplot(2,1,2)
hold on
[~,hf]=contourf(t_RTWB,-P_RTWB,TGfs,c200);
hf.LineStyle ='none';
plot(t_RTWB,-Pfs,'-k')
ylim([-600 -50])
colorbar
caxis([min(c200) max(c200)])
title('New EB1 temp. field used in the analysis (linear regression on upper RTWB1 microcat data for Mar-May 2017 ')
xlim(xlimit)
datetick('keeplimits')
print(gcf,'-dpng',fullfile(pathosnap, '\data\moor\proc','hydro_grid_merged','otherfigure','_EB1_corrected_temp.png'));