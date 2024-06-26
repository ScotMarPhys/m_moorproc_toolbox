% mc_caldip_check reads microcat caldip data from raw rodb file
% and compares data with the lowered CTD data
%
% DAS created this file based on mc_call_caldip

% LOH May 2017 - adapted for SBE37 ODO and CTD with two oxygen sensors
% SCU June 2018 AR30-04 

global MOORPROC_G
clearvars -except MOORPROC_G
close all;

cruise = MOORPROC_G.cruise; % used for microcat data
YEAR = MOORPROC_G.YEAR;
% Time origin to turn jday into dday -- same function used by rodbload
jd0 = julian(YEAR,1,1,0);

% -----------------------------------------------------------------
% User supplied information
% Select cast number and choose which set of instruments from CTD
cast = input('Which cast number? ','s');
ctdnum = sprintf('%03d',str2double(cast));
ctdsen = input('Which CTD sensors (1 or 2 [or blank to use already-selected primary])?) ','s');
oxysen = input('Which CTD oxygen (1 or 2 [or blank to use already-selected primary])?' , 's');

% --- get paths for data input and output ---
pd = moor_inoutpaths('microcat_cal_dip',cast);
if ~exist(fileparts(pd.stage2log),'dir')
    warning('creating directory for log file')
    try
        mkdir(fileparts(pd.stage2log))
    end
end

% ----------------- load CTD DATA   ----------------------------------
cvars = 'time press temp1 cond1 oxygen1 temp2 cond2 temp cond oxygen ';
h = m_read_header(pd.ctdfile); if sum(strcmp(h.fldnam,'oxygen2')); cvars = [cvars 'oxygen2 ']; end
d = mload(pd.ctdfile,[cvars ' ']);
if strcmp(cruise,'d382')
    % Correction for Di382
    d.cond1=d.cond1*10;  
    d.cond2=d.cond2*10;  
end
%and rename primary(s)
d.temp = d.(['temp' ctdsen]);    
dnum = m_commontime(d,'time',h,'datenum');
d.dday = dnum - datenum(MOORPROC_G.YEAR,1,1);
d.cond = d.(['cond' ctdsen]);
d.oxygen = d.(['oxygen' oxysen]);
ctdtimetx = datestr(dnum(1));
ptittxt = sprintf('Cast %s start %s SBE-CTD sensor set %s',cast,ctdtimetx,ctdsen);


% Selecttime period that we will analyse
if strcmp(cruise,'en705') && strcmp(cast,'2')
    %due to microcat setup computer time error data logging didn't start
    %until upcast, so just use deepest stop
    mxpi = 3046;
    imp = d.press > mxpi-10 & d.press < mxpi+10 & d.dday>d.dday(d.press==max(d.press));
    pwarn = 1;
    lct = 'CTD stats for deepest 10 dbar where microcats were logging';
else
    %within 10 dbar of deepest
    mxpi = max(d.press);
    imp = d.press > mxpi-10;
    pwarn = 0;
    lct = 'CTD stats for deepeset 10 dbar';
end
trc = [min(d.dday(imp)) max(d.dday(imp))];
tm1p = trc(1)+diff(trc)/4;
tm2p = trc(2)-diff(trc)/4;
imp2 = d.dday > tm1p & d.dday < tm2p;
meanctdpr = nanmean(d.press(imp2));

% CTD stats during period
meanctdpr = nanmean(d.press(imp2));
stdctdpr = nanstd(d.press(imp2));
ctd1_cond_mn = nanmean(d.cond1(imp2));
ctd1_cond_st = nanstd(d.cond1(imp2));
ctd1_temp_mn = nanmean(d.temp1(imp2));
ctd1_temp_st = nanstd(d.temp1(imp2));
ctd1_oxy_mn = nanmean(d.oxygen1(imp2));
ctd1_oxy_st = nanstd(d.oxygen1(imp2));
ctd2_cond_mn = nanmean(d.cond2(imp2));
ctd2_cond_st = nanstd(d.cond2(imp2));
ctd2_temp_mn = nanmean(d.temp2(imp2));
ctd2_temp_st = nanstd(d.temp2(imp2));
ctd2_oxy_mn = nanmean(d.oxygen2(imp2));
ctd2_oxy_st = nanstd(d.oxygen2(imp2));
ctd_cond_mn = nanmean(d.cond(imp2));
ctd_cond_st = nanstd(d.cond(imp2));
ctd_temp_mn = nanmean(d.temp(imp2));
ctd_temp_st = nanstd(d.temp(imp2));
ctd_oxy_mn = nanmean(d.oxygen(imp2));
ctd_oxy_st = nanstd(d.oxygen(imp2));

% --- get mooring information from infofile ---
[zins,id,sn]= rodbload(pd.infofile,'z:instrument:serialnumber');

% --- vector of serial numbers ---
ii = find(id >= 332 & id <= 337);
vec = sn(ii);
id2=id(ii);
nvec = length(ii);
zmic = zins(ii);

% Open output file for text and set plot name
%pd.stage2log = [outpath,'microcat_check',cast,'.log'];
ilogf = fopen(pd.stage2log,'w');
if ilogf==-1
    error('could not open log file %s for writing',pd.stage2log)
end
if ~exist(fileparts(pd.stage2fig),'dir')
    warning('creating directory for figures')
    try
        mkdir(fileparts(pd.stage2fig))
    catch
        error('no directory for %s',pd.stage2fig)
    end
end

% --- read data loop --
for i = 1:nvec
    % display( [,num2str(vec(i)),])

    infile = fullfile(pd.stage1path,sprintf(pd.stage1form,vec(i)));
    % --- load rodb data ---
    if id2(i)~=335 % checks if not ODO microcat
        [yy,mm,dd,hh,c,t,p] = rodbload(infile,'yy:mm:dd:hh:c:t:p');
        o=c*nan;
    else
        [yy,mm,dd,hh,c,t,p,o] = rodbload(infile,'yy:mm:dd:hh:c:t:p:o2');
    end
    if (i > 7 && i<=14)  lstr='--'; elseif i>14  lstr ='-.'; else lstr = '-'; end
    % Time variable
    dday = julian(yy,mm,dd,hh)-jd0;
    % interpolate CTD onto microcat for a rough and ready mean diff
    pi = interp1(d.dday, d.press, dday);
    ti = interp1(d.dday, d.temp, dday);
    ci = interp1(d.dday, d.cond, dday);
    oi = interp1(d.dday, d.oxygen, dday);

    %% PLOT DATA AT BOTTLE STOPS
    dobotstop=1; 
    if dobotstop %***come back to this, only giving 0-1 dots per figure?
        dp=gradient(p(i,:));
        % i=find(dp<0.02 & dp>-0.02 & abs(ctdp(1,j)-p(1,j)) < 10);
        if id2(i)==337
            iis=find(dp<0.2 & dp>-0.2);
            nsubplot =3;
        elseif id2(i)==335
            iis=find(dp<1 & dp>-1);
            nsubplot =4;
        end
        figure(10+i);clf;
        subplot(nsubplot,1,1);hold on ; grid on;
        plot(pi(iis),pi(iis)-p(iis),'ko');
        xlim([0 max(pi)]); ylim([-5 5]);
        xlabel('press [db]');ylabel('dp [db]');
        title(['s/n : ',num2str(vec(i))]);

        % subplot(nsubplot,1,2);hold on ; grid on;
        % plot(ctdp(j,i),(ctdt2(j,i)-t(j,i)),'ko');
        % xlim([0 max(ctdp(j,:))]); ylim([-0.01 0.01]);
        % xlabel('press [db]');ylabel('dt [C]');
        % 
        % 
        % subplot(nsubplot,1,3);hold on ; grid on;
        % plot(ctdp(j,i),(ctdc2(j,i)-c(j,i)),'ko');
        % xlim([0 max(ctdp(j,:))]); ylim([-0.01 0.01]);
        % xlabel('press [db]');ylabel('dc [mS/cm]');
        % 
        % if id2(j)==335
        %     subplot(nsubplot,1,4);hold on ; grid on;
        %     plot(ctdp(j,i),(ctdo2(j,i)-o2(j,i)),'ko');
        %     xlim([0 max(ctdp(j,:))]); ylim([-20 0]);
        %     xlabel('press [db]');ylabel('do2 [umol/l]');
        % end

end

    % Select data from period when CTD kept at maximum depth
    impt = dday>tm1p & dday < tm2p;
    if sum(impt)==0
        warning('no data from %4.4d near pressure %d',vec(i),mxpi)
    end
    nimpt(i,1) = sum(impt);           % no of data selcted
    pstd(i,1) = nanstd(p(impt) - pi(impt));
    cstd(i,1) = nanstd(c(impt) - ci(impt));
    tstd(i,1) = nanstd(t(impt) - ti(impt));
    ostd(i,1) = nanstd(o(impt) - oi(impt));
    pdifx(i,1) = nanmean(p(impt) - pi(impt));
    cdifx(i,1) = nanmean(c(impt) - ci(impt));
    tdifx(i,1) = nanmean(t(impt) - ti(impt));
    odifx(i,1) = nanmean(o(impt) - oi(impt));
    vectx(i,:) = sprintf('%5i',vec(i));

    % DIfferences of CTD sensors for comparison
    if i == 1
        tix = interp1(d.dday, d.temp2-d.temp1, dday);
        cix = interp1(d.dday, d.cond2-d.cond1, dday);
        oix = interp1(d.dday, d.oxygen2-d.oxygen1, dday);
        c_ctd_m = nanmean(cix(impt));
        c_ctd_s = nanstd(cix(impt));
        t_ctd_m = nanmean(tix(impt));
        t_ctd_s = nanstd(tix(impt));
        o_ctd_m = nanmean(oix(impt));
        o_ctd_s = nanstd(oix(impt));
    end
    % Select data close to microcats nominal deployment depth
    pbin = 2;
    pbstep = 200;
    ipx = find(pi >zmic(i)-(pbstep/pbin) & pi < zmic(i)+pbstep/pbin);
    [nh,px] = hist(pi(ipx),floor(pbstep/pbin));
    nmhx = find(nh == max(nh));
    ptest(i) = px(nmhx(1));
    ipx2 = find(pi > ptest(i)-2*pbin & pi < ptest(i)+2*pbin);
    ipx3 = ipx2(5:end-2);
    nimpt(i,2) = length(ipx3);
    pstd(i,2) = nanstd(p(ipx3) - pi(ipx3)); % pi is ctd on microcat time base
    cstd(i,2) = nanstd(c(ipx3) - ci(ipx3));
    tstd(i,2) = nanstd(t(ipx3) - ti(ipx3));
    ostd(i,2) = nanstd(o(ipx3) - oi(ipx3));
    pdifx(i,2) = nanmean(p(ipx3) - pi(ipx3));
    cdifx(i,2) = nanmean(c(ipx3) - ci(ipx3));
    tdifx(i,2) = nanmean(t(ipx3) - ti(ipx3));
    odifx(i,2) = nanmean(o(ipx3) - oi(ipx3));
    %     disp(['proceeding to next file '])
end % for i = 1:length(vec)

disp('data loaded')
% Quick look all data and search for outliers
toutlie(1:nvec) = {' '};
if nvec~=1
    for kk = 1:3
        if kk == 1
            xxv = cdifx(:,1); % microcat - ctd
            sxv = cstd(:,1);
            txv = 'Conductivity';
            t2xv = 'C';t3xv = 'c';
        elseif kk == 2
            xxv = tdifx(:,1);
            sxv = tstd(:,1);
            txv = 'Temperature';
            t2xv = 'T';t3xv = 't';
        elseif kk == 3
            xxv = pdifx(:,1);
            sxv = pstd(:,1);
            txv = 'Pressure';
            t2xv = 'P';t3xv = 'p';
        end
        xxmean = nanmean(xxv);
        xxstd = nanstd(xxv);
        xxerr = median(sxv)/sqrt(nvec-1);

        % Remove outliers do again if one makes a big difference to variance
        ioutlie = 1; n = 0; itok = false(size(xxv));
        while ioutlie > 0 && ~isnan(xxmean+xxstd) && n<100
            yymean = xxmean;
            yystd = xxstd;
            itok = abs(xxv-xxmean) < 2.32*xxstd;
            xxmean = mean(xxv(itok));
            xxstd = std(xxv(itok));
            if xxstd > 0.7*yystd
                ioutlie = -1;
            end
            n = n+1;
        end

        % Repeat process for std
        sxmean = nanmean(sxv);
        sxstd = nanstd(sxv);
        % Remove outliers do again if one makes a big difference to variance
        ioutlie = 1; n = 0; itok2 = false(size(xxv)); itok3 = itok2;
        while ioutlie > 0 && ~isnan(sxmean+sxstd) && n<100
            symean = sxmean;
            systd = sxstd;
            itok2 = abs(sxv-sxmean) < 2.32*sxstd;
            itok3 = sxv-sxmean < 2.32*sxstd;   % Because not a problem if is lower than usual
            sxmean = mean(sxv(itok2));
            sxstd = std(sxv(itok2));
            if sxstd > 0.7*systd
                ioutlie = -1;
            end
            n = n+1;
        end

        % Now write restults to a string
        for i = 1:nvec
            if ~itok(i)
                toutlie(i) = {[char(toutlie(i)) t2xv]};
            end
            if ~itok3(i)
                toutlie(i) = {[char(toutlie(i)) t3xv]};
            end
        end

        fprintf(1,'\n%s \n',txv)
        fprintf(1,'Mean of all differences %8.5f std %8.5f Std err %8.5f No of outliers %i \n', ...
            xxmean,xxstd,xxerr,nvec-sum(itok))
        if sum(~itok >= 1)
            fprintf(1,' %i  ',vec(~itok))
        end
        fprintf(1,' \n')
    end
end
% Display CTD characteristics


fprintf(ilogf,'\n%s: \n',lct);
fprintf(ilogf,'CTD mean temp (1,2) %8.4f %8.4f \n',ctd1_temp_mn,ctd2_temp_mn);
fprintf(ilogf,'CTD sdev temp (1,2) %8.5f %8.5f \n',ctd1_temp_st,ctd2_temp_st);
fprintf(ilogf,'CTD sensor mean and sdev diff tenperature  %8.5f %8.5f  \n', ...
    t_ctd_m,t_ctd_s);
fprintf(ilogf,'CTD mean cond (1,2) %8.4f %8.4f \n',ctd1_cond_mn,ctd2_cond_mn);
fprintf(ilogf,'CTD sdev cond (1,2) %8.5f %8.5f \n',ctd1_cond_st,ctd2_cond_st);
fprintf(ilogf,'CTD sensor mean and sdev diff conductivity %8.5f %8.5f  \n', ...
    c_ctd_m,c_ctd_s);
fprintf(ilogf,'CTD mean oxy (1,2) %8.4f %8.4f \n',ctd1_oxy_mn,ctd2_oxy_mn);
fprintf(ilogf,'CTD sdev oxy (1,2) %8.5f %8.5f \n',ctd1_oxy_st,ctd2_oxy_st);
fprintf(ilogf,'CTD sensor mean and sdev diff oxygen %8.5f %8.5f  \n ', ...
    o_ctd_m,o_ctd_s);
fprintf(ilogf,'\n%s\n',ptittxt);
fprintf(ilogf,'\n Number of MicroCATs = %i \n ',nvec);
fprintf(ilogf,'Stats for cast %s dday %7.3f to %7.3f with CTD sensor set %s \n ', ...
    cast,tm1p,tm2p,ctdsen);
fprintf(ilogf,'------------------------------------ \n');
if pwarn
    fprintf(ilogf,'At maximum overlap depth of CTD cast and MicroCATdata press = %6.1f db std = %4.1f \n',meanctdpr,stdctdpr );
else
    fprintf(ilogf,'At maximum depth of CTD cast press = %6.1f db std = %4.1f \n',meanctdpr,stdctdpr );
end
fprintf(ilogf,'Serial  No.       Conductivity     Temperature     Pressure       Oxygen    Outliers \n');
fprintf(ilogf,'Number Samples  Mean Dif  St.D.  Mean Dif  St.D.  MeanDif  St.D. MeanDif  St.D. \n');
jj = 1;
for i = 1:length(vec)
    fprintf(ilogf,' %5i %5i  %7.4f %8.5f  %7.4f %8.5f  %4.1f %4.1f   %4.1f %4.1f %s \n', ...
        vec(i),nimpt(i,jj),cdifx(i,jj),cstd(i,jj),tdifx(i,jj),tstd(i,jj), ...
        pdifx(i,jj),pstd(i,jj), odifx(i,jj), ostd(i,jj)  ,char(toutlie(i)));
end
fprintf(ilogf,'Outliers only calculated for temp, cond and press \n');
jj = 2;
fprintf(ilogf,'------------------------------------ \n');
fprintf(ilogf,'At nominal depth of instrument \n');
fprintf(ilogf,'Serial  No.   Nom.  Depth of   Conductivity     Temperature     Pressure      Oxygen\n');
fprintf(ilogf,'Number Samples Depth  Test    Dif  St.D.       Dif  St.D.    Dif  St.D.    Dif  St.D.\n');
for i = 1:length(vec)
    fprintf(ilogf,' %5.0f %4.0f %5.0f %5.0f  %7.4f  %8.5f  %7.4f  %8.5f  %4.1f  %4.1f  %4.1f  %4.1f\n', ...
        vec(i),nimpt(i,jj),zmic(i),ptest(i),cdifx(i,jj),cstd(i,jj),tdifx(i,jj),tstd(i,jj), ...
        pdifx(i,jj),pstd(i,jj), odifx(i,jj), ostd(i,jj));
end
%%
doplot=1;
if doplot
    figure(1);clf;grid on;grid minor
    subplot(4,1,1);
    bar([cdifx(:,1), 15*cstd(:,1)])
    ylim([-.025 0.025])

    set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx)
    title([ptittxt ' (b: Cond diff; r: Cond StD X 15)'])
    subplot(4,1,2);grid on;grid minor
    bar([tdifx(:,1), 3*tstd(:,1)])
    ylim([-0.005 0.005])

    set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx)
    title('b: Temp diff; r: Temp StD X 3')
    subplot(4,1,3);grid on;grid minor
    bar([pdifx(:,1), 20*pstd(:,1)])
    ylim([-40 40])

    set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx)
    title('b:Pres diff; r: Pres StD X 20')
    subplot(4,1,4);grid on;grid minor
    bar([pdifx(:,2), 10*pstd(:,2)])
    ylim([-25 25])

    set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx)
    title('At nominal pressure (b:Pres diff; r: Pres StD X 10)')


    % Finally save plotfile
    set(gcf,'PaperUnits','centimeters','PaperPosition',[-2 0 27 18 ])
    print('-depsc', pd.stage2fig)
    % Display the results
end


%%
figure(2);clf
subplot(411);hold on;grid on;grid minor;
errorbar(cdifx(:,1),cstd(:,1))
ylim([-.03 0.03])
%ylim([min(cdifx(:,1)-max(cstd(:,1))) max(cdifx(:,1)+max(cstd(:,1)))])
set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx,'XTickLabelRotation',45)

title([{ptittxt,'cond diff @ pmax'}])
ylabel('[mS/cm]');

subplot(4,1,2);hold on;grid on;grid minor;
errorbar(tdifx(:,1),tstd(:,1))
%ylim([-0.004 0.001])
ylim([min(tdifx(:,1)-max(tstd(:,1))) max(tdifx(:,1)+max(tstd(:,1)))])

set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx,'XTickLabelRotation',45)
title('Temp diff @ pmax')
ylabel('[^oC]')

subplot(4,1,3);hold on;grid on;grid minor;
errorbar(pdifx(:,1),pstd(:,1))
ylim([-5 5])
%ylim([min(pdifx(:,1)-max(pstd(:,1))) max(pdifx(:,1)+max(pstd(:,1)))])

set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx,'XTickLabelRotation',45)
title('Press diff @ pmax')
ylabel('[dbar]')

subplot(4,1,4);hold on;grid on;grid minor;
errorbar(pdifx(:,2),pstd(:,2))
%ylim([-10 10])
yl = [min(pdifx(:,2)-max(pstd(:,2))) max(pdifx(:,2)+max(pstd(:,2)))];
if isnan(sum(yl))
    yl = [-10 10];
end
ylim(yl)

set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx,'XTickLabelRotation',45)
title({'Pdiff @ nominal deployment pressure'})
ylabel('[dbar]')

% Finally save plotfile
set(gcf,'PaperUnits','centimeters','PaperPosition',[-2 0 27 18 ])
print('-dpng', pd.stage2fig)

type(pd.stage2log)