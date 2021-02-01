% mc_caldip_check_jc103 reads microcat caldip data from raw rodb file
% and compares data with the lowered CTD data
%
% DAS created this file based on mc_call_caldip

% LOH May 2017 - adapted for SBE37 ODO and CTD with two oxygen sensors


% -----------------------------------------------------------------
% --- This is the information that needs to be modified for -------
% --- different users, directory trees, and moorings --------------
% ----------------------------------------------------------------
clearvars -except pathosnap 
close all;



cruise_ctd = 'dy078'; % used for ctd data dir
cruise = cruise_ctd; % used for microcat data
cruise2=cruise_ctd; % used for ctd data name
%   Time origin for jday
jd0 = julian(2017,1,0,0);

basedir = [pathosnap filesep 'data' filesep];
ctddir = [pathosnap filesep 'data' filesep cruise_ctd filesep]; 

% -----------------------------------------------------------------
% User supplied information
% Select cast number and choose which set of instruments from CTD
cast = input('Which cast number? ','s');
ctdsen = input('Which CTD sensors (1 or 2?) ','s');
ctdnum = sprintf('%03d',str2num(cast));

% --- set paths for data input and output ---
outpath   = [basedir 'moor/proc_calib/' cruise '/cal_dip/microcat/cast' cast '/'];
infofile  = [basedir 'moor/proc_calib/' cruise '/cal_dip/cast',cast,'info.dat'];
%ctdinfile = [ctddir  'ctd_' cruise2 '_',ctdnum,'_raw.nc'];
ctdinfile = [ctddir  'ctd_' cruise2 '_',ctdnum,'_psal.nc'];


% ----------------- load CTD DATA   ----------------------------------

%  [d h]=mload(ctdinfile,'timeJ','prDM','t090C','c0S_slash_m',' ');
 [d h]=mload(ctdinfile,'time','press','temp1','cond1','temp2','cond2','oxygen1','oxygen2',' ','q');
 HH = h.data_time_origin(4)+h.data_time_origin(5)/60+h.data_time_origin(6)/3600;
 jtime=h.data_time_origin(1:3);
 jtime=[jtime HH];
 d.timeJ=julian(jtime)+d.time/86400-jd0;
 ctdtimetx = datestr(h.data_time_origin);
 ptittxt = sprintf('Cast %s start %s CTD sensor set %s',cast,ctdtimetx,ctdsen);

if strcmp(cruise,'d382')
    % Correction for Di382
    d.cond1=d.cond1*10;  
    d.cond2=d.cond2*10; 
end

% Selecttime period that we will analyse
mxpi = max(d.press);
imp = d.press > mxpi-10;
tm1p = min(d.timeJ(imp))+0.25*(max(d.timeJ(imp))- min(d.timeJ(imp)));
tm2p = max(d.timeJ(imp))-0.25*(max(d.timeJ(imp))- min(d.timeJ(imp)));
imp2 = d.timeJ > tm1p & d.timeJ < tm2p;
meanctdpr = nanmean(d.press(imp2));

% CTD stats during period
meanctdpr = nanmean(d.press(imp2));
stdctdpr = nanstd(d.press(imp2));
ctd1_cond_mn = nanmean(d.cond1(imp2));
ctd1_cond_st = nanstd(d.cond1(imp2));
ctd2_cond_mn = nanmean(d.cond2(imp2));
ctd2_cond_st = nanstd(d.cond2(imp2));
ctd1_temp_mn = nanmean(d.temp1(imp2));
ctd1_temp_st = nanstd(d.temp1(imp2));
ctd2_temp_mn = nanmean(d.temp2(imp2));
ctd2_temp_st = nanstd(d.temp2(imp2));
ctd1_oxy_mn = nanmean(d.oxygen1(imp2));
ctd1_oxy_st = nanstd(d.oxygen1(imp2));
ctd2_oxy_mn = nanmean(d.oxygen2(imp2));
ctd2_oxy_st = nanstd(d.oxygen2(imp2));

% --- get mooring information from infofile ---
'here'
infofile
[zins,id,sn]= rodbload(infofile,'z:instrument:serialnumber');
'done'

% --- vector of serial numbers ---
ii = find(id >= 332 & id <= 337);
vec = sn(ii);
id2=id(ii);
nvec = length(ii);
zmic = zins(ii);

% Open output file for text and set plot name
outlogf = [outpath,'microcat_check',cast,'.log'];
ilogf = fopen(outlogf,'w');
outplot = [outpath,'microcat_check',cast,'_plot'];

% --- read data loop --
for i = 1:nvec
% display( [,num2str(vec(i)),])
   outfile = [outpath, 'cast', cast ,'_',sprintf('%4.4d',vec(i)),'.raw'];
% --- load rodb data ---
   if id2(i)~=335 % checks if not ODO microcat
            [yy,mm,dd,hh,c,t,p] = rodbload(outfile,'yy:mm:dd:hh:c:t:p');
            o=c*nan;
   else
            [yy,mm,dd,hh,c,t,p,o] = rodbload(outfile,'yy:mm:dd:hh:c:t:p:o2');
   end
%  if (i > 6 & i<=12)  lstr='--'; elseif i>12  lstr ='-.'; else lstr = '-'; end 
   if (i > 7 & i<=14)  lstr='--'; elseif i>14  lstr ='-.'; else lstr = '-'; end 
%   disp(['Checking ',num2str(i),': s/n:',num2str(sn(i))])
%    pause
% Time variable
   jd = julian(yy,mm,dd,hh)-jd0;
% interpolate CTD onto microcat for a rough and ready mean diff
     if ctdsen == '1'
       pi = interp1(d.timeJ, d.press, jd);
       ti = interp1(d.timeJ, d.temp1, jd);
       ci = interp1(d.timeJ, d.cond1, jd);
       oi = interp1(d.timeJ, d.oxygen1, jd);       
     elseif ctdsen == '2'
       pi = interp1(d.timeJ, d.press, jd);
       ti = interp1(d.timeJ, d.temp2, jd);
       ci = interp1(d.timeJ, d.cond2, jd);
       oi = interp1(d.timeJ, d.oxygen2, jd);        
     end

% Select data from period when CTD kept at maximum depth
     impt = jd>tm1p & jd < tm2p;
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
       tix = interp1(d.timeJ, d.temp2-d.temp1, jd);
       cix = interp1(d.timeJ, d.cond2-d.cond1, jd);
       oix = interp1(d.timeJ, d.oxygen2-d.oxygen1, jd);      
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

'end'
% Quick look all data and search for outliers
toutlie(1:nvec) = {' '};
for kk = 1:3
    if kk == 1 
        xxv = cdifx(:,1); % microcat - ctd
        sxv = cstd(:,1);
        txv = 'Conductivity';
        t2xv = 'C';t3xv = 'c';
    elseif kk == 2 
        xxv = tdifx(:,1);
        sxv = tstd(:,1);
        txv = 'Temperateure';
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
ioutlie = 1;
while ioutlie > 0
    yymean = xxmean;
    yystd = xxstd;
    itok = abs(xxv-xxmean) < 2.32*xxstd;
    xxmean = mean(xxv(itok));
    xxstd = std(xxv(itok));
    if xxstd > 0.7*yystd
       ioutlie = -1;
    end
end

% Repeat process for std 
  sxmean = nanmean(sxv);
  sxstd = nanstd(sxv);
% Remove outliers do again if one makes a big difference to variance
ioutlie = 1;
while ioutlie > 0
    symean = sxmean;
    systd = sxstd;
    itok2 = abs(sxv-sxmean) < 2.32*sxstd;
    itok3 = sxv-sxmean < 2.32*sxstd;   % Because not a problem if is lower than usual
    sxmean = mean(sxv(itok2));
    sxstd = std(sxv(itok2));
    if sxstd > 0.7*systd
       ioutlie = -1;
    end
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
  fprintf(1,'Mean of all differences %8.5f std %8.5f Std err %8.5f No of outlliers %i \n', ...
    xxmean,xxstd,xxerr,nvec-sum(itok))
 if sum(~itok >= 1)
    fprintf(1,' %i  ',vec(~itok))
  end
  fprintf(1,' \n')
end

% Display CTD characteristics


fprintf(ilogf,'\nCTD mean temp (1,2) %8.4f %8.4f \n',ctd1_temp_mn,ctd2_temp_mn);
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
fprintf(ilogf,'Stats for cast %s jday %7.3f to %7.3f with CTD sensor set %s \n ', ...
    cast,tm1p,tm2p,ctdsen);
fprintf(ilogf,'------------------------------------ \n');
fprintf(ilogf,'At maximum depth of CTD cast press = %6.1f db std = %4.1f \n',meanctdpr,stdctdpr );
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

figure
subplot(4,1,1)
bar([cdifx(:,1), 15*cstd(:,1)])
ylim([-.025 0.025])
grid on;grid minor
set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx)
title([ptittxt ' (b: Cond diff; r: Cond StD X 15)'])
subplot(4,1,2)
bar([tdifx(:,1), 3*tstd(:,1)])
ylim([-0.005 0.005])
grid on;grid minor
set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx)
title('b: Temp diff; r: Temp StD X 3')
subplot(4,1,3)
bar([pdifx(:,1), 20*pstd(:,1)])
ylim([-40 40])
grid on;grid minor
set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx)
title('b:Pres diff; r: Pres StD X 20')
subplot(4,1,4)
bar([pdifx(:,2), 10*pstd(:,2)])
ylim([-25 25])
grid on;grid minor
set(gca,'xtick',1:1:length(vectx),'XTicklabel',vectx)
title('At nominal pressure (b:Pres diff; r: Pres StD X 10)')


% Finally save plotfile
set(gcf,'PaperUnits','centimeters','PaperPosition',[-2 0 27 18 ])
print('-depsc', outplot)
% Display the results
eval(['!cat ',outlogf]);

% --- That's it! ---
% keyboard















