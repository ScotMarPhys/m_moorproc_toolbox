% MC_CALL_CALDIP_JC103 is a script that performs stage1 processing
% on microcat data from CTD calibration casts (caldips).  It
% converts microcat data from raw to rodb format for an entire
% caldip, and plots it with CTD data.
%
% It calls microcat2rodb_3 (to convert microcat_data), rodbload.m,
% timeaxis.m, auto_filt.m, julian.m
%
% quick-look for microcat calibration profiles and conversion to rodb
%
% 10/07 scu, D324: now reads and plots ctd ctu file using pload
% 28/10/08 ZBS: modified to print combined graphs for d334
% 22/10/09 PGW: updated for D344
% 26/03/10 SCU: Added mload to import netcdf CTD 1hz file
% 27/03/10 scu/zbs, oc459: updated for oc459, copied this file from
%   mc_call_calib2_oc459.m, added a descriptive header
% 19/12/10 efw, d359: updated for d359, copied this file from
%   mc_call_caldip_oc459.m.
% V2 February 2012, GDM: updated for RB1201,
% - now has accurate legends irregardless of discrepancies in number of
% mcat files and no of files info.dat is expecting
% - including changing parse_cnv for the files on the Ron Brown
%
% V3 on D382 by DAS added lines so that when data read in any records outisde
% start and end times in the info.dat file are ignored.
% DR added oxygen sensor compatibility on JC103 prior to GDM creating V4
% V4 on JC103 by GDM more efficient loading to stop flashing figures;
% -----------------------------------------------------------------
% --- This is the information that needs to be modified for -------
% --- different users, directory trees, and moorings --------------
% ----------------------------------------------------------------

clearvars -except pathosnap 
close all;

if exist('mload')~=2 % if mstart routines not already in the path
    pathmsetup = [pathosnap '/../../Cruises/PE400_2015/mstar/pe400/data/mexec_processing_scripts_v2/m_setup_pe400_mdrive.m'];
    run(pathmsetup)
end

%cruise = 'dy053';cruise_ctd = 'dy053'; cast = '22'; dateoffsetmc = 2016; % origine of the julian day in the matlab cnv files
cruise = 'pe400';cruise_ctd = 'pe400'; cast = '26'; dateoffsetmc = 2015; 

doctd = 0; % whether to load and plot CTD data
jd0 = julian(dateoffsetmc,1,0,0); % set to current year


ctdnum = sprintf('%03d',str2num(cast));

basedir = [pathosnap filesep 'data' filesep];
ctddir = [pathosnap filesep 'data' filesep cruise_ctd filesep]; 


% -----------------------------------------------------------------


% --- set paths for data input and output ---

inpath    = [basedir 'moor/raw/' cruise '/microcat_cal_dip/cast',cast,'/'];
outpath   = [basedir 'moor/proc_calib/' cruise '/cal_dip/microcat/cast' cast '/'];
infofile  = [basedir 'moor/proc_calib/' cruise '/cal_dip/cast',cast,'info.dat'];
ctdinfile = [ctddir  'ctd_' cruise_ctd '_' ctdnum '_psal.nc'];


if doctd == 1;
    [d h]=mload(ctdinfile,'/');
    
    HH = h.data_time_origin(4)+h.data_time_origin(5)/60+h.data_time_origin(6)/3600;
    jtime=h.data_time_origin(1:3);
    jtime=[jtime HH];
    d.timeJ=julian(jtime)+d.time/86400-jd0;
    
end

% --- get mooring information from infofile ---
[id,sn]= rodbload(infofile,'instrument:serialnumber');

% --  Start and end times
[stdt,stti,endt,enti]= rodbload(infofile,'StartDate:StartTime:EndDate:EndTime');
jdstt = julian([stdt;stti(1)+stti(2)/60]')-jd0;
jdend = julian([endt;enti(1)+enti(2)/60]')-jd0;

% --- vector of serial numbers ---
ii = find(id >= 332 & id <= 337);
vec = sn(ii);
id2=id(ii);


% --- create log file ---
fidlog = fopen([outpath,'microcat2rodb.log'],'w');

legend_handle = ['[';'[';'['];
legend_string = [];

mx_length = 0;
mx_length_o2 = 0;
% --- read data loop --
for i = 1:length(vec)
    % display( [,num2str(vec(i)),])
    fprintf(fidlog,'\n\n');
    
    infile = [inpath,sprintf('%4.4d',vec(i)),'cal2.asc'];
    
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%4.4d',vec(i)),'cal.asc'];
    end
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%3.3d',vec(i)),'cal.asc'];
    end
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%4.4d',vec(i)),'CAL.asc'];
    end
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%3.3d',vec(i)),'CAL.asc'];
    end
    if exist(infile) ~= 2
        infile = [inpath,'cal',sprintf('%4.4d',vec(i)),'.asc'];
    end
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%4.4d',vec(i)),'_cal_dip2.asc'];
    end
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%4.4d',vec(i)),'_cal_dip_data2.asc'];
    end
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%4.4d',vec(i)),'_test.asc'];
    end
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%4.4d',vec(i)),'_cal_dip.asc'];
    end
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%4.4d',vec(i)),'_cal_dip_data.cnv'];
    end
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%4.4d',vec(i)),'_cal_dip_data.asc'];
    end    
    
    outfile = [outpath, 'cast', cast ,'_',sprintf('%4.4d',vec(i)),'.raw'];
    
    % --- convert from raw to rodb format ---
    
    microcat2rodb_4(infile,outfile,infofile,fidlog,'w',dateoffsetmc)
    
    % --- load rodb data ---
    if id2(i)~=335 % checks if not ODO microcat
        [yy,mm,dd,hh,c,t,p] = rodbload(outfile,'yy:mm:dd:hh:c:t:p');
        jd = julian(yy,mm,dd,hh)-jd0;
        ixjd = jd >jdstt & jd < jdend;
        c(~ixjd) = NaN;
        t(~ixjd) = NaN;
        p(~ixjd) = NaN;
        eval(['jd' sprintf('%d',sn(i)) '= jd;']);
        eval(['c' sprintf('%d',sn(i)) '= c;']);
        eval(['t' sprintf('%d',sn(i)) '= t;']);
        eval(['p' sprintf('%d',sn(i)) '= p;']);
        if length(jd)>mx_length;mx_length=length(jd);end;
    else % treats as ODO
        [yy,mm,dd,hh,c,t,p,ot,o2] = rodbload(outfile,'yy:mm:dd:hh:c:t:p:ot:o2');
        jd = julian(yy,mm,dd,hh)-jd0;
        ixjd = jd >jdstt & jd < jdend;
        c(~ixjd) = NaN;
        t(~ixjd) = NaN;
        p(~ixjd) = NaN;
        o2(~ixjd) = NaN;
        ot(~ixjd) = NaN;
        eval(['jd' sprintf('%d',sn(i)) '= jd;']);
        eval(['c' sprintf('%d',sn(i)) '= c;']);
        eval(['t' sprintf('%d',sn(i)) '= t;']);
        eval(['p' sprintf('%d',sn(i)) '= p;']);
        eval(['ot' sprintf('%d',sn(i)) '= ot;']);
        eval(['o2' sprintf('%d',sn(i)) '= o2;']);
        
        if length(jd)>mx_length_o2;mx_length_o2=length(jd);end;
        
    end
    
end
if mx_length_o2>mx_length; mx_length = mx_length_o2; end;

jd = ones(length(sn),mx_length)+nan;
c = ones(length(sn),mx_length)+nan;
t = ones(length(sn),mx_length)+nan;
p = ones(length(sn),mx_length)+nan;
ot = ones(length(sn),mx_length)+nan;
o2 = ones(length(sn),mx_length)+nan;

ctdc = ones(length(sn),mx_length)+nan;
ctdt = ones(length(sn),mx_length)+nan;
ctdp = ones(length(sn),mx_length)+nan;
ctdot = ones(length(sn),mx_length)+nan;
ctdo2 = ones(length(sn),mx_length)+nan;

for i = 1:length(vec)
    ll = eval(['length(jd' sprintf('%d',sn(i)) ')']);
    jd(i,1:ll)=eval(['jd' sprintf('%d',sn(i)) ';']);
    t(i,1:ll)=eval(['t' sprintf('%d',sn(i)) ';']);
    c(i,1:ll)=eval(['c' sprintf('%d',sn(i)) ';']);
    p(i,1:ll)=eval(['p' sprintf('%d',sn(i)) ';']);
    if exist(['o2' sprintf('%d',sn(i))])
        ot(i,1:ll)=eval(['ot' sprintf('%d',sn(i)) ';']);
        o2(i,1:ll)=eval(['o2' sprintf('%d',sn(i)) ';']);
    end
    % interp here if doctd
    if doctd
        ni = ~isnan(c(i,:)); ctdc(i,ni)=interp1(d.timeJ, d.cond, jd(i,ni));
        ni = ~isnan(t(i,:)); ctdt(i,ni)=interp1(d.timeJ, d.temp, jd(i,ni));
        ni = ~isnan(p(i,:)); ctdp(i,ni)=interp1(d.timeJ, d.press, jd(i,ni));
        ni = ~isnan(ot(i,:)); ctdot(i,ni)=interp1(d.timeJ, d.temp, jd(i,ni));
    end
end

pdiff = ['mean p diff = ' num2str(nanmedian(p'-ctdp')) ...
    ' --- ' num2str(vec(i))];
cdiff = ['mean c diff = ' num2str(nanmedian(c'-ctdc')) ...
    ' --- ' num2str(vec(i))];
tdiff = ['mean t diff = ' num2str(nanmedian(t'-ctdt')) ...
    ' --- ' num2str(vec(i))];
if id2(i)==335
    odiff = ['mean o2 diff = ' num2str(nanmedian(o2'-ctdo2')) ...
        ' --- ' num2str(vec(i))];
end

disp(pdiff)
disp(cdiff)
disp(tdiff)
if id2(i)==335
    disp(odiff)
end

fprintf(fidlog,'%s \n','for mcats wi serial nos ')
fprintf(fidlog,'%s \n',num2str(sn));
fprintf(fidlog,'%s \n','mean p diff = ')
fprintf(fidlog,'%s \n',pdiff);
fprintf(fidlog,'%s \n','mean c diff = ')
fprintf(fidlog,'%s \n',cdiff);
fprintf(fidlog,'%s \n','mean t diff = ')
fprintf(fidlog,'%s \n',tdiff);
if id2(i)==335
    fprintf(fidlog,'%s \n','mean o2 diff = ')
    fprintf(fidlog,'%s \n',odiff);
end

outfig = [outpath, 'cast', cast ,'_all'];
figure(34); hold on; box on;
% line below is not working correctly so setting up a loop for plotting instead
%set(gcf,'DefaultAxesLineStyleOrder','-|-|-|-|-|-|-|-.|-.|-.|-.|-.|-.|-.|--|--|--|--|--|--|--|:|:|:|:|:|:|:')
a=size(jd);
colours='bgrcmykbgrcmykbgrcmykbgrcmyk';
linestyles={'-','-','-','-','-','-','-','-.','-.','-.','-.','-.','-.','-.','--','--','--','--','--','--','--',':',':',':',':',':',':',':'};
for i=1:a(1)
    plot(jd(i,:),c(i,:),'color',colours(i),'linestyle',linestyles{i});
end
legend(num2str(sn),'location','eastoutside')
ylabel('conductivity')
xlabel('yearday (relative to 2014/1/0 00:00)')
title(['CAST ' cast ' Calibration Dip'])
if doctd
    plot(d.timeJ,d.cond,'k-') % bim
    title(['CAST ' cast  ' Calibration Dip (-k=CTD)'])
end
xlim([jdstt-0.01 jdend+0.01]);
orient tall
print(gcf,'-depsc',[outfig '_cond.ps'])
saveas(gcf,[outfig '_cond.fig'],'fig')

figure(35);hold on; box on;
for i=1:a(1)
    plot(jd(i,:),t(i,:),'color',colours(i),'linestyle',linestyles{i});
end
legend(num2str(sn),'location','eastoutside')
ylabel('temperature')
xlabel('yearday (relative to 2012/1/0 00:00)')
title(['CAST ' cast ' Calibration Dip'])
if doctd
    plot(d.timeJ,d.temp,'color',[.3 .3 .3]); % testing for sensor drift
    title(['CAST ' cast ' Calibration Dip (-k=CTD)'])
end
xlim([jdstt-0.01 jdend+0.01]);
orient tall
print(gcf,'-depsc',[outfig '_temp.ps'])
saveas(gcf,[outfig '_temp.fig'],'fig')

figure(36);hold on; box on;
for i=1:a(1)
    plot(jd(i,:),p(i,:),'color',colours(i),'linestyle',linestyles{i});
end
legend(num2str(sn),'location','eastoutside')
ylabel('pressure')
xlabel('yearday (relative to 2012/1/0 00:00)')
title(['CAST ' cast ' Calibration Dip'])
if doctd
    plot(d.timeJ,d.press,'k-'); % bim
    title(['CAST ' cast ' Calibration Dip (-k=CTD)'])
end
xlim([jdstt-0.01 jdend+0.01]);
orient tall
print(gcf,'-depsc',[outfig '_pres.ps'])
saveas(gcf,[outfig '_pres.fig'],'fig')


disp(['number of MicroCATs processed = ' num2str(length(vec))])



