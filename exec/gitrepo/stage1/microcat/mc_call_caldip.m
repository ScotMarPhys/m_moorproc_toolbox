% MC_CALL_CALDIP_V4 is a script that performs stage1 processing
% on microcat data from CTD calibration casts (caldips).  It
% converts microcat data from raw to rodb format for an entire
% caldip, and plots it with CTD data.
%
% This was created from mc_call_caldip_jc103_v4.m so we don't keep creating
% cruise specific ones based on name - they are archived by cruise anyway.
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
% 06/05/2017 LOH: Adapted for OSNAP, (dy078 cruise)
%
%
% V3 on D382 by DAS added lines so that when data read in any records outisde
% start and end times in the info.dat file are ignored.
% DR added  sensor compatibility on JC103 prior to GDM creating V4
% V4 on JC103 by GDM more efficient loading to stop flashing figures;
% DR removed cruise name from function naming so just update this file each
% time
%
% modified on AR30 July 2018 by LOH to display 02 data from SBE37-ODO at bottle stops
% modified on AR30 July 2018 by LOH to load ctd data in native cnv format without header , 
% dotctd has to be set to 99 (1 if mstar format)
% modified on AR30 July 2018 by LOH: Add if/else conditions in the script
% to be able to process specific cases when microcat starting jday is >365

% modified Clare Johnson, Nov 2018 to include basic calibration for CTD
% oxygen sensor (comparison with bottle stops only)

% Modified for dy120 by S Jones, Oct 2020, using mc_call_ar30_cj.m as
% reference (see above)

% under verison control (GitHub) so no longer need for comments. ]
% L. Drysdale 2020

% -----------------------------------------------------------------
% --- This is the information that needs to be modified for -------
% --- different users, directory trees, and moorings --------------
% ----------------------------------------------------------------

clearvars -except MEXEC MEXEC_A MEXEC_G moordatadir ctddir YEAR this_cruise cruise_ctd;
close all;

do_microcat2rodb = 1; % if 0, Does not write rodb format files. Useful for fast rerun of plots

cruise = this_cruise;

cast = '36'; 

doctd = 1;% 1; % whether to load and plot CTD data: 1 if mstar format, 99 if native cnv file (without header)
jd0 = julian(YEAR,1,0,0); % set to current year

ctdnum = sprintf('%03d',str2num(cast));

% --- set paths for data input and output ---
 
inpath    = [moordatadir '/raw/' cruise '/microcat_cal_dip/cast',cast,'/'];
outpath   = [moordatadir '/proc_calib/' cruise '/cal_dip/microcat/cast' cast '/'];
infofile  = [moordatadir '/proc_calib/' cruise '/cal_dip/cast',cast,'info.dat'];
ctdinfile = [ctddir '/ctd_' cruise '_' ctdnum '_psal.nc'];

%addpath('/noc/users/pstar/di359/data/mexec_processing_scripts/',path);
%jd0 = julian(2014,1,0,0);
if doctd == 1
    [d, h]=mload(ctdinfile,'/');
    HH = h.data_time_origin(4)+h.data_time_origin(5)/60+h.data_time_origin(6)/3600;
    jtime=h.data_time_origin(1:3);
    jtime=[jtime HH];
    d.timeJ=julian(jtime)+d.time/86400-jd0;
    if ~isfield(d,'oxygen2')
        d.oxygen2 = d.oxygen1;
    end
elseif doctd == 99 % using cnv file instead of .nc
    warning('hardwired cast numbers/times')
    if contains(cast,'1')
       start_date_cast = datenum(2018,07,02,18,17,06);  
    elseif contains(cast,'2')
       start_date_cast = datenum(2018,07,03,17,17,42);  
    elseif contains(cast,'3')
       start_date_cast = datenum(2018,07,04,17,19,06);         
    elseif contains(cast,'5')
       start_date_cast = datenum(2018,07,07,03,08,39); 
       start_date_cast_9141 = datenum(2016,7,4,16,00,01);
    elseif contains(cast,'6')
       start_date_cast = datenum(2018,07,07,13,13,15);
    end
    read_ctd_cnv(ctdinfile,start_date_cast)
end

% --- get mooring information from infofile ---
[zp,id,sn]= rodbload(infofile,'z:instrument:serialnumber');

% --  Start and end times
[stdt,stti,endt,enti]= rodbload(infofile,'StartDate:StartTime:EndDate:EndTime');
jdstt = julian([stdt;stti(1)+stti(2)/60]')-jd0;
jdend = julian([endt;enti(1)+enti(2)/60]')-jd0;

% --- vector of serial numbers ---
ii = find(id >= 332 & id <= 337);
invalid_sn = sn(setdiff(1:length(sn),ii));
%ylf dy146 replace vec with sn as they are identical
id2 = id(ii);
zp = zp(ii);
sn = sn(ii);

% --- create log file ---
if exist(outpath,'dir')==0
    mkdir(outpath)
end
fidlog = fopen([outpath,'microcat2rodb.log'],'w');

legend_handle = ['[';'[';'['];
legend_string = [];

mx_length = 0;
valid_sn = false(1,length(sn));
% --- read data loop --
for i = 1:length(sn)
    % display( [,num2str(sn(i)),])
    
    fprintf(fidlog,'\n\n');
    
    %ylf dy146 condensed this part
    % try to find infile in list of possibilities
    infiles = {[inpath,sprintf('%4.4d',sn(i)),'cal2.asc'];
        [inpath,sprintf('%4.4d',sn(i)),'cal.asc'];
        [inpath,sprintf('%3.3d',sn(i)),'cal.asc'];
        [inpath,sprintf('%4.4d',sn(i)),'CAL.asc'];
        [inpath,sprintf('%3.3d',sn(i)),'CAL.asc'];
        [inpath,'cal',sprintf('%4.4d',sn(i)),'.asc'];
        [inpath,sprintf('%4.4d',sn(i)),'_cal_dip2.asc'];
        [inpath,sprintf('%4.4d',sn(i)),'_cal_dip_data2.asc'];
        [inpath,sprintf('%4.4d',sn(i)),'_test.asc'];
        [inpath,sprintf('%4.4d',sn(i)),'_data.asc'];
        [inpath,sprintf('%4.4d',sn(i)),'_cal_dip.asc'];
        [inpath,sprintf('%4.4d',sn(i)),'_cal_dip_data.asc'];
        [inpath,sprintf('%4.4d',sn(i)),'_cal_dip_data.cnv'];
        [inpath,sprintf('%4.4d',sn(i)),'_Cal_Dip_Data.cnv']};
    while ~isempty(infiles)
        infile = infiles{1};
        if exist(infile,'file')==2
            break
        else
            infiles(1) = [];
        end
    end
    outfile = [outpath, 'cast', cast ,'_',sprintf('%4.4d',sn(i)),'.raw'];

    % --- convert from raw to rodb format ---
    
    if exist(infile,'file')
        valid_sn(i)=true;
        %------------------------------------------------------------
        % specific cases when microcat starting jday is >365
        if(contains(cruise,'ar30') & sn(i) == 11327 & contains(cast,'6'))
            dateoffsetmc = 2017;
        elseif  (contains(cruise,'ar30') & sn(i) == 9141 & contains(cast,'5'))
             dateoffsetmc = 2016;          
        elseif sn(i)==6322
            disp('ADDING AN OFFSET')
            dateoffsetmc=60/86400;
        else
            dateoffsetmc = 0;
        end
        %--------------------------------------------------------------
        if do_microcat2rodb
            microcat2rodb(infile,outfile,infofile,fidlog,'y',dateoffsetmc)
        else
            disp(['*** WARNING **** NOT WRITING RODB FILES']);
        end

        % --- load rodb data ---
        %ylf dy146 split conditional part into two to not repeat common
        %part
        if id2(i)~=335 % checks if not ODO microcat
            [yy,mm,dd,hh,c,t,p] = rodbload(outfile,'yy:mm:dd:hh:c:t:p');
        else % treats as ODO
            [yy,mm,dd,hh,c,t,p,ot,o2] = rodbload(outfile,'yy:mm:dd:hh:c:t:p:ot:o2');
        end
        jd = julian(yy,mm,dd,hh)-jd0;
        %ylf dy146 files can be different sizes anyway, so no need to keep
        %all the NaNs we would put at the ends (~ixjd), right?
        ixjd = jd >jdstt & jd < jdend;
        mx_length = max(mx_length,length(ixjd));
        data(i).jd = jd(ixjd);
        data(i).c = c(ixjd);
        data(i).t = t(ixjd);
        data(i).p = p(ixjd);
        if id2(i)==335
            data(i).ot = ot(ixjd);
            data(i).os = o2(ixjd);
        end
    end

end
invalid_sn=[invalid_sn; sn(~valid_sn)]; %combine with the ones we saved earlier that don't have MC or ODO id
sn=sn(valid_sn);
zp = zp(valid_sn);
data = data(valid_sn);
snl = [num2str(sn) repmat(' (',length(sn),1) num2str(zp) repmat(')',length(sn),1)];

jd = ones(length(sn),mx_length)+nan;
c = ones(length(sn),mx_length)+nan;
t = ones(length(sn),mx_length)+nan;
p = ones(length(sn),mx_length)+nan;
ot = ones(length(sn),mx_length)+nan;
o2 = ones(length(sn),mx_length)+nan;

ctdc = ones(length(sn),mx_length)+nan;
ctdt = ones(length(sn),mx_length)+nan;
ctdc2 = ones(length(sn),mx_length)+nan;
ctdt2 = ones(length(sn),mx_length)+nan;
ctdp = ones(length(sn),mx_length)+nan;
ctdot = ones(length(sn),mx_length)+nan;
ctdo2 = ones(length(sn),mx_length)+nan;

for i = 1:length(sn)
    ll = length(data(i).jd);
    jd(i,1:ll) = data(i).jd;
    t(i,1:ll) = data(i).t;
    c(i,1:ll) = data(i).c;
    p(i,1:ll) = data(i).p;
    if isfield(data(i),'o2')
        o2(i,1:ll) = data(i).o2;
        ot(i,1:ll) = data(i).ot;
    end
    % interp here if doctd
    if doctd==1 || doctd==99
        ni = ~isnan(jd(i,:));
        ctdc2(i,ni)=interp1(d.timeJ, d.cond, jd(i,ni));
        ctdt2(i,ni)=interp1(d.timeJ, d.temp, jd(i,ni));
        ctdp(i,ni)=interp1(d.timeJ, d.press, jd(i,ni));
        ctdot(i,ni)=interp1(d.timeJ, d.temp, jd(i,ni));
        ctdo2(i,ni)=interp1(d.timeJ, d.oxygen, jd(i,ni));
    end
end

pdiff = ['mean p diff = ' num2str(nanmedian(p'-ctdp'))];
cdiff = ['mean c diff = ' num2str(nanmedian(c'-ctdc2'))];
tdiff = ['mean t diff = ' num2str(nanmedian(t'-ctdt2'))];
if sum(id2==335)>0
    odiff = ['mean o2 diff = ' num2str(nanmedian(o2'-ctdo2'))];
end

disp(['instrument SN = ' num2str(sn')]);
disp(pdiff)
disp(cdiff)
disp(tdiff)
if sum(id2==335)>0
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
a=size(jd);
colours = repmat([0 0 .8; 0 1 0; 1 0 0; 0 1 1; 1 .2 .9; .7 .5 0; 0 .4 0; 0 0 0],4,1);
linestyles={'-','-','-','-','-','-','-','-','-.','-.','-.','-.','-.','-.','-.','-.','--','--','--','--','--','--','--','--',':',':',':',':',':',':',':',':'};
for i=1:a(1)
    plot(jd(i,:),c(i,:),'color',colours(i,:),'linestyle',linestyles{i});
end
legend(snl,'location','eastoutside')
ylabel('conductivity')
xlabel(['yearday (relative to ' num2str(YEAR) '/1/0 00:00)'])
title(['CAST ' cast ' Calibration Dip'])
if doctd
    hl = plot(d.timeJ,d.cond1-.02,d.timeJ,d.cond1+.02,d.timeJ,d.cond1,d.timeJ,d.cond2); 
    set(hl(1:2),'color',[.8 .8 .8],'linestyle','--');
    set(hl(3),'color',[0 0 0]);set(hl(4),'color',[.4 .4 .4]); % ylf jc145 plot both
    title(['CAST ' cast  ' Calibration Dip (-k=CTD1,gray=CTD2)'])
end
xlim(sort([jdstt-.01 jdend+.01])); grid
orient tall
print(gcf,'-depsc',[outfig '_cond.ps'])
saveas(gcf,[outfig '_cond.fig'],'fig')

figure(35);hold on; box on;
for i=1:a(1)
    plot(jd(i,:),t(i,:),'color',colours(i,:),'linestyle',linestyles{i});
end
legend(snl,'location','eastoutside')
ylabel('temperature')
xlabel(['yearday (relative to ' num2str(YEAR) '/1/0 00:00)'])
title(['CAST ' cast ' Calibration Dip'])
if doctd
    hl = plot(d.timeJ,d.temp1+.005,d.timeJ,d.temp1-.005,d.timeJ,d.temp1,d.timeJ,d.temp2); 
    set(hl(1:2),'color',[.8 .8 .8],'linestyle','--')
    set(hl(3),'color',[0 0 0]);set(hl(4),'color',[.4 .4 .4]); % ylf jc145 plot both
    title(['CAST ' cast  ' Calibration Dip (-k=CTD1,gray=CTD2)'])
end
jdstt
jdend
xlim([jdstt-0.01 jdend+0.01]); grid
orient tall
figname=[outfig,'_temp'];
print(gcf,'-depsc',[figname,'.ps'])
figname=[figname,'.fig'];
saveas(gcf,figname,'fig')
%%
figure(36);hold on; box on;
for i=1:a(1)
    plot(jd(i,:),p(i,:),'color',colours(i,:),'linestyle',linestyles{i});
end
legend(snl,'location','eastoutside')
ylabel('pressure')
xlabel(['yearday (relative to ' num2str(YEAR) '/1/0 00:00)'])
title(['CAST ' cast ' Calibration Dip'])
if doctd
    hl = plot(d.timeJ,d.press+5,d.timeJ,d.press-5,d.timeJ,d.press,'k-'); % bim
    set(hl(1:2),'color',[.8 .8 .8],'linestyle','--')
    set(hl(3),'color',[0 0 0])
    title(['CAST ' cast ' Calibration Dip (-k=CTD)'])
end
xlim([jdstt-0.01 jdend+0.01]); grid
orient tall
print(gcf,'-depsc',[outfig '_pres.ps'])
saveas(gcf,[outfig '_pres.fig'],'fig')

if find(id2==335)

    figure(38);hold on; box on;
    for i=1:a(1)
        plot(jd(i,:),o2(i,:),'color',colours(i,:),'linestyle',linestyles{i});
    end
    legend(snl,'location','eastoutside')
    ylabel('oxygen')
    xlabel(['yearday (relative to ' num2str(YEAR) '/1/0 00:00)'])
    title(['CAST ' cast ' Calibration Dip'])
    if doctd
        hl = plot(d.timeJ,d.oxygen1-20,d.timeJ,d.oxygen1+20,d.timeJ,d.oxygen1,d.timeJ,d.oxygen2); 
        set(hl(1:2),'color',[.8 .8 .8],'linestyle','--')
        set(hl(3),'color',[0 0 0]); set(hl(4),'color',[.4 .4 .4]); % ylf jc145 plot both
        title(['CAST ' cast ' Calibration Dip (-k=CTD)'])
    end
    xlim([jdstt-0.01 jdend+0.01]); grid
    orient tall
    print(gcf,'-depsc',[outfig '_oxy.ps'])
    saveas(gcf,[outfig '_oxy.fig'],'fig')
end

disp(['number of MicroCATs processed = ' num2str(length(sn))])
if ~isempty(invalid_sn)
    disp('WARNING: Some MicroCATs listed in info.dat file were not processed')
    disp('Serial numbers:')
    sprintf('%d\n',invalid_sn)
end
