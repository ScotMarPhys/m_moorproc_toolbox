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
% -----------------------------------------------------------------
% --- This is the information that needs to be modified for -------
% --- different users, directory trees, and moorings --------------
% ----------------------------------------------------------------

cruise = 'jc103';
cruise_ctd = 'jc103';
cast = '16';

ctdnum = sprintf('%03d',str2num(cast));
doctd = 1; % whether to load and plot CTD data
if exist('/Volumes/rpdmoc/rapid/data/exec/jc103/stage1/microcat/mc_call_caldip_jc103_v3.m','file')
    % using DR Mac with mount to banba on JC103
    basedir = '/Volumes/rpdmoc/rapid/data/';
    ctddir = '/Volumes/rpdmoc/jc103/data/ctd/';
    
else
    basedir = '/local/users/pstar/rpdmoc/rapid/data/';
    ctddir = '/local/users/pstar/cruise/data/ctd/';
end

year = 2014;

% -----------------------------------------------------------------


% --- set paths for data input and output ---

inpath    = [basedir 'moor/raw/' cruise '/microcat_cal_dip/cast',cast,'/'];
outpath   = [basedir 'moor/proc_calib/' cruise '/cal_dip/microcat/cast' cast '/'];
infofile  = [basedir 'moor/proc_calib/' cruise '/cal_dip/cast',cast,'info.dat'];
ctdinfile = [ctddir  'ctd_' cruise_ctd '_' ctdnum '_psal.nc'];

%addpath('/noc/users/pstar/di359/data/mexec_processing_scripts/',path);
jd0 = julian(year,1,0,0);
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

% --- initialize figures
figure(34);clf;hold on
col = 'brgkmcybrgkmcybrgkmc';
figure(35);clf;hold on
col = 'brgkmcybrgkmcybrgkmc';
figure(36);clf;hold on
col = 'brgkmcybrgkmcybrgkmc';
if find(id2==335)
    figure(37);clf;hold on
    col = 'brgkmcybrgkmcybrgkmc';
    figure(38);clf;hold on
    col = 'brgkmcybrgkmcybrgkmc';
end

% --- create log file ---
%fidlog = fopen([outpath,'microcat2rodb.log'],'a');
fidlog = fopen([outpath,'microcat2rodb.log'],'w');

% if id2(i)==335
    legend_handle = ['[';'[';'[';'[';'['];
% else
%     legend_handle = ['[';'[';'['];
% end
legend_string = [];
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
        infile = [inpath,sprintf('%4.4d',vec(i)),'_cal_dip_data.asc'];
    end
    if exist(infile) ~= 2
        infile = [inpath,sprintf('%4.4d',vec(i)),'_cal_dip_data.cnv'];
    end
    
    outfile = [outpath, 'cast', cast ,'_',sprintf('%4.4d',vec(i)),'.raw'];
    
    % --- convert from raw to rodb format ---
    
    microcat2rodb_3(infile,outfile,infofile,fidlog,'w',0)
    
    % --- load rodb data ---
    if ~id2(i)==335 % checks if not ODO microcat
        [yy,mm,dd,hh,c,t,p] = rodbload(outfile,'yy:mm:dd:hh:c:t:p');
    else % treats as ODO
        [yy,mm,dd,hh,c,t,p,ot,o2] = rodbload(outfile,'yy:mm:dd:hh:c:t:p:ot:o2');
    end
    %  if (i > 6 & i<=12)  lstr='--'; elseif i>12  lstr ='-.'; else lstr = '-'; end
    if (i > 7 & i<=14)  lstr='--'; elseif i>14  lstr ='-.'; else lstr = '-'; end
    
    disp(['plotting ',num2str(i),': s/n:',num2str(sn(i))])
    %    pause
    % correct timing errors for a number of microcats
    jd = julian(yy,mm,dd,hh)-jd0;
    
    ixjd = jd >jdstt & jd < jdend;
    c(~ixjd) = NaN;
    t(~ixjd) = NaN;
    p(~ixjd) = NaN;
    if id2(i)==335
        o2(~ixjd) = NaN;
        ot(~ixjd) = NaN;
    end
    
    macfigure(34)
    eval(['hc' num2str(i) '=plot(jd,c,[col(i),lstr])']); grid on
    macfigure(35)
    eval(['ht' num2str(i) '=plot(jd,t,[col(i),lstr])']); grid on
    macfigure(36)
    eval(['hp' num2str(i) '=plot(jd,p,[col(i),lstr])']); grid on
    if id2(i)==335
        macfigure(37)
        eval(['hw' num2str(i) '=plot(jd,ot,[col(i),lstr])']); grid on
        macfigure(38)
        eval(['ho' num2str(i) '=plot(jd,o2,[col(i),lstr])']); grid on
    end
    %    plot(jd,p,[col(i),'.']);
    if doctd % interpolate CTD onto microcat for a rough and ready mean diff
        pi = interp1(d.timeJ, d.press, jd);
        ti = interp1(d.timeJ, d.temp1, jd);
        % ci = interp1(d.timeJ, d.cond1*10., jd); % D382 had C in S/m
        ci = interp1(d.timeJ, d.cond1, jd); % JC103 has it in mS/cm
        if find(id2==335) % include oxygen
            oi = interp1(d.timeJ, d.oxygen_sbe2, jd); % secondary sensor is preferred sensor on JC103
        end
        dp = diff(pi) ;
        
        idp = find(dp(floor(length(dp)/2)-36:end-36) < 0.5);
        
        pdiff = ['mean p diff = ' num2str(nanmean(abs(p(idp) - pi(idp)))) ...
            ' --- ' num2str(vec(i))];
        cdiff = ['mean c diff = ' num2str(nanmean(abs(c(idp) - ci(idp)))) ...
            ' --- ' num2str(vec(i))];
        tdiff = ['mean t diff = ' num2str(nanmean(abs(t(idp) - ti(idp)))) ...
            ' --- ' num2str(vec(i))];
        if id2(i)==335
            odiff = ['mean o2 diff = ' num2str(nanmean(abs(o2(idp) - oi(idp)))) ...
                ' --- ' num2str(vec(i))];
        end
        
        disp(pdiff)
        disp(cdiff)
        disp(tdiff)
        if id2(i)==335
            disp(odiff)
        end
        
        %         fprintf(fidlog,'%s \n',pdiff);
        %         fprintf(fidlog,'%s \n',cdiff);
        %         fprintf(fidlog,'%s \n',tdiff);
        if id2(i)==335
            fprintf(fidlog,'%s \n',odiff);
        end
    end
    
    if i<length(vec)&~isempty(eval(['hc' num2str(i)]));
        %        legend_handle(1,:) = [legend_handle(1,:) ['hc' num2str(i) ',']];
        %        legend_handle(2,:) = [legend_handle(2,:) ['hp' num2str(i) ',']];
        %        legend_handle(3,:) = [legend_handle(3,:) ['ht' num2str(i) ',']];
%         if find(id2==335)
%             if ~exist(['ho' num2str(i)])
%                 legend_handle = [legend_handle ...
%                     ['hc' num2str(i) ',';'hp' num2str(i) ',';'ht' num2str(i) ',']];
%             else
                legend_handle = [legend_handle ...
                    ['hc' num2str(i) ',';'hp' num2str(i) ',';'ht' num2str(i) ','; 'hw' num2str(i) ',';'ho' num2str(i) ',']];
%             end
%         else
%             legend_handle = [legend_handle ...
%                 ['hc' num2str(i) ',';'hp' num2str(i) ',';'ht' num2str(i) ',']];
%         end
        legend_string = [legend_string '''' num2str(vec(i)) ''',' ];
    elseif i==length(vec)&~isempty(eval(['hc' num2str(i)]));
        %        legend_handle(1,:) = [legend_handle(1,:) ['hc' num2str(i) ']']];
        %        legend_handle(2,:) = [legend_handle(2,:) ['hp' num2str(i) ']']];
        %        legend_handle(3,:) = [legend_handle(3,:) ['ht' num2str(i) ']']];
%         if find(id2==335)
%             if ~exist(['ho' num2str(i-1)])
%                 legend_handle = [legend_handle ...
%                     ['hc' num2str(i) ']';'hp' num2str(i) ']';'ht' num2str(i) ']']];
%             else
                legend_handle = [legend_handle ...
                    ['hc' num2str(i) ']';'hp' num2str(i) ']';'ht' num2str(i) ']'; 'hw' num2str(i) ']';'ho' num2str(i) ']']];
%             end
%         else
%             legend_handle = [legend_handle ...
%                 ['hc' num2str(i) ']';'hp' num2str(i) ']';'ht' num2str(i) ']']];
%         end
        legend_string = [legend_string '''' num2str(vec(i)) '''' ];
        
    end;
    disp(['proceeding to next file '])
    
end % for i = 1:length(vec)

fclose(fidlog)


% --- tidy-up graphics, add CTD data ---

outfig = [outpath, 'cast', cast ,'_all'];
legend_string = [legend_string ', ''location'', ''EastOutside'''];
figure(34)
% num_legend(vec(:)',[],4)
eval(['legend(' legend_handle(1,:) ',' legend_string ')'])
ylabel('conductivity')
xlabel('yearday (relative to 2012/1/0 00:00)')
title(['CAST ' cast ' Calibration Dip'])
if doctd
    % plot(ctd_jd - jd0 ,ctd_cond*10,'k-')
    %  plot(d.timeJ,d.c0S_slash_m*10.,'k-')
    %     plot(d.timeJ+13,d.cond,'k-') % bim
    plot(d.timeJ,d.cond,'k-') % bim
    plot(d.timeJ,d.cond2,'color',[.3 .3 .3]); % testing for sensor drift
    title(['CAST ' cast  ' Calibration Dip (-k=CTD)'])
end
xlim([jdstt-0.01 jdend+0.01]);
orient tall
print(gcf,'-depsc',[outfig '_cond.ps'])
saveas(gcf,[outfig '_cond.fig'],'fig')
%print(['-f' num2str(gcf)],[outfig '_cond.fig'])

figure(35)
% num_legend(vec(:)',[],4)
eval(['legend(' legend_handle(2,:) ',' legend_string ')'])
ylabel('temperature')
xlabel('yearday (relative to 2012/1/0 00:00)')
title(['CAST ' cast ' Calibration Dip'])
if doctd
    % plot(ctd_jd - jd0,ctd_temp,'k-')
    %     plot(d.timeJ+13,d.temp1,'k-'); % bim
    plot(d.timeJ,d.temp2,'color',[.3 .3 .3]); % testing for sensor drift
    %  plot(d.timeJ,d.t090C,'k-');
    title(['CAST ' cast ' Calibration Dip (-k=CTD)'])
    title(['CAST ' cast ' Calibration Dip (-k=CTD)'])
end
xlim([jdstt-0.01 jdend+0.01]);
orient tall
print(gcf,'-depsc',[outfig '_temp.ps'])
saveas(gcf,[outfig '_temp.fig'],'fig')
%print(['-f' num2str(gcf)],[outfig '_temp.fig'])

figure(36)
% num_legend(vec(:)',[],4)
eval(['legend(' legend_handle(3,:) ',' legend_string ')'])
ylabel('pressure')
xlabel('yearday (relative to 2012/1/0 00:00)')
title(['CAST ' cast ' Calibration Dip'])
if doctd
    % plot(ctd_jd - jd0 ,ctd_press,'k-')
    %     plot(d.timeJ+13,d.press,'k-'); % bim
    plot(d.timeJ,d.press,'k-'); % bim
    title(['CAST ' cast ' Calibration Dip (-k=CTD)'])
end
xlim([jdstt-0.01 jdend+0.01]);
orient tall
print(gcf,'-depsc',[outfig '_pres.ps'])
saveas(gcf,[outfig '_pres.fig'],'fig')
%print(['-f' num2str(gcf)],[outfig '_pres.fig'])

if find(id2==335)
    figure(37)
    % num_legend(vec(:)',[],4)
eval(['legend(' legend_handle(4,:) ',' legend_string ')'])
    ylabel('oxygen temp.')
    xlabel('yearday (relative to 2012/1/0 00:00)')
    title(['CAST ' cast ' Calibration Dip'])
    if doctd
        % plot(ctd_jd - jd0 ,ctd_press,'k-')
        %         plot(d.timeJ+13,d.temp1,'k-'); % bim
        plot(d.timeJ,d.temp1,'k-'); % bim
        title(['CAST ' cast ' Calibration Dip (-k=CTD)'])
    end
    xlim([jdstt-0.01 jdend+0.01]);
    orient tall
    print(gcf,'-depsc',[outfig '_oxy_temp.ps'])
    saveas(gcf,[outfig '_oxy_temp.fig'],'fig')
    %print(['-f' num2str(gcf)],[outfig '_pres.fig'])
    
    figure(38)
    % num_legend(vec(:)',[],4)
eval(['legend(' legend_handle(5,:) ',' legend_string ')'])
    ylabel('oxygen')
    xlabel('yearday (relative to 2012/1/0 00:00)')
    title(['CAST ' cast ' Calibration Dip'])
    if doctd
        % plot(ctd_jd - jd0 ,ctd_press,'k-')
        plot(d.timeJ,d.oxygen_sbe2,'k-'); % bim
        %         plot(d.timeJ+13,d.oxygen_sbe1,'k-'); % bim
        title(['CAST ' cast ' Calibration Dip (-k=CTD)'])
    end
    xlim([jdstt-0.01 jdend+0.01]);
    orient tall
    print(gcf,'-depsc',[outfig '_oxy.ps'])
    saveas(gcf,[outfig '_oxy.fig'],'fig')
    %print(['-f' num2str(gcf)],[outfig '_pres.fig'])
end


disp(['number of MicroCATs processed = ' num2str(length(vec))])

figure(34)
%timeaxis([2008,1,0])
figure(35)
%timeaxis([2008,1,0])
figure(36)
%timeaxis([2008,1,0])
figure(37)
figure(38)

fclose all

