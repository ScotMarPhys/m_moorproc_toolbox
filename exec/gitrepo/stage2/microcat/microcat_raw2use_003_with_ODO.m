% MICROCAT_RAW2USE_003 is a script that performs stage2 processing
% on microcat data:
%      1. eliminates lauch and recovery periods
%      2. saves data to rodb file
%      3. creates data overview figures
% 
% uses timeaxis.m, auto_filt.m, julian.m 

% 11.01.01 Kanzow 
% 13.08.02 Kanzow : debugged
% Paths changed for D344 18 October 2009 P Wright
% 22.03.2010 ZB Szuts: modified for Oceanus 459
% 22/Nov/2015 : bim: DY039, updated to accept ODO oxygen temp and O2     

global MOORPROC_G
clearvars -except MOORPROC_G

% only mooring name and dates need to be modified, rest set in MOORPROC_G by
% startup{cruise}.m

moor = input('mooring deployment (e.g. ebh2_15_2022) to process:   ','s');
%moor = 'ebh2_15_2022';
% the start and end times of the time axis for plotting
%plot_interval = [2023 03 05; 2023 07 21];
%plot_interval = [2022 02 24; 2024 03 30];

cruise = MOORPROC_G.cruise;
operator = MOORPROC_G.operator;

% get paths to files
pd = moor_inoutpaths('microcat',moor);
basedir = MOORPROC_G.moordatadir;


% -------------------------------------------------------------------

mc_id    = [333 335 337] ;             % microcat id numbers


% --- get mooring information from infofile 

[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr] = ...
    rodbload(pd.infofile,['instrument:serialnumber:z:Start_Time:Start_Date:'...
                    'End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring']);

ii = find(id >= mc_id(1) & id<=mc_id(3));
sn = sn(ii);
z  = z(ii);
id = id(ii);

%sn = sn(1:2); z = z(1:2); id = id(1:2);

 
[z,zx] = sort(z);  % sort instruments by their depth
sn     = sn(zx);
id     = id(zx);



fid_stat = fopen(pd.stage2log,'w');
fprintf(fid_stat,'Processing steps taken by %s:\n',mfilename);
fprintf(fid_stat,'  1. eliminate lauching and recovery period\n');
fprintf(fid_stat,'  2. save data to rdb file\n');
fprintf(fid_stat,'\n Operated by:%s on %s\n',operator,datestr(clock));
fprintf(fid_stat,['        MicroCAT in Mooring ',moor,'\n\n\n']);
fprintf(fid_stat,'     ID    Depth   Start         End      Cycles  Spikes  Gaps   Mean     STD     Max     Min\n');

% ---- despike parameters
despike = 1; %***prompt
if despike %***track parameters (and whether used) for each cruise/record
    T_range = [-15 +15];
    C_range = [-30 +30];
    P_range = [-100 2000];
    O_range = [-50 550];
    % change from one point to the next can be no larger than 18*standard deviation of time series (after range editing)
    dT_range = 18; 
    dC_range = 18; 
    dP_range = 18; 
    dO_range = 18;
    nloop    = 3;
    %*** consider adding (even more conservative) parameters for oxygen?
end

dummy    = -9999;


%-----------------------------------------
% --- preprocessing loop -------------------
% ----------------------------------------

inst = 1;

jd_s  = julian(s_d(1),s_d(2),s_d(3),s_t(1)+s_t(2)/60);  % start time
jd_e  = julian(e_d(1),e_d(2),e_d(3),e_t(1)+e_t(2)/60);  % end time

for proc = 1 : length(sn)
    disp('plotting')

    infile  = fullfile(pd.stage1path,sprintf(pd.stage1form,sn(proc)));
    if exist(infile,'file')

        outfile = fullfile(pd.stage2path,sprintf(pd.stage2form,sn(proc)));

        inst = inst +1;

        if id(proc) == 335
            % load in ODO
            % OT - Oxygen temperature
            % O2 - Oxygen
            [YY,MM,DD,HH,C,T,P,OT,O2] = rodbload(infile,'YY:MM:DD:HH:C:T:P:OT:O2');
        else
            % standard SMP
            [YY,MM,DD,HH,C,T,P] = rodbload(infile,'YY:MM:DD:HH:C:T:P');
        end

        %------------------------------------------
        %----- cut off launching and recovery period
        %------------------------------------------
        disp('cut off launching and recovery period')

        jd = julian(YY,MM,DD,HH);
        ii = find(jd <= jd_e & jd >= jd_s );
        YY = YY(ii);   MM = MM(ii);   DD = DD(ii);
        HH = HH(ii);   c = C(ii);     t = T(ii);
        jd = jd(ii);
        if length(P) > 1,  p = P(ii); end
        if isempty(YY)
            warning('no data found in deployment period for file %d, %s; skipping',proc,infile)
            continue
        end
        % ODO - added bim DY039
        if id(proc) == 335, ot = OT(ii); o2= O2(ii); end

        cycles     = length(ii);
        Start_Date = [YY(1) MM(1) DD(1)];
        Start_Time = HH(1);
        End_Date = [YY(cycles) MM(cycles) DD(cycles)];
        End_Time = HH(cycles);

        if despike %***also plot with and without! 
            %------------------------------------------
            %--- despike ------------------------------
            %------------------------------------------
            disp('ddspike')

            [t,tdx,tndx] = ddspike(t,T_range,dT_range,nloop,'y',dummy);
            [c,cdx,cndx] = ddspike(c,C_range,dC_range,nloop,'y',dummy);
            if length(p) > 1
                [p,pdx,pndx] = ddspike(p,P_range,dP_range,nloop,'y',dummy);
            end
            if id(proc)==335
               [o2,o2dx,o2ndx] = ddspike(o2,O_range,dO_range,nloop,'y',dummy);
               [ot,otdx,otndx] = ddspike(ot,T_range,dT_range,nloop,'y',dummy);
            end

            % -----------------------------------------
            % ---  basic statistics -------------------
            % -----------------------------------------
            tstat = t ~= dummy;
            cstat = c ~= dummy;
            tstat = t(tstat);
            cstat = c(cstat);
            if id(proc)==335
                o2stat = o2 ~= dummy;
                o2stat = o2(o2stat);
            end
        end

        tm = meannan(t);
        cm = meannan(c);

        tsd= stdnan(t);
        csd= stdnan(c);

        tmx = max(t);
        cmx = max(c);
        tmn = min(t);
        cmn = min(c);

        if length(P) > 1
            if despike
                pstat = find(p ~= dummy);
            end
            pm  = meannan(p);
            psd = stdnan(p);
            pmx = max(p);
            pmn = min(p);
        end

        if id(proc) == 335 % ODO
            otm  = meannan(ot); otsd = stdnan(ot);
            otmx = max(ot); otmn = min(ot);
            o2m  = meannan(o2); o2sd = stdnan(o2);
            o2mx = max(o2); o2mn = min(o2);
        end

        %------------------------------------------
        %---- fill time gaps  with dummy
        %------------------------------------------

        disp(' fill time gaps  with dummy')

        djd = diff(jd);           % time step
        sr  = median(djd);        % sampling interval
        ii  = find(djd > 1.5*sr);  % find gaps
        gap = round(djd(ii)/sr)-1;
        addt= [];

        for i = 1 : length(gap)
            addt = [addt; [[1:gap(i)]*sr + jd(ii(i))]'];

        end

        [jd,xx] = sort([jd; addt]);   % add time
        ngap    = length(addt);       % number of time gaps
        gt      = gregorian(jd);
        YY = gt(:,1);   MM = gt(:,2);   DD = gt(:,3);
        if size(gt,2) == 6
            HH=hms2h(gt(:,4),gt(:,5),gt(:,6));
        else
            HH= gt(:,4);
        end


        t = [t;dummy*ones(ngap,1)]; t = t(xx);
        c = [c;dummy*ones(ngap,1)]; c = c(xx);
        if length(P) > 1
            p = [p;dummy*ones(ngap,1)]; p = p(xx);
        end
        %ODO
        if id(proc) == 335
            ot = [ot;dummy*ones(ngap,1)]; ot = ot(xx);
            o2 = [o2;dummy*ones(ngap,1)]; o2 = o2(xx);
        end

        %-----------------------------------------------------
        %  write output to logfile ---------------------------
        %-----------------------------------------------------

        disp(' write output to logfile')


        fprintf(fid_stat,'T   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d   %d         %d   %5.2f   %5.2f   %5.2f   %5.2f \n',...
            sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,tm,tsd,tmx,tmn');

        fprintf(fid_stat,'C   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d   %d         %d   %5.2f   %5.2f   %5.2f   %5.2f \n',...
            sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,cm,csd,cmx,cmn');

        if length(P) > 1
            fprintf(fid_stat,'P   %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d  %d        %d    %5.1f   %5.2f   %5.2f   %5.2f \n',...
                sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,pm,psd,pmx,pmn');
        end
        if id(proc) == 335 % ODO
            fprintf(fid_stat,'OT  %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d  %d        %d    %5.2f   %5.2f   %5.2f   %5.2f \n',...
                sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,otm,otsd,otmx,otmn');
            fprintf(fid_stat,'O2  %5.5d  %4.4d  %2.2d/%2.2d/%2.2d   %2.2d/%2.2d/%2.2d  %d        %d    %5.2f   %5.2f   %5.2f   %5.2f \n',...
                sn(proc),z(proc),Start_Date,End_Date,cycles,ngap,o2m,o2sd,o2mx,o2mn');
        end
        fprintf(fid_stat,'\n');

        %-----------------------------------
        %--- write data to rodb format -----
        %-----------------------------------

        disp(['writing data to ',outfile])

        rodboutvars = ['Latitude:Longitude:Columns:Start_Date:Start_Time:'...
            'SerialNumber:Mooring:WaterDepth:Instrdepth:'...
            'End_Date:End_Time'];
        if length(P) <= 1 % no pressure measurement
            sub =2;
            fort = '%4.4d   %2.2d   %2.2d   %8.5f   %6.4f   %6.4f';
            cols = 'YY:MM:DD:HH:T:C';
            rodbsave(outfile,rodboutvars,fort,...
                lat,lon,cols,Start_Date,Start_Time,sn(proc),mr,wd,...
                z(proc),End_Date,End_Time,[YY MM DD HH t c]);

        elseif id(proc) == 335 % ODO
            sub  = 5;
            fort = '%4.4d   %2.2d   %2.2d   %8.5f   %6.4f   %6.4f  %5.1f %8.5f %6.3f';
            cols = 'YY:MM:DD:HH:T:C:P:OT:O2';
            rodbsave(outfile,rodboutvars,fort,...
                lat,lon,cols,Start_Date,Start_Time,sn(proc),mr,wd,...
                z(proc),End_Date,End_Time,[YY MM DD HH t c p ot o2]);
        else
            sub  = 3;
            fort = '%4.4d   %2.2d   %2.2d   %8.5f   %6.4f   %6.4f  %5.1f';
            cols = 'YY:MM:DD:HH:T:C:P';
            rodbsave(outfile,rodboutvars,fort,...
                lat,lon,cols,Start_Date,Start_Time,sn(proc),mr,wd,...
                z(proc),End_Date,End_Time,[YY MM DD HH t c p]);
        end

        %%%%%%%%%% Graphics %%%%%%%%%%%%%%%%

        if exist('plot_interval','var') && ~isempty(plot_interval)
            jd1 = julian(plot_interval(1,:));
            jd2 = julian(plot_interval(2,:));
        else
            jd1 = jd_s-7;
            jd2 = jd_e+7;
        end

        figure(1);clf
        subplot(sub,1,1); m =~isnan(t)&t>dummy;

        plot(jd(m)-jd1,t(m))
        title(['MicroCAT s/n: ',num2str(sn(proc)),'; Target Depth: ',num2str(z(proc))])
        ylabel('Temperature [deg C]')
        grid on
        xlim([0 jd2-jd1])

        subplot(sub,1,2); m = ~isnan(c)&c>dummy;

        plot(jd(m)-jd1,c(m))
        ylabel('Conductivity [mS/cm]')
        grid on
        xlim([0 jd2-jd1])


        if sub == 3

            subplot(sub,1,3); m = ~isnan(p)&p>dummy;
            plot(jd(m)-jd1,p(m))

            ylabel('Pressure [dbar]')
            grid on
            xlim([0 jd2-jd1])

        end

        if sub == 5

            subplot(sub,1,3); ii = find(~isnan(p)&p>dummy);
            plot(jd(ii)-jd1,p(ii))

            ylabel('Pressure [dbar]')
            grid on
            xlim([0 jd2-jd1])

            subplot(sub,1,4); ii = find(~isnan(ot)&ot>dummy);
            plot(jd(ii)-jd1,ot(ii))

            ylabel('Oxygen temperature [degC]')
            grid on
            xlim([0 jd2-jd1])

            subplot(sub,1,5); ii = find(~isnan(o2)&o2>dummy);
            plot(jd(ii)-jd1,o2(ii))

            ylabel('Oxygen [umol/kg]')
            grid on
            xlim([0 jd2-jd1])

        end
        orient tall
        print(gcf,'-dpng','-r300',fullfile(pd.stage2figpath,sprintf(pd.stage2form,sn(proc))));

        sampling_rate = 1/median(diff(jd));
        
        % ==== TSD modified here (11/Jul/2024) to exclude dummy for filtering 
        % Caveat: this simple solution assumes that the gaps are small in
        %         time, which is not real. If there is a long gap, the
        %         filtering will propagate to next real value. The best
        %         solution would be to break the time series in chunks
        %         everytime we have a gap, and run the filter for each
        %         chunk.

        af  = find(t>-500);                                 % Find Dummy (-9999)
        tf  = auto_filt(t(af), sampling_rate, 1/2,'low',4); % Filter without Dummy
        tf2 = nan(length(t),1);                             % Create a new vector with NaNs
        tf2(af) = tf; tf = tf2; clear tf2 af                % Substitute the filtered values into the new vector
            
        af  = find(c>-500);                                 % Find Dummy (-9999)
        cf  = auto_filt(c(af), sampling_rate, 1/2,'low',4); % Filter without Dummy
        cf2 = nan(length(c),1);                             % Create a new vector with NaNs
        cf2(af) = cf; cf = cf2; clear cf2 af                % Substitute the filtered values into the new vector
        
        af  = find(p>-500);                                 % Find Dummy (-9999)
        pf  = auto_filt(p(af), sampling_rate, 1/2,'low',4); % Filter without Dummy
        pf2 = nan(length(p),1);                             % Create a new vector with NaNs
        pf2(af) = pf; pf = pf2; clear pf2 af                % Substitute the filtered values into the new vector

    	if id(proc) == 335
            % ==== TSD modified here to exclude dummy from filtering 
            af  = find(ot>-500);                                  % Find Dummy (-9999)
            otf  = auto_filt(ot(af), sampling_rate, 1/2,'low',4); % Filter without Dummy
            otf2 = nan(length(ot),1);                             % Create a new vector with NaNs
            otf2(af) = otf; otf = otf2; clear otf2 af             % Substitute the filtered values into the new vector
            
            af   = find(o2>-500);                                 % Find Dummy (-9999)
            o2f  = auto_filt(o2(af), sampling_rate, 1/2,'low',4); % Filter without Dummy
            o2f2 = nan(length(o2),1);                             % Create a new vector with NaNs
            o2f2(af) = o2f; o2f = o2f2; clear o2f2 af             % Substitute the filtered values into the new vector
    	end

        figure(2);clf
        ax(1)=subplot(sub,1,1);
        ii = find(~isnan(t)&t>dummy);

        plot(jd-jd1,tf)
        title(['2-day low-pass; MicroCAT s/n: ',num2str(sn(proc)),'; Target Depth: ',num2str(z(proc))])
        ylabel('Temperature [deg C]')
        grid on
        xlim([0 jd2-jd1])
        xlabel('Days since deployment');

        ax(2)=subplot(sub,1,2); ii = find(~isnan(c)&c>dummy);

        plot(jd-jd1,cf)
        ylabel('Conductivity [mS/cm]')
        grid on
        xlim([0 jd2-jd1])
        xlabel('Days since deployment');

        if sub == 3

            ax(3)=subplot(sub,1,3);

            plot(jd-jd1,pf)
            ylabel('Pressure [dbar]')
            grid on
            xlim([0 jd2-jd1])
            xlabel('Days since deployment');

        end
        if sub == 5

            ax(3)=subplot(sub,1,3);
            plot(jd(ii)-jd1,pf(ii))

            ylabel('Pressure [dbar]')
            grid on
            xlim([0 jd2-jd1])
            xlabel('Days since deployment');

            ax(4)=subplot(sub,1,4);
            plot(jd(ii)-jd1,otf(ii))

            ylabel('Oxygen temperature [degC]')
            grid on
            xlim([0 jd2-jd1])
            xlabel('Days since deployment');


            ax(5)=subplot(sub,1,5);
            plot(jd(ii)-jd1,o2f(ii))

            ylabel('Oxygen [umol/kg]')
            grid on
            xlim([0 jd2-jd1])
            xlabel('Days since deployment');

        end

        orient tall
        print(gcf,'-dpng','-r300',fullfile(pd.stage2figpath,[sprintf(pd.stage2form,sn(proc)) '_lowpass']));
        disp('pause (press any key to continue)'); pause

    end % if exist(infile)

end % for proc = 1 : length(sn),



