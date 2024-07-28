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
fprintf(fid_stat,'  2. Despike data and fill time gaps with dummy\n');
fprintf(fid_stat,'  3. save data to rdb file\n');
fprintf(fid_stat,'\n Operated by:%s on %s\n',operator,datestr(clock));
fprintf(fid_stat,['        MicroCAT in Mooring ',moor,'\n\n\n']);

manual_despiking='No spikes were removed manually \n';

% ---- ODO temperature sensor check parameters
dT_odo_tol = 0.01; % Tolerance for diff between sbe63 and sbe37 temp sensors
dT_odo_pct_bad_tol = 1; % How many percent bad points (diff above tolerance) will trigger warning
% ---- despike parameters
despike = 1; %***prompt
if despike %***track parameters (and whether used) for each cruise/record
    T_range = [-15 +15];
    C_range = [-30 +30];
    P_range = [-100 2000];
    O_range = [-50 550];
    % change from one point to the next can be no larger than 18*standard deviation of time series (after range editing)
    dT_range = 10; 
    dC_range = 10; 
    dP_range = 10; 
    dO_range = 10;
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

        %-----------------------------------------------
        %----- ODOs only: check offset between sbe37 
        %----- and sbe63 temp sensors, and ask if values
        %----- should be recalculated using sbe37 temp
        %-----------------------------------------------
        if id(proc) == 335
            
            dT_odo = OT - T;

            figure;
            plot(jd,dT_odo)
            hold on; grid on;
            title(['Temperature difference sbe63-sbe37; MicroCAT-ODO s/n: ',num2str(sn(proc)),'; Target Depth: ',num2str(z(proc))])
            xlabel('Days since deployment'); ylabel('Temperature difference [deg C]');

            ind_dT_odo_bad = find(abs(dT_odo)>=dT_odo_tol);
            pct_bad_dT_odo = length(ind_dT_odo_bad)/length(dT_odo)*100;

            if pct_bad_dT_odo>dT_odo_pct_bad_tol/100
                disp(['WARNING : large difference between sbe37 and sbe63 temperatures detected for ODO s/n ' num2str(sn(proc)) ' (' num2str(pct_bad_dT_odo,'%2.0f') '% of points over ' num2str(dT_odo_tol) ' degC difference theshold)'])
                redo_odo = input('Do you want to recaculate oxygen concentrations using sbe37 temperature (y/n)? [note: instrument raw data file will need to contain SBE63 phase delay voltage]    ','s');
                if redo_odo=='y'
                    % Backup original oxygen values for comparison
                    O2_original = O2;
                    clear O2;
                    % Recalculate oxygen
                    infile_raw = fullfile(pd.rawpath,[num2str(sn(proc)) '_data.cnv']);
                    O2 = microcat_odo_recalculate_oxy(infile_raw);
                    % Check plot
                    figure;
                    plot(jd,O2_original,'r')
                    hold on; grid on;
                    plot(jd,O2,'g')
                    title(['Oxygen concentration, original (red) vs recalculated (green); MicroCAT-ODO s/n: ',num2str(sn(proc)),'; Target Depth: ',num2str(z(proc))])
                    xlabel('Days since deployment'); ylabel('O_2 [/mumol/kg]');
                end
            end

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

        end

        close(989)



        % ----------------------------------------
        % ---  Manual De-spike -------------------
        % ----------------------------------------
        
        if exist('plot_interval','var') && ~isempty(plot_interval)
            jd1 = julian(plot_interval(1,:));
            jd2 = julian(plot_interval(2,:));
        else
            jd1 = jd_s-7;
            jd2 = jd_e+7;
        end

        % manual editing with brush
        n = 3; varnames = {'p','t','c'};
        x = jd-jd1; xl = x([1 end]);
        if id(proc)==335
            n = 5; varnames = [varnames,'ot','o2'];
        end
        %set up variables and initial plot
        editfig = figure;
        estr = '';
        for no = 1:n
            dsplotdata.(varnames{no}) = eval(varnames{no});
            dsplotdata.(varnames{no})(dsplotdata.(varnames{no})<-999) = nan;
            numspikes.(varnames{no}) = 0;
            ha(no) = subplot(n,1,no);
            hl(no) = plot(x,dsplotdata.(varnames{no})); xlim(xl); grid on
            ylabel(varnames{no})
            estr = [estr sprintf('%d = %s,',no,varnames{no})];
        end
        estr = estr(1:end-1);
        editing=input('Manual editing for plotting required? 0=no, 1=yes: ');
        while editing==1
            delete(hl); try delete(hle); catch; end
            %replot -- this is necessary if looping more than once
            for no = 1:n
                axes(ha(no))
                hl(no) = plot(x,dsplotdata.(varnames{no})); 
                xlim(xl); grid on; ylabel(varnames{no})
            end
            zoomx = input('use figure buttons to zoom then type z to zoom all to current xlims, or enter to skip  ','s');
            if ~isempty(zoomx) && strcmp(zoomx,'z')
                xl = get(gca,'xlim');
                set(ha(:),'xlim',xl)
            end
            varnum = input(sprintf('which variable do you want to edit (%s or enter to stop editing)? ',estr));
            if isempty(varnum)
                editing = 0; continue
            end
            fprintf([' Manually select data on subplot %d (from top) using data brushing tool \n...' ...
                'and create variable "brushedData" from bad data. \n ...' ...
                'Need to right-click after data selection and *export brushed* as "brushedData".\n ...' ...
                'To select more than one period, hold shift and select the other periods.\n'],varnum);
            clear brushedData
            brush(editfig); pause
            if exist('brushedData','var')
                m = ismember(x,brushedData(:,1));
                dat = dsplotdata.(varnames{varnum});
                dat(ismember(x,brushedData(:,1))) = nan;
                dsplotdata.(varnames{varnum}) = dat;
                numspikes.(varnames{varnum}) = size(brushedData,1);
                axes(ha(varnum)); hold on
                hle = plot(x,dsplotdata.(varnames{varnum}),'k');
                legend('before','edited')
            end
            editing = input('Make more edits?  0=no, 1=yes:  ');
            manual_despiking = 'Manual despiking applied \n';
        end
        for no = 1:n
            dsplotdata.(varnames{no})(isnan(dsplotdata.(varnames{no}))) = dummy;
            eval([varnames{no} ' = dsplotdata.(varnames{no});']);
        end
        close(editfig); clear dsplotdata


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
	    fprintf(fid_stat,'Serialnumber %d \n',sn(proc));
        fprintf(fid_stat,'Infile %s \n',infile);
        fprintf(fid_stat,'Outfile %s \n',outfile);
	    fprintf(fid_stat,'Operation interval: %s  to  %s\n', ...
            datestr(gregorian(jd(1))),datestr(gregorian(jd(end)) ));
        fprintf(fid_stat,'%s',manual_despiking);
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
        print(gcf,'-dpng','-r300',fullfile(pd.stage2figpath,[sprintf(pd.stage2form,sn(proc)),'.png']));

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
        print(gcf,'-dpng','-r300',fullfile(pd.stage2figpath,[sprintf(pd.stage2form,sn(proc)) '_lowpass.png']));
        disp('pause (press any key to continue)'); pause
        close all
    end % if exist(infile)

end % for proc = 1 : length(sn),

