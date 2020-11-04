% stage1 routine to load Seaguard data to rodb format.
% written by Darren Rayner - cruise D344 - 28/10/09

% call the routine with seaguard_call.m so don't need to adjust paths in
% this routine.
function seaguard2rodb(infofile,cruise,cal_dip);
    %inputs:-
    % infofile = full path and filename of info.dat file
    % cruise = cruise name e.g. 'd344'
    % cal_dip = logical expression. 1=cal_dip data, 0=mooring data

    Instrument='Seaguard'; % for rodb header
    cols       = 'YY:MM:DD:HH:U:V:CS:CD:CSSD:MSS:HDG:PIT:ROL:T:C:TC:P:TP:IPOW';    % column info for rodb header
    fort        = '%4.4d %2.2d %2.2d %7.5f %6.3f %6.3f %6.3f %5.1f %6.2f %6.3f %5.1f %4.1f %4.1f %6.3f %6.3f %6.3f %5.1f %6.3f %5.3f'; %data output format

    % check if infofile exists
    if exist(infofile) ~= 2
       disp(['infofile:  ',infofile,' does not exist'])
       return
    end
    %---------------------------------------------------------------
    % load info.dat file
    %---------------------------------------------------------------
        infovar ='Mooring:Latitude:Longitude:Waterdepth:id:serialnumber:z:StartDate:StartTime:EndDate:EndTime'; 

        [mo,la,lo,wd,id,sn,z,sdate,stime,edate,etime]=rodbload(infofile,infovar); 
        if isempty(id) | isnan(id)
          infovar ='Mooring:Latitude:Longitude:Waterdepth:instrument:serialnumber:z:StartDate:StartTime:EndDate:EndTime'; 
          [mo,la,lo,wd,id,sn,z,sdate,stime,edate,etime]=rodbload(infofile,infovar); 
        end
        if iscell(mo)
          mo = deal(mo{:}); % convert cell array
        end
        ii=find(id==301);
        z=z(ii); sn=sn(ii); 
    % end of reading info.dat file
    %----------------------------------------
    a=strfind(infofile, 'proc');
    if cal_dip==1
        outpath=[infofile(1:a-1) 'proc_calib/' cruise '/cal_dip/seaguard/cast' mo(end) '/'];
    elseif cal_dip==0
        outpath=[infofile(1:a-1) 'proc/',mo,'/seaguard/'];
    end
    log=[outpath 'seaguard_stage1.log'];
    fid_log=fopen(log,'a');  % opens log file for appending to
    if cal_dip==1
        in_dir=[infofile(1:a-1) 'raw/' cruise '/seaguard_cal_dip/cast' mo(end) '/'];
    elseif cal_dip==0
        in_dir=[infofile(1:a-1) 'raw/' cruise '/seaguard/'];
    end

    for i=1:length(sn)
        outfile=[outpath mo '_' num2str(sn(i)) '.raw'];
        %------------------------------------------------------------------------------------------------
        % determine filenames of different sensors from the folder containing the
        % unique instrument data eg. folder 114 corresponds to data from instrument
        % 114 - with 5 data files per instrument exported from the Seaguard Studio software
        in_dir2=[in_dir num2str(sn(i)) '/'];
        c_file=ls([in_dir2 'Conductivity*']);
        c_file=['''' c_file(1:end-1) '''']; % cut off CR character at end of ls response
        t_file=ls([in_dir2 'Temperature*']);
        t_file=['''' t_file(1:end-1) '''']; % cut off CR character at end of ls response
        p_file=ls([in_dir2 'Pressure*']);
        p_file=['''' p_file(1:end-1) '''']; % cut off CR character at end of ls response
        cm_file=ls([in_dir2 'DCS*']);
        cm_file=['''' cm_file(1:end-1) '''']; % cut off CR character at end of ls response
        sys_file=ls([in_dir2 'System*']);
        sys_file=['''' sys_file(1:end-1) '''']; % cut off CR character at end of ls response

        % load data from files determined above
        c_data=eval(['csvread(' c_file ',1,0)']);
        t_data=eval(['csvread(' t_file ',1,0)']);
        cm_data=eval(['csvread(' cm_file ',1,0)']);
        p_data=eval(['csvread(' p_file ',1,0)']);
        sys_data=eval(['csvread(' sys_file ',1,0)']);

        % allocate correct columns to correct variable names
        C=c_data(:,3);
        TC=c_data(:,4);
        T=t_data(:,3);
        P=p_data(:,3)/10; % pressure is in kPa so need to divide by 10
        TP=p_data(:,4);
        time=cm_data(:,2);
        CS=cm_data(:,3);
        CD=cm_data(:,4);
        U=cm_data(:,5);
        V=cm_data(:,6);
        HDG=cm_data(:,7);
        PIT=cm_data(:,8);
        ROL=cm_data(:,9);
        CSSD=cm_data(:,10);
        MSS=cm_data(:,11);
        IPOW=sys_data(:,3);

        % calculate time in rodb format from seaguard format
        jd=time+1721059; % 1721059 corresponds to 1/1/0000
                            % Seaguard date is in days AD.
        gtime=gregorian(jd);

    if isunix
        [gash, operator]=system('whoami');  % This line will not work if run from a PC. May need to edit it out.
    else
        operator = 'unknown';
    end
    fprintf(fid_log,'Transformation of Seaguard csv data files to rodb format \n');
    fprintf(fid_log,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));

    fprintf(fid_log,'Mooring   %s \n',mo);
    fprintf(fid_log,'Latitude  %6.3f \n',la);
    fprintf(fid_log,'Longitude %6.3f \n',lo);
    fprintf(fid_log,'Serial Number: %d\n\n',sn(i));
    fprintf(fid_log,'System File %s\n',sys_file);
    fprintf(fid_log,'Curren sensor File %s\n',cm_file);
    fprintf(fid_log,'Conductivity File %s\n',c_file);
    fprintf(fid_log,'Temperature File %s\n',t_file);
    fprintf(fid_log,'Pressure File %s\n',p_file);
    fprintf(fid_log,'Outfile: %s\n\n',outfile);
    fprintf(fid_log,'Start Date: %4.0d/%02.0f/%02.0f\n',gtime(1,1), gtime(1,2), gtime(1,3));
    fprintf(fid_log,'Start Time: %02.0f:%02.0f:%02.0f\n',gtime(1,4), gtime(1,5), gtime(1,6));
    fprintf(fid_log,'End Date: %4.0f/%02.0f/%02.0f\n',gtime(end,1), gtime(end,2), gtime(end,3));
    fprintf(fid_log,'End Time: %02.0f:%02.0f:%02.0f\n',gtime(end,4), gtime(end,5), gtime(end,6));
    fprintf(fid_log,'Number of Samples: %d\n',length(jd));

    % ==============================================
    % START OF SAVE TO RODB FORMAT
    disp(['writing data to ',outfile]) 
    'YY:MM:DD:HH:U:V:CS,CD,CSSD,MSS,HDG,PIT,ROL,T:C:TC:P:TP:IPOW';    % column info for rodb header
    gdate=gregorian(jd);
    gdate=[gdate(:,1:3) hms2h(gdate(:,4:6))];
    data = [gdate U V CS CD CSSD MSS HDG PIT ROL T C TC P TP IPOW]; 
    
           
    rodbsave(outfile,...
           'Latitude:Longitude:Columns:SerialNumber:Mooring:WaterDepth:Instrdepth:StartDate:StartTime:EndDate:EndTime',...
             fort,...
             la,lo,cols,sn(i),mo,wd,z,sdate,stime,edate,etime,...
             data);

    end
end