% function seaphox_raw2use_01(moor,'procpath',procpath,'outpath',outpath,'plot_interval',plot_interval)
%
% required inputs: moor - mooring name e.g 'eb3_1_200406'
%
% optional inputs: procpath - if not using standard paths for info.dat file
%                  outpath - if not using standard paths 
%                  cruise - if running on historical cruise, otherwise will
%                        take cruise name from mfilename path. Also used if
%                        running for CTD cal dip so knows where to loog for
%                        data files
%                  plot_interval - matrix of start and end dates for plot
%                                  e.g. [2004 02 01 00; 2005 06 01 00]
%                                  dates are:- yyyy mm dd hh
%
% features
%      1. eliminate launching and recovery period
%      2. save data to rodb file
%      3. create data overview sheet
%
% uses timeaxis.m, auto_filt.m, julian.m, rodbload.m, rodbsave.m

% 18/11/2015 - DR modified from nortek_raw2use_02

function seaphox_raw2use_01(moor,varargin)
  
% if nargin==0
%     help seaphox_raw2use_01;
%     return
% end
% 
% if strfind(mfilename('fullpath'),'/Volumes/rpdmoc/rapid/data/exec/')
%     % using Mac with mount to rpdmoc either on a cruise or at NOC
%     basedir = '/Volumes/rpdmoc/rapid/data/';
% elseif strfind(mfilename('fullpath'),'/Volumes/pstar/rpdmoc/rapid/data/exec/')
%     % using Mac with mount to rpdmoc either on a cruise or at NOC
%     basedir = '/Volumes/pstar/rpdmoc/rapid/data/';
% elseif strfind(mfilename('fullpath'),'/local/users/pstar/rpdmoc/rapid/data/exec/')
%     % using seagoing workstation
%     basedir = '/local/users/pstar/rpdmoc/rapid/data/';
% else % using NOC network
%     basedir = '/noc/mpoc/rpdmoc/rapid/data/';
% end


global MOORPROC_G

basedir = MOORPROC_G.moordatadir;
cruise= MOORPROC_G.cruise;

if nargin==0
    help seaphox_raw2use_01
    return
end

%if nargin==0
%    help(mfilename)
%    return
%end

if strcmpi(moor(1:4),'cast')
    % indicates is a CTD cal dip
    cal_dip=1;
else
    cal_dip=0;
end

% check for optional arguments
a=find(strcmp('procpath',varargin));
if a>0
    procpath=char(varargin(a+1));
else
    procpath=[basedir '/proc'];
end

a=find(strcmp('cruise',varargin));
if a>0
    cruise=char(varargin(a+1));
else
    script=mfilename('fullpath');
    cruise=script(strfind(script,'/data/exec/')+11:strfind(script,'/stage2/seaphox/')-1);
end

if cal_dip==1
        inpath=[basedir '/proc_calib/' cruise '/cal_dip/seaphox/' moor '/'];
        infofile =[basedir '/proc_calib/' cruise '/cal_dip/' moor 'info.dat'];
    else
        inpath=[basedir '/proc/' moor '/seaphox/'];
        infofile =[basedir '/proc/' moor '/' moor 'info.dat'];
end

a = find(strcmp('outpath', varargin));
if a>0
    outpath=char(varargin(a+1));
else
    outpath = inpath;
end

a = find(strcmp('plot_interval', varargin));
if a>0
    plot_interval=varargin(a+1);
    plot_interval2=zeros(2,4); 
    plot_interval2(1,1:4)=plot_interval{1}(1,1:4); 
    plot_interval2(2,1:4)=plot_interval{1}(2,1:4);
    plot_interval=plot_interval2;
else
    plot_interval=0;
end

[gash, operator]=system('whoami');  % This line will not work if run from a PC. May need to edit it out.

[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

vec=find(id==375); % Seaphox id number

sn=sn(vec);
z = z(vec);

fid_stat= fopen([outpath moor '_seaphox_stage2.log'],'a');
fprintf(fid_stat,['Processing steps taken by ' mfilename ':\n']);
fprintf(fid_stat,'  1. eliminate launch and recovery period\n');
fprintf(fid_stat,'  2. resave data to rodb file\n');
fprintf(fid_stat,'\n Operated by:%s on %s\n',operator,datestr(clock)); 
fprintf(fid_stat,['        SeapHOx in Mooring ',moor,'\n\n\n']);

dummy    = -9999;

% Determine plot_interval if not input to function
if plot_interval==0
    plot_interval = zeros(2,4);
    plot_interval(1,1) = s_d(1); plot_interval(1,2) = s_d(2)-1; plot_interval(1,3) = 1; plot_interval(1,4) = 0;
    plot_interval(2,1) = e_d(1); plot_interval(2,2) = e_d(2)+1; plot_interval(2,3) = 1; plot_interval(2,4) = 0;
    if plot_interval(1,2)==0
        plot_interval(1,2)=12; plot_interval(1,1)=plot_interval(1,1)-1;
    end
    if plot_interval(2,2)==13
        plot_interval(2,2)=1; plot_interval(2,1)=plot_interval(2,1)+1;
    end
end


%-----------------------------------------
% --- preprocessing loop -------------------
% ----------------------------------------

jd_s  = datenum(s_d(1),s_d(2),s_d(3),s_t(1),s_t(2)/60,0);  % start time
jd_e  = datenum(e_d(1),e_d(2),e_d(3),e_t(1),e_t(2)/60,0);  % end time


for proc = 1 : length(vec)
    
    indep  = z(proc);
    infile  = [inpath, moor,'_',num2str(sn(proc)),'.raw'];
  
    if exist(infile)==0
        % check if sprintf formatted serial number to 3 digits instead of 1
        infile2 = [inpath,moor,'_',sprintf('%03d',sn(proc)),'.raw'];
        if exist(infile2)==0
            disp(['infile: ' infile ' does not exist.'])
            disp(['infile: ' infile2 ' does not exist.'])
        else
            infile=infile2;
        end
    end
    if exist(infile)==0 %do nothing as now shouldn't end up in this part of loop as covered above
    elseif exist(infile)   > 0
        rodbfile= [moor,'_',sprintf('%03d',sn(proc)),'.use']; 
        outfile = [outpath,rodbfile];
        fprintf(fid_stat,'Serialnumber %d \n',sn(proc));
        fprintf(fid_stat,'Infile %s \n',infile);
        fprintf(fid_stat,'Outfile %s \n',outfile);
        % seaphox data is slightly different format for V1 versus V2
        % instruments, so need to check  which variables present in file
        fid=fopen(infile);
        found_it=0;
        while found_it==0
            headerline=fgetl(fid);
            if contains(headerline,'Columns')
                found_it=1;
                if contains(headerline,'BATT1:BATT2')
                    % then is file produced from V1 SeapHOx (V2s don't log
                    % battery voltage)
                    version=1;
                else
                    version=2;
                end
                fclose(fid);
            end
        end
        
        if version==1
            columns = 'YY:MM:DD:HH:PH:PH_V:T:P:O:S:BATT1:BATT2';
            variables='[YY,MM,DD,HH,ph,ph_v,t,p,o,s,batt1,batt2]';
        elseif version==2
            columns = 'YY:MM:DD:HH:PH:PH_V:T:C:P:O:S';
            variables='[YY,MM,DD,HH,ph,ph_v,t,c,p,o,s]';
        end
        
        eval([variables '=rodbload(infile,columns);'])
        %[YY,MM,DD,HH,ph,ph_v,t,p,o,s,batt1,batt2] = ...
        %    rodbload(infile,columns);

        %------------------------------------------ 
        %----- cut off launching and recovery period
        %------------------------------------------
        disp('cut off launching and recovery period')

        jd  = datenum(YY,MM,DD,HH,0,0);
        %gregorian(jd);
%         jd_e
%         gregorian(jd_e)
%         jd_s
%         gregorian(jd_s)
        
        ii  = find(jd <= jd_e & jd >= jd_s );

        YY=YY(ii);MM=MM(ii);DD=DD(ii);HH=HH(ii);
        ph=ph(ii); ph_v=ph_v(ii);
        t=t(ii);p=p(ii);o=o(ii);s=s(ii);
        if version==1
            batt1=batt1(ii);batt2=batt2(ii);
        elseif version==2;
            c=c(ii);
        end
        jd  = jd(ii); 

        cycles     = length(ii);
        Y = [YY(1) MM(1) DD(1)];
        Start_Time = HH(1);
        End_Date = [YY(cycles) MM(cycles) DD(cycles)];
        End_Time = HH(cycles);     


        %------------------------------------------
        %---- fill time gaps  with dummy
        %------------------------------------------
        disp(' fill time gaps  with dummy')

        djd = diff(jd);           % time step  
        sr  = median(djd);        % sampling interval
        ii  = find(djd > 1.5*sr);  % find gaps
        gap = round(djd(ii)/sr)-1;
        addt= []; 

        

        %-----------------------------------------------------
        %  write output to logfile ---------------------------
        %-----------------------------------------------------

        fprintf(fid_stat,'Operation interval: %s  to  %s\n', ... 
             datestr(gregorian(jd(1))),datestr(gregorian(jd(end)) ));
        fprintf(fid_stat,'\n');

        %-----------------------------------  
        %--- write data to rodb format -----
        %-----------------------------------

        disp(['writing data to ',outfile]) 
        if version==1
            fort =['%4.4d %2.2d %2.2d  %6.4f  %7.5f %10.8f  %6.4f %7.3f %6.3f  %6.4f %5.3f %5.3f'];
            rodbsave(outfile,...
          'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
           fort,...
          lat,lon,columns,s_d,Start_Time,sn(proc),mr,wd,indep,End_Date,End_Time,...
          [YY,MM,DD,HH,ph,ph_v,t,p,o,s,batt1,batt2]);
        elseif version==2
            fort =['%4.4d %2.2d %2.2d  %6.4f  %6.4f %8.6f  %6.4f %6.4f %7.3f  %5.3f %6.4f'];
            rodbsave(outfile,...
          'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
           fort,...
          lat,lon,columns,s_d,Start_Time,sn(proc),mr,wd,indep,End_Date,End_Time,...
          [YY,MM,DD,HH,ph,ph_v,t,c,p,o,s]);
        end
        
        
      
      
        
        %%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%
        %jd0 = julian(-1,1,1,24);
        %jd1 = julian(plot_interval(1,:))-jd0;
        %jd2 = julian(plot_interval(2,:))-jd0; 
        jd1 = datenum(plot_interval(1,1),plot_interval(1,2),plot_interval(1,3),plot_interval(1,4),0,0);
        jd2 = datenum(plot_interval(2,1),plot_interval(2,2),plot_interval(2,3),plot_interval(2,4),0,0);
        sampling_rate = 1/median(diff(jd));

        STR = ['pH                        ';
               'MC Temperature [deg C]    ';
               'MC Pressure [dbar]        ';
               'MC Oxygen [umol/kg]       ';
               'MC Salinity               ';];
        VAR1= {'ph','t ','p ','o ','s '};
        panels=5;
        if version==1
            VAR2= {'ph_v','batt1','batt2'};
            STR2 =['pH voltage [V]            ';
                   'voltage [V]               ';];
               panels2=3;
        elseif version==2
            VAR2= {'ph_v'};
            STR2 =['pH voltage [V]            '];
            panels2=1;
        end
        
                
        figure(1);clf
        for sub = 1 : panels
          eval(['var1 = ',VAR1{sub},';'])
          var2=[];
          var3=[];
%           if strfind(STR(sub,:),'pH')
%               var1(var1>9)=NaN;
%               var1(var1<6)=NaN;
%           end
          ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR(sub,:),sub,[jd1 jd2],'n',panels);
          xlim([jd1 jd2])
          datetick('x','dd-mmm','keepticks','keeplimits')
        end

        subplot(5,1,1)
        title(['SeapHOx s/n: ',num2str(sn(proc)), ...
             '; Target Depth: ',num2str(indep)])

        orient tall
        eval(['print -depsc2 -tiff ',outfile,'.eps']) 


        figure(2);clf

        for sub = 1 : 5
          eval(['var1 = ',VAR1{sub},';'])
          var2=[];
          var3=[];
          % if pH variable then remove spikes first - filtering doesn't
          % work when have burst sampling with varying sampling rate
          % between individual samples and the whole burst
          if strfind(STR(sub,:),'pH')
              var1(var1>9)=NaN;
              var1(var1<6)=NaN;
          end
          ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR(sub,:),sub,[jd1 jd2],'y',panels);
          xlim([jd1 jd2])
          datetick('x','dd-mmm','keepticks','keeplimits')

        end

        subplot(5,1,1)
        title(['SeapHOx s/n: ',num2str(sn(proc)), ...
               '; Target Depth: ',num2str(indep)])
        orient tall
        eval(['print -depsc ',outfile,'.filtered.eps']) 
        
        
        % plot of diagnostics info
        figure(3);clf
        subplot(3,1,1)
        eval(['var1 = ',VAR2{1},';'])
        var2=[];
        var3=[];
        ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR2(1,:),1,[jd1 jd2],'n',panels2);
        xlim([jd1 jd2])
        datetick('x','dd-mmm','keepticks','keeplimits')
        
        legend({'ph voltage'})
        
        if version==1;
            subplot(3,1,2)
            eval(['var1 = ',VAR2{2},';'])
            eval(['var2 = ',VAR2{3},';'])
            var3=[];
            ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR2(2,:),2,[jd1 jd2],'n',panels2);
            xlim([jd1 jd2])
            datetick('x','dd-mmm','keepticks','keeplimits')

            legend({'main batt.','isolated batt.'})
        end
        
        suptitle(['SeapHOx s/n: ',num2str(sn(proc)), ...
                 '; Target Depth: ',num2str(indep)])

          orient tall
        eval(['print -depsc2 -tiff ',outfile,'_ancillary.eps'])
        end % if exist(infile)==0
  end  % for proc=1:length(combo_sn)+length(individual_sn) loop
end % function
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function ok = plot_timeseries(jd,var1,var2,var3,sr,str,sub,jdlim,filt,panels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot time series


 % jd0 = julian(-1,1,1,24);

  i1    = find(~isnan(var1) & var1~=0);
 
  if strcmp(filt,'y')
    var1  = auto_filt(var1(i1),sr,1/2,'low',4);
  else
    var1  = var1(i1);
  end 

  if ~isempty(var2)if jdlim(2)-jdlim(1)>365/12;
    datetick('x',12)
  elseif jdlim(2)-jdlim(1)<1;
      datetick('x','HH:MM','keeplimits')
  else
      datetick('x','dd-mmm')
  end
    i2    = find(~isnan(var2) & var2~=0);
    if strcmp(filt,'y')
      var2  = auto_filt(var2(i2),sr,1/2,'low',4);
    else
      var2  = var2(i2);
    end
  end

  if ~isempty(var3)
    i3    = find(~isnan(var3) & var3~=0);
    if strcmp(filt,'y')
      var3  = auto_filt(var3(i3),sr,1/2,'low',4);
    else
      var3  = var3(i3);
    end
  end
  
  subplot(panels,1,sub);
         
  %plot(jd(i1)-jd0,var1)
  plot(jd(i1),var1)
  hold on

  if ~isempty(var2)
    %plot(jd(i2)-jd0,var2,'r')
    plot(jd(i2),var2,'r')
  end
 
  if ~isempty(var3)
    %plot(jd(i3)-jd0,var3,'g')
    plot(jd(i3),var3,'g')
  end
  
  ylabel(str)
  grid on
  xlim([jdlim])
%   if jdlim(2)-jdlim(1)>365/12;
%     datetick('x',12)
%   elseif jdlim(2)-jdlim(1)<1;
%       datetick('x','HH:MM','keeplimits')
%   else
%       datetick('x','dd-mmm','keepticks','keeplimits')
%   end
% 
%   xlim([jdlim(1) jdlim(2)])
      
  ok=1;

  end