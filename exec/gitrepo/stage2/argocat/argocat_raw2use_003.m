% function argocat_raw2use_003(moor,'procpath',procpath,'outpath',outpath)
%
% basic preprocessing for combined Sontek Aargonaut MD data
% either combined with Microcat or stand alone
% 
% required inputs: moor - mooring name e.g 'eb3_1_200406'
%
% optional inputs: procpath - path to proc directory if not using standard
%                             paths
%                  outpath - path to output processed files to - otherwise
%                            outputs to directory function run from
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

% 14.04.05 Kanzow 
% 25.04.05 Kanzow; debugged    
% 13/12/06 DR - converted to function with path inputs
% 13/12/06 DR - modified to read files with noise, snr and pgp columns
% 13/12/06 DR - modified to work with unpaired argonauts as well
% 12/11/12 DR - updated to produce plots of diagnostic data too (velocity
%               standard deviation, beam signal strength, noise and signal
%               to noise ratio.) This highlights beam failures better.

function argocat_raw2use_003(moor,varargin)
  
% --- get mooring information from infofile 

if nargin==0
    help argocat_raw2use_003;
    return
end

% check for optional arguments
a=strmatch('procpath',varargin,'exact');
if a>0
    procpath=char(varargin(a+1));
else
    procpath='/noc/users/pstar/rpdmoc/rapid/data/moor/proc/';
end
a=strmatch('outpath',varargin,'exact');
if a>0
    outpath=char(varargin(a+1));
else
    outpath = ['/noc/users/pstar/rpdmoc/rapid/data/moor/proc/' moor '/arg/'];
end
a=strmatch('plot_interval',varargin,'exact');
if a>0
    plot_interval=varargin(a+1);
else
    plot_interval=0;
end

[gash, operator]=system('whoami');  % This line will not work if run from a PC. May need to edit it out.

% -- set path for data input
inpath  = [procpath '/' moor '/'];

% --- get moring information from infofile 
infofile =[procpath '/' moor '/' moor 'info.dat'];

%plot_interval = [2004 02 01 00;   % start time of time axis on plot
%		 2005 06 01 00];  % end time of time axis on plot

[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

combo=find(id==366337); % Argonauts coupled with Microcats 
individual=find(id==366); % uncoupled Argonauts

individual_sn=sn(individual);
individual_z = z(individual);
combo_sn=sn(combo);
combo_z  =  z(combo);

for i=1:length(individual)
    for j=1:length(combo)
        if individual_sn(i)==sn(combo(j))
            individual_sn(i)=NaN;
        end
    end
end

a=find(~isnan(individual_sn));
individual_sn=individual_sn(a);
individual_z=individual_z(a);

fid_stat= fopen([outpath moor '_Argonaut_stage2_log'],'w');
fprintf(fid_stat,['Processing steps taken by ' mfilename ':\n']);
fprintf(fid_stat,'  1. eliminate lauching and recovery period\n');
fprintf(fid_stat,'  2. save data to rdb file\n');
fprintf(fid_stat,'\n Operated by:%s on %s\n',operator,datestr(clock)); 
fprintf(fid_stat,['        Argonaut/MicroCATs in Mooring ',moor,'\n\n\n']);

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

jd_s  = julian(s_d(1),s_d(2),s_d(3),s_t(1)+s_t(2)/60);  % start time
jd_e  = julian(e_d(1),e_d(2),e_d(3),e_t(1)+e_t(2)/60);  % end time


all_z=[combo_z individual_z];
all_sn=[combo_sn individual_sn];

for proc = 1 : length(combo_sn)+length(individual_sn)
    if find(intersect(all_sn(proc),combo_sn))
        columns = 'YY:MM:DD:HH:T:TCAT:P:PCAT:C:U:V:W:HDG:PIT:ROL:USD:VSD:WSD:USS:VSS:WSS:HDGSD:PITSD:ROLSD:IPOW:UNOISE:VNOISE:WNOISE:USNR:VSNR:WSNR:PGP';
    elseif find(intersect(all_sn(proc),individual_sn))
        columns = 'YY:MM:DD:HH:T:P:U:V:W:HDG:PIT:ROL:USD:VSD:WSD:USS:VSS:WSS:HDGSD:PITSD:ROLSD:IPOW:UNOISE:VNOISE:WNOISE:USNR:VSNR:WSNR:PGP';
    end
    
  indep  = all_z(proc);

  infile  = [inpath,'arg/',moor,'_',num2str(all_sn(proc)),'.raw'];
  
  if exist(infile)==0
      disp(['infile: ' infile ' does not exist.'])
      
  elseif exist(infile)   > 0 
 
    rodbfile= [moor,'_',num2str(all_sn(proc)),'.use']; 
    outfile = [outpath,rodbfile];

    fprintf(fid_stat,'Serialnumber %d \n',all_sn(proc));
    fprintf(fid_stat,'Infile %s \n',infile);
    fprintf(fid_stat,'Outfile %s \n',outfile);
    
    if find(intersect(all_sn(proc),combo_sn))
        [YY,MM,DD,HH,t,tc,p,pc,c,u,v,w,...
         hdg,pit,rol,usd,vsd,wsd,uss,vss,wss,hdgsd,pitsd,rolsd,ipow,unoise,...
         vnoise,wnoise,usnr,vsnr,wsnr,pgp,serial_mc,serial_arg] = ...
         rodbload(infile,[columns,':MicrocatSN:SerialNumber']);
    elseif find(intersect(all_sn(proc),individual_sn))
        [YY,MM,DD,HH,t,p,u,v,w,...
         hdg,pit,rol,usd,vsd,wsd,uss,vss,wss,hdgsd,pitsd,rolsd,ipow,unoise,...
         vnoise,wnoise,usnr,vsnr,wsnr,pgp,serial_arg] = ...
         rodbload(infile,[columns,':SerialNumber']);
    end
    
    %------------------------------------------ 
    %----- cut off launching and recovery period
    %------------------------------------------
    disp('cut off launching and recovery period')
 
    jd  = julian(YY,MM,DD,HH);
    ii  = find(jd <= jd_e & jd >= jd_s );

    YY=YY(ii);MM=MM(ii);DD=DD(ii);HH=HH(ii);
    t=t(ii);p=p(ii);
    u=u(ii);v=v(ii);w=w(ii);

    if find(intersect(all_sn(proc),combo_sn))
        tc=tc(ii);c=c(ii);pc=pc(ii);
    end 
        
    hdg=hdg(ii);pit=pit(ii);rol=rol(ii);usd=usd(ii);vsd=vsd(ii);wsd=wsd(ii);
    uss=uss(ii);vss=vss(ii);wss=wss(ii);hdgsd=hdgsd(ii);pitsd=pitsd(ii);
    rolsd=rolsd(ii);ipow=ipow(ii);unoise=unoise(ii);vnoise=vnoise(ii);
    wnoise=wnoise(ii);usnr=usnr(ii);vsnr=vsnr(ii);wsnr=wsnr(ii);pgp=pgp(ii);

    jd  = jd(ii); 

    cycles     = length(ii);
    Start_Date = [YY(1) MM(1) DD(1)];
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

    if 0 % copied from preivous v003. surely this loop doesn't even start?

        for i = 1 : length(gap), 
          addt = [addt; [[1:gap(i)]*sr + jd(ii(i))]'];

        end 

        [jd,xx] = sort([jd; addt]);   % add time
        ngap    = length(addt);       % number of time gaps         
        gt      = gregorian(jd);
        YY=gt(:,1); MM=gt(:,2); DD=gt(:,3); 
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

    end % if 0

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
    
    if find(intersect(all_sn(proc),combo_sn))
        fort =['%4.4d %2.2d %2.2d  %6.4f    %7.4f %7.4f %6.2f %6.2f %7.4f  ',...
              ' %6.2f %6.2f %6.2f   %4.1f %4.1f %4.1f   %6.2f %6.2f %6.2f  ',...
              ' %3.3d %3.3d %3.3d   %4.1f %4.1f %4.1f   %4.2f   ',...
              ' %3.3d %3.3d %3.3d   %3.1f %3.1f %3.1f   %3.3d'];
        rodbsave(outfile,...
          'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time:MicrocatSN',...
           fort,...
          lat,lon,columns,Start_Date,Start_Time,all_sn(proc),mr,wd,indep,End_Date,End_Time,serial_mc,...
          [ YY,MM,DD,HH,t,tc,p,pc,c,u,v,w,...
            hdg,pit,rol,usd,vsd,wsd,uss,vss,wss,hdgsd,pitsd,rolsd,ipow,unoise,vnoise,wnoise,usnr,vsnr,wsnr,pgp]);
   
    elseif find(intersect(all_sn(proc),individual_sn))
        fort =['%4.4d %2.2d %2.2d  %6.4f    %7.4f %6.2f ',...
              ' %6.2f %6.2f %6.2f   %4.1f %4.1f %4.1f   %6.2f %6.2f %6.2f  ',...
              ' %3.3d %3.3d %3.3d   %4.1f %4.1f %4.1f   %4.2f   ',...
              ' %3.3d %3.3d %3.3d   %3.1f %3.1f %3.1f   %3.3d'];
        rodbsave(outfile,...
          'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
           fort,...
          lat,lon,columns,Start_Date,Start_Time,all_sn(proc),mr,wd,indep,End_Date,End_Time,...
          [ YY,MM,DD,HH,t,p,u,v,w,...
            hdg,pit,rol,usd,vsd,wsd,uss,vss,wss,hdgsd,pitsd,rolsd,ipow,unoise,vnoise,wnoise,usnr,vsnr,wsnr,pgp]);
   
    end
    
      
    

  %%%%%%%%%% Graphics %%%%%%%%%%%%%%%%

    jd0 = julian(-1,1,1,24);
    jd1 = julian(plot_interval(1,:))-jd0;
    jd2 = julian(plot_interval(2,:))-jd0; 

      
    sampling_rate = 1/median(diff(jd));

    if find(intersect(all_sn(proc),combo_sn))
        STR = ['Temperature [deg C]   ';
             'Pressure [dbar]       ';  
             'Conductivity [mS/cm]  ';
             'Zonal Velocity [cm/s] ';
             'Merid. Velocity [cm/s]'];
        STR2 = {'velocity stdev [cm/s]';
                'beam sig. str. [counts]';
                'beam noise [counts]';
                'beam SNR'};
        VAR1= ['t';'p';'c';'u';'v']; 
        VAR2= ['tc';'pc'];        
        VAR3= ['usd';'vsd';'wsd';'uss';'vss';'wss';'unoise';'vnoise';'wnoise';'usnr';'vsnr';'wsnr'];
        panels=5;
        panels2=4;
        
        figure(1);clf
        for sub = 1 : 5
          eval(['var1 = ',VAR1(sub),';'])
          if sub < 3 
            eval(['var2 = ',VAR2(sub,:),';'])
          else
            var2 = [];
          end  
        ok = plot_timeseries(jd,var1,var2,[],sampling_rate,STR(sub,:),sub,[jd1 jd2],'n',panels);
        end

        subplot(5,1,1)
           title(['ArgoCAT s/n: ',num2str(all_sn(proc)), ...
                 '; Target Depth: ',num2str(indep)])

           orient tall
        eval(['print -depsc ',outfile,'.eps']) 


        figure(2);clf

        for sub = 1 : 5
          eval(['var1 = ',VAR1(sub),';'])
          if sub < 3
            eval(['var2 = ',VAR2(sub,:),';'])
          else
            var2 = [];
          end  
        ok = plot_timeseries(jd,var1,var2,[],sampling_rate,STR(sub,:),sub,[jd1 jd2],'y',panels);
        end

        subplot(5,1,1)
          title(['ArgonCAT s/n: ',num2str(all_sn(proc)), ...
                 '; Target Depth: ',num2str(indep)])
          orient tall
        eval(['print -depsc ',outfile,'.filtered.eps']) 
        
        % plot of diagnostics info
        figure(3);clf
        m=0;
        for sub = 1:4
          eval(['var1 = ',VAR3{sub+m},';'])
          eval(['var2 = ',VAR3{sub+m+1},';'])
          eval(['var3 = ',VAR3{sub+m+2},';'])
          ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR2(sub,:),sub,[jd1 jd2],'n',panels);
          legend({'u','v','w'})
          hold on
          m=m+2;
        end

        subplot(4,1,1)
          title(['Argonaut s/n: ',num2str(all_sn(proc)), ...
                 '; Target Depth: ',num2str(indep)])

          orient tall
        eval(['print -depsc ',outfile,'_diagnostics.eps'])
        
    elseif find(intersect(all_sn(proc),individual_sn))
        STR = ['Temperature [deg C]   ';
             'Pressure [dbar]       ';  
             'Zonal Velocity [cm/s] ';
             'Merid. Velocity [cm/s]'];
        STR2 = {'velocity stdev [cm/s]';
               'beam sig. str. [counts]';
               'beam noise [counts]';
               'beam SNR'};
        VAR1= ['t';'p';'u';'v']; 
        var2=[];
        VAR3= {'usd';'vsd';'wsd';'uss';'vss';'wss';'unoise';'vnoise';'wnoise';'usnr';'vsnr';'wsnr'};
        panels=4;
        panels2=4;
        
        figure(1);clf
        for sub = 1 : 4
          eval(['var1 = ',VAR1(sub),';'])
          ok = plot_timeseries(jd,var1,var2,[],sampling_rate,STR(sub,:),sub,[jd1 jd2],'n',panels);
        end

        subplot(4,1,1)
           title(['Argonaut s/n: ',num2str(all_sn(proc)), ...
                 '; Target Depth: ',num2str(indep)])

           orient tall
        eval(['print -depsc ',outfile,'.eps']) 


        figure(2);clf

        for sub = 1 : 4
          eval(['var1 = ',VAR1(sub),';'])
          ok = plot_timeseries(jd,var1,var2,[],sampling_rate,STR(sub,:),sub,[jd1 jd2],'y',panels);
        end

        subplot(4,1,1)
          title(['Argonaut s/n: ',num2str(all_sn(proc)), ...
                 '; Target Depth: ',num2str(indep)])

          orient tall
        eval(['print -depsc ',outfile,'.filtered.eps'])
        
        % plot of diagnostics info
        figure(3);clf
        m=0;
        for sub = 1:4
          eval(['var1 = ',VAR3{sub+m},';'])
          eval(['var2 = ',VAR3{sub+m+1},';'])
          eval(['var3 = ',VAR3{sub+m+2},';'])
          ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR2(sub,:),sub,[jd1 jd2],'n',panels);
          legend({'u','v','w'})
          hold on
          m=m+2;
        end

        subplot(4,1,1)
          title(['Argonaut s/n: ',num2str(all_sn(proc)), ...
                 '; Target Depth: ',num2str(indep)])

          orient tall
        eval(['print -depsc ',outfile,'_diagnostics.eps'])
    end
  end  % for proc=1:length(combo_sn)+length(individual_sn) loop
  %keyboard
end % function
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function ok = plot_timeseries(jd,var1,var2,var3,sr,str,sub,jdlim,filt,panels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot time series


  jd0 = julian(-1,1,1,24);

  i1    = find(~isnan(var1) & var1~=0);
 
  if strcmp(filt,'y')
    var1  = auto_filt(var1(i1),sr,1/2,'low',4);
  else
    var1  = var1(i1);
  end 

  if ~isempty(var2)
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
         
  plot(jd(i1)-jd0,var1)
  hold on

  if ~isempty(var2)
    plot(jd(i2)-jd0,var2,'r')
  end
 
  if ~isempty(var3)
    plot(jd(i3)-jd0,var3,'g')
  end
  
  ylabel(str)
  grid on
  xlim([jdlim])
  datetick('x',12)
      
  ok=1;

