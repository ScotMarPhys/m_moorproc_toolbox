% function nortek_raw2use_cal_dip_01(moor,'procpath',procpath,'outpath',outpath)
%
% basic preprocessing for Nortek Aquadopp data
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

% 15/01/2007 - DR modified from argocat_raw2use_003

function nortek_raw2use_cal_dip_01(moor,varargin)

cruise ='d382';

if nargin==0
    help nortek_raw2use_cal_dip_01;
    return
end

% check for optional arguments
a=strmatch('procpath',varargin,'exact');
if a>0
    procpath=char(varargin(a+1));
else
   % procpath='/noc/ooc/rpdmoc/rapid/data/moor/proc/';
     %procpath='/Volumes/RB1201/rapid/data/moor/proc/';
     %procpath='/Users/hydrosea5/Desktop/RB1201/rapid/data/moor/proc';
     
     procpath='/noc/users/pstar/rpdmoc/rapid/data/moor/proc_calib/'; % path on D382 using oceanus workstation
end
a=strmatch('outpath',varargin,'exact');
if a>0
    outpath=char(varargin(a+1));
else
    outpath = [procpath '/' cruise '/cal_dip/nor/' moor '/'];
end
a=strmatch('plot_interval',varargin,'exact');
if a>0
    plot_interval=varargin(a+1);
else
    plot_interval=0;
end

[gash, operator]=system('whoami');  % This line will not work if run from a PC. May need to edit it out.

% -- set path for data input
%inpath  = [procpath '/' moor '/'];
inpath  = [procpath cruise '/cal_dip/nor/' moor '/'];

% --- get moring information from infofile 
infofile =[procpath cruise '/cal_dip/' moor 'info.dat'];
%infofile=['/Users/hydrosea5/Desktop/RB1201/rapid/data/moor/proc/wb1_8_201113/wb1_8_201113info.dat']
%infofile =[procpath '/' moor '/' moor 'info.dat']; % aurelie


[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

vec=find((id==368|id==370)); % Nortek id number

sn=sn(vec);
z = z(vec);

fid_stat= fopen([outpath moor '_Nortek_stage2.log'],'a');
fprintf(fid_stat,['Processing steps taken by ' mfilename ':\n']);
fprintf(fid_stat,'  1. eliminate lauching and recovery period\n');
fprintf(fid_stat,'  2. save data to rdb file\n');
fprintf(fid_stat,'\n Operated by:%s on %s\n',operator,datestr(clock)); 
fprintf(fid_stat,['        Norteks in Mooring ',moor,'\n\n\n']);

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

%all_z=[combo_z individual_z];
%all_sn=[combo_sn individual_sn];

for proc = 1 : length(vec)
    columns = 'YY:MM:DD:HH:T:P:U:V:W:HDG:PIT:ROL:USS:VSS:WSS:IPOW:CS:CD';
    indep  = z(proc);
    infile  = [inpath,'/',moor,'_',num2str(sn(proc)),'.raw'];
  
    if exist(infile)==0
        disp(['infile: ' infile ' does not exist.'])
      
    elseif exist(infile)   > 0 
        rodbfile= [moor,'_',num2str(sn(proc)),'.use']; 
        outfile = [outpath,rodbfile];
        fprintf(fid_stat,'Serialnumber %d \n',sn(proc));
        fprintf(fid_stat,'Infile %s \n',infile);
        fprintf(fid_stat,'Outfile %s \n',outfile);

        [YY,MM,DD,HH,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = ...
            rodbload(infile,[columns]);

        %------------------------------------------ 
        %----- cut off launching and recovery period
        %------------------------------------------
        disp('cut off launching and recovery period')

        jd  = julian(YY,MM,DD,HH);
        %gregorian(jd);
%         jd_e
%         gregorian(jd_e)
%         jd_s
%         gregorian(jd_s)
        
        ii  = find(jd <= jd_e & jd >= jd_s );

        YY=YY(ii);MM=MM(ii);DD=DD(ii);HH=HH(ii);
        t=t(ii);p=p(ii);
        u=u(ii);v=v(ii);w=w(ii);
        hdg=hdg(ii);pit=pit(ii);rol=rol(ii);
        uss=uss(ii);vss=vss(ii);wss=wss(ii);
        ipow=ipow(ii); cs=cs(ii); cd=cd(ii);
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
        fort =['%4.4d %2.2d %2.2d  %6.4f  %4.2f %7.3f  %4.1f %4.1f %4.1f  '...
            '%4.1f %4.1f %4.1f  %2d %2d %2d  %2.1f  %4.1f %5.2f'];
        
        rodbsave(outfile,...
          'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
           fort,...
          lat,lon,columns,s_d,Start_Time,sn(proc),mr,wd,indep,End_Date,End_Time,...
          [ YY,MM,DD,HH,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd]);
      
      
      % aurelie: i have replace Start_Date by s_d and it works !!!!
        
        %%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%
%who
        jd0 = julian(-1,1,1,24);
        jd1 = julian(plot_interval(1,:))-jd0;
        jd2 = julian(plot_interval(2,:))-jd0; 
        sampling_rate = 1/median(diff(jd));

        STR = ['Temperature [deg C]   ';
               'Pressure [dbar]       ';  
               'Zonal Velocity [cm/s] ';
               'Merid. Velocity [cm/s]';
               'Vert. Velocity [cm/s] ']; 
        VAR1= ['t';'p';'u';'v';'w']; 
        panels=5;
        figure(1);clf
        for sub = 1 : 5
          eval(['var1 = ',VAR1(sub),';'])
          var2=[];
          ok = plot_timeseries(jd,var1,var2,sampling_rate,STR(sub,:),sub,[jd1 jd2],'n',panels);
        end

        subplot(5,1,1)
        title(['Nortek s/n: ',num2str(sn(proc)), ...
             '; Target Depth: ',num2str(indep)])

        orient tall
        eval(['print -depsc ',outfile,'.eps']) 


        figure(2);clf

        for sub = 1 : 5
          eval(['var1 = ',VAR1(sub),';'])
          ok = plot_timeseries(jd,var1,var2,sampling_rate,STR(sub,:),sub,[jd1 jd2],'y',panels);
        end

        subplot(5,1,1)
        title(['Nortek s/n: ',num2str(sn(proc)), ...
               '; Target Depth: ',num2str(indep)])
        orient tall
        eval(['print -depsc ',outfile,'.filtered.eps']) 
        
        end % if exist(infile)==0
  end  % for proc=1:length(combo_sn)+length(individual_sn) loop
end % function
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function ok = plot_timeseries(jd,var1,var2,sr,str,sub,jdlim,filt,panels)
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

  subplot(panels,1,sub);
         
  plot(jd(i1)-jd0,var1)
  hold on

  if ~isempty(var2)
    plot(jd(i2)-jd0,var2,'r')
  end
 
  ylabel(str)
  grid on
  xlim([jdlim])
  datetick('x',12)
      
  ok=1;

  end
