% MICROCAT_RAW2USE_003 is a script that performs stage2 processing
cruise = 'dy120';
operator  = 'lad';
moor = 'rteb1_05_2018';
% the start and end times of the time axis for plotting
plot_interval = [2018 09 07; 2020 14 10];
% ----------------- set path for data input and output --------------
%inpath   = [basedir,'moor/raw/',cruise,'/microcat/'];
inpath  = [basedir,'data/moor/proc/',moor,'/seaphox/'];
outpath  = [basedir 'data/moor/proc/' moor '/seaphox/'];
infofile = [basedir 'data/moor/proc/' moor '/' moor 'info.dat'];

[gash, operator]=system('whoami');  % This line will not work if run from a PC. May need to edit it out.

[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

vec=find(id==375); % Seaphox id number

sn=sn(vec);
z = z(vec);

fid_stat= fopen([outpath moor '_Seaphox_stage2.log'],'a');
fprintf(fid_stat,['Processing steps taken by ' mfilename ':\n']);
fprintf(fid_stat,'  1. eliminate launch and recovery period\n');
fprintf(fid_stat,'  2. resave data to rodb file\n');
fprintf(fid_stat,'\n Operated by:%s on %s\n',operator,datestr(clock)); 
fprintf(fid_stat,['        SeapHOx in Mooring ',moor,'\n\n\n']);

dummy    = -9999;

%-----------------------------------------
% --- preprocessing loop -------------------
% ----------------------------------------

jd_s  = datenum(s_d(1),s_d(2),s_d(3),s_t(1),s_t(2)/60,0);  % start time
jd_e  = datenum(e_d(1),e_d(2),e_d(3),e_t(1),e_t(2)/60,0);  % end time


for proc = 1 : length(vec)
    
columns = 'YY:MM:DD:HH:PH:PH_V:T:P:O:S:BATT1:BATT2';
indep  = z(proc);
infile  = [inpath, moor,'_',num2str(sn(proc)),'.raw'];
  
    if exist(infile)==0
        disp(['infile: ' infile ' does not exist.'])

    elseif exist(infile)   > 0 
        rodbfile= [moor,'_',num2str(sn(proc)),'.use']; 
        outfile = [outpath,rodbfile];
        fprintf(fid_stat,'Serialnumber %d \n',sn(proc));
        fprintf(fid_stat,'Infile %s \n',infile);
        fprintf(fid_stat,'Outfile %s \n',outfile);
    end
    
[YY,MM,DD,HH,ph,ph_v,t,p,o,s,batt1,batt2] = ...
        rodbload(infile,columns);

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
batt1=batt1(ii);batt2=batt2(ii);
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
fort =['%4.4d %2.2d %2.2d  %6.4f  %7.5f %10.8f  %6.4f %7.3f %6.3f  %6.4f %5.3f %5.3f'];

rodbsave(outfile,...
  'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
   fort,...
  lat,lon,columns,s_d,Start_Time,sn(proc),mr,wd,indep,End_Date,End_Time,...
  [YY,MM,DD,HH,ph,ph_v,t,p,o,s,batt1,batt2]);


%%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%
jd0 = julian(-1,1,1,24);
jd1 = julian(plot_interval(1,:))-jd0;
jd2 = julian(plot_interval(2,:))-jd0; 
sampling_rate = 1/median(diff(jd));

STR = ['pH                        ';
       'MC Temperature [deg C]    ';
       'MC Pressure [dbar]        ';
       'MC Oxygen [ml/l]          ';
       'MC Salinity               ';]
STR2 =['pH voltage [V]            ';
       'voltage [V]               ';]
VAR1= {'ph','t ','p ','o ','s '};
VAR2= {'ph_v','batt1','batt2'};
panels=5;
panels2=3;

figure(1);clf
for sub = 1 : 5
  eval(['var1 = ',VAR1{sub},';'])
  var2=[];
  var3=[];
  ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR(sub,:),sub,[jd1 jd2],'n',panels);          
end

subplot(5,1,1)
title(['SeapHOx s/n: ',num2str(sn(proc)), ...
     '; Target Depth: ',num2str(indep)])

orient tall
print(gcf,'-dpng',[outpath moor,'_',num2str(sn(proc))]) 


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
end

subplot(5,1,1)
title(['SeapHOx s/n: ',num2str(sn(proc)), ...
       '; Target Depth: ',num2str(indep)])
orient tall
print(gcf,'-dpng',[outpath moor,'_',num2str(sn(proc)) 'filtered']) 


% plot of diagnostics info
figure(3);clf
subplot(3,1,1)
eval(['var1 = ',VAR2{1},';'])
var2=[];
var3=[];
ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR2(1,:),1,[jd1 jd2],'n',panels2);
legend({'ph voltage'})

subplot(3,1,2)
eval(['var1 = ',VAR2{2},';'])
eval(['var2 = ',VAR2{3},';'])
var3=[];
ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR2(2,:),2,[jd1 jd2],'n',panels2);
legend({'main batt.','isolated batt.'})

suptitle(['SeapHOx s/n: ',num2str(sn(proc)), ...
         '; Target Depth: ',num2str(indep)])

orient tall
print(gcf,'-dpng',[outpath moor,'_',num2str(sn(proc)) 'ancilliary']) 

end % if exist(infile)==0
  
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
  if jdlim(2)-jdlim(1)>365/12;
    datetick('x',12)
  elseif jdlim(2)-jdlim(1)<1;
      datetick('x','HH:MM','keeplimits')
  else
      datetick('x','dd-mmm')
  end
      
  ok=1;
  end
