% 1) dedrift microcats using post- and pre-cruise  ctd calibration 
% 2) set suspicious data to dummy 
%
% prior to this programme, insitu_cal_osnap.m has to be run for pre- and post-cruise calibration
% calibration offset then have to be entered in 'microcat_cond.csv' /
% 'microcat_temp.csv' / microcat_pres.csv 
%
%
%
% kanzow, 23.08.05 
% edited by Loic Houpert 21/10/2015
% edited by lewis Drysdale Dec 2020
% edited by Loic Houpert Jan 2021

close all
clearvars  -except pathosnap pathgit
warning off

% path of the mooring data define in the startup file under osnap/
%moor = 'rtwb1_05_2018';
moor = 'ib5_01_2018';

%=========================================================================
% Apply calibration coefficients to series, removes bad data. If required, applies
% constant offsets, and conductivity pressure correction
p_applycal.operator  = 'LH';
p_applycal.mooring  = moor;   
p_applycal.sensortyp = 'microcat';   % arg / microcat / rbr / idr
p_applycal.delim = ',';

% input directories & files 
p_applycal.mooring_dir         = [pathosnap '/data/moor/proc/'];
p_applycal.mooring_outdir      = [pathosnap '/data/moor/proc/'];
p_applycal.coef_dir            = [pathgit '/data/moor/cal_coef/']; 
p_applycal.external_ctd_dir    = [pathosnap '/cruise_data/'];
p_applycal.ctd_ref_cruises     = {''};%{'pe400'}; %{'kn221-02';'pe399'}; % references cruises for the QC
p_applycal.distrangectd        = 100e3; % distance of the reference ctd from the mooring
p_applycal.strformat.mctemptxt = repmat('%s',1,32);
p_applycal.strformat.mctempnum = ['%f%f%s%s%f%s%s%f%f%s%s' repmat('%f%s%s',1,7)];
p_applycal.strformat.mcsaltxt  = repmat('%s',1,32);
p_applycal.strformat.mcsalnum  = ['%f%f%s%s%f%s%s%f%f%s%s' repmat('%f%s%s',1,7)];
p_applycal.strformat.mcprestxt = repmat('%s',1,65);
p_applycal.strformat.mcpresnum = ['%f' repmat('%s%f%s%s',1,16)];

loclegend = 'north';

%-------------------------------------------------------------------------
% Before running the script, check that the ctd data are present in 
% external_ctd_dir . If not, script to convert mstar ctd format to rodb 
% format are in users/loh/mcatpostcruisecalib/ctd_matfile_rodb/
% ---------------------------------------------------------------------------

% ---------------------------------------------------------------------------
operator  = p_applycal.operator; 
mooring   = p_applycal.mooring ;
sensortyp = p_applycal.sensortyp;   % arg / microcat / rbr / idr
strformat = p_applycal.strformat;

delim = p_applycal.delim;

distrangectd = p_applycal.distrangectd;
mooring_dir         = p_applycal.mooring_dir;
mooring_outdir      = p_applycal.mooring_outdir;
coef_dir            = p_applycal.coef_dir;
external_ctd = [];
external_ctd_file = [];

for ijk=1:length(p_applycal.ctd_ref_cruises)
    external_ctd        = [external_ctd ; ...
            [p_applycal.external_ctd_dir p_applycal.ctd_ref_cruises{ijk} '/ctdref/']];
    external_ctd_file   = [external_ctd_file ; ...
            [p_applycal.ctd_ref_cruises{ijk} '_pos.mat']];
end

% extract calibration cruise indices from database (.xls spreadsheet)
fid_cruise = fopen([coef_dir,'microcat_calib_cruise.csv']);
cruise     = textscan(fid_cruise,'%s%d%d','delimiter',delim,'HeaderLines',2);
fclose(fid_cruise)

moorI                = find((strcmp(cruise{1},mooring))>0);
cruise               = [cruise{2} cruise{3}];
cruiseI              = cruise(moorI,:);

%  extract deployment can recovery cruise names from database

[num,txt,raw]=xlsread([coef_dir,'cruise_id.xls']);
deplcr = txt(find(num==cruiseI(1))+1);
reccr  = txt(find(num==cruiseI(2))+1);

% ------ output  data ---- 

%%ext      = '.microcat';
dum      = -9999; 
pref     = 0;
t90_68   = 1.00024;  % convert its90 to its68 for cond. to sal.

% ------ conversion ------------
cols      = 'YY:MM:DD:HH:T:C'; % column info for rodb header
colsp     = 'YY:MM:DD:HH:T:C:P'; % column info for rodb header (mc with pressure sensor)
colsarg   = 'YY:MM:DD:HH:T:TCAT:P:PCAT:C:U:V:W';
fort      = '%4.4d  %2.2d  %2.2d  %7.5f   %6.4f  %6.4f'; %data output format
fortp     = '%4.4d  %2.2d  %2.2d  %7.5f   %6.4f  %6.4f  %5.1f'; %data output format(mc with pressure sensor)  

% ------ sensor specific settings ------

if strcmp('microcat',sensortyp)
  typ_id = [333 337];
elseif strcmp('arg',sensortyp)
  typ_id = [366];  
elseif strcmp('rbr',sensortyp)
  typ_id = [330];
elseif strcmp('idr',sensortyp)
    typ_id = [339];
end
ext = ['.',sensortyp];

% ---- microcat raw data ----
head    = ['Mooring:SerialNumber:WaterDepth:InstrDepth:Start_Date:Start_Time:End_Date:End_Time:Latitude:Longitude:Columns'];

% ---- load calib offsets (pre and post cruise) ----
% Conductivity / Temperature / Depth

if strcmp('microcat',sensortyp) | strcmp('idr',sensortyp) | strcmp('rbr',sensortyp)
 
  fidt  = fopen([coef_dir,'microcat_temp.csv'],'r');
  ttext = textscan(fidt,strformat.mctemptxt,'delimiter',delim);
  fseek(fidt,0,'bof');
  tnum  = textscan(fidt,strformat.mctempnum,'delimiter',delim,'HeaderLines',5);
  
  fidc  = fopen([coef_dir,'microcat_cond.csv'],'r');
  ctext = textscan(fidc,strformat.mcsaltxt ,'delimiter',delim);
  fseek(fidc,0,'bof');
  cnum  = textscan(fidc,strformat.mcsalnum,'delimiter',delim,'HeaderLines',5);


  fidp  = fopen([coef_dir,'microcat_pres.csv'],'r');
  ptext = textscan(fidp,strformat.mcprestxt,'delimiter',delim);
  fseek(fidp,0,'bof');
  pnum  = textscan(fidp,strformat.mcpresnum,'delimiter',delim,'HeaderLines',5);

elseif strcmp('arg',sensortyp)
    
  [cnum,ctext] = xlsread([coef_dir,'microcat_cond.xls']);
  [tnum,ttext] = xlsread([coef_dir,'microcat_temp.xls']);
  [pnum,ptext] = xlsread([coef_dir,'microcat_pres.xls']);
  [targnum,targtext] = xlsread([coef_dir,'argonaut_temp.xls']);
  [pargnum,pargtext] = xlsread([coef_dir,'argonaut_pres.xls']);
  targin             =   targnum(:,1);  % Instruments
  pargin             =   pargnum(:,1);  % Instruments

end

% ---- pick instruments ---------------

cin    =   cnum{1};   % Instruments
tin    =   tnum{1};   % Instruments
pin   = pnum{1};
plen  = length(ptext);
clen  = length(ctext);
tlen  = length(ttext);

% ------ pre cruise calibration ----------------

if cruiseI(1) < 0
    preC(1:length(cin),1) = NaN;
    preT(1:length(tin),1) = NaN;
    preP(1:length(pin),1) = NaN;
    preCcomment(1:length(cin),1) = {' '};
    preTcomment(1:length(tin),1) = {' '};
    prePcomment(1:length(pin),1) = {' '};  
else

    
% ------ Codcutivity ------
  for i = 1 :clen
    cvalI(i)=strcmp(ctext{i}{2},num2str(cruiseI(1)));
  end
  ccol = find(cvalI == 1);
% ------ Temperature ------
  for i = 1 :tlen
    tvalI(i)=strcmp(ttext{i}{2},num2str(cruiseI(1)));
  end
  tcol = find(tvalI == 1);
% ------ Pressure ------
  for i = 1 :plen
    pvalI(i)=strcmp(ptext{i}{2},num2str(cruiseI(1)));
  end
  pcol = find(pvalI == 1);
 
  preC  =  cnum{ccol}; %pre cruise offsets
  preCcomment =  cnum{ccol+2};
  
  preT  =  tnum{tcol}; %pre cruise offsets
  preTcomment =  tnum{tcol+2};
  
    if ~isempty(find(strcmp(ptext{pcol(1)-1},mooring)>0))
        preP  = pnum{pcol(1)};
        prePcomment =  pnum{pcol(1)+2};
    elseif ~isempty(find(strcmp(ptext{pcol(2)-1},mooring)>0))
        preP  = pnum{pcol(2)}; 
        prePcomment =  pnum{pcol(2)+2};
    end
end

% ----- post cruise calibration -----------------

if cruiseI(2) < 0

  postC(1:length(cin),1) = NaN;
  postT(1:length(tin),1) = NaN;
  postP(1:length(pin),1) = NaN;
    postCcomment(1:length(cin),1) = {' '};
    postTcomment(1:length(tin),1) = {' '};
    postPcomment(1:length(pin),1) = {' '};

else
    
  for i = 1 :clen
    cvalI(i)=strcmp(ctext{i}{2},num2str(cruiseI(2)));
  end
  ccol = find(cvalI == 1);
  
  for i = 1 :tlen
    tvalI(i)=strcmp(ttext{i}{2},num2str(cruiseI(2)));
  end
  tcol = find(tvalI == 1);
  
  for i = 1 :plen
     pvalI(i)=strcmp(ptext{i}{2},num2str(cruiseI(2)));
  end
  pcol = find(pvalI == 1);

  postC  = cnum{ccol}; %post cruise offsets
  postCcomment =  cnum{ccol+2};  
  
  postT  = tnum{tcol}; %post cruise offsets
  postTcomment =  tnum{tcol+2};  
     
  if ~isempty(find(strcmp(ptext{pcol(1)-1},mooring)>0))
  
    postP  = pnum{pcol(1)};
    postPcomment =  pnum{pcol(1)+2};
  
  elseif ~isempty(find(strcmp(ptext{pcol(2)-1},mooring)>0))
      postP  = pnum{pcol(2)}; 
      postPcomment =  pnum{pcol(2)+2};
  end
end

% General trend: to be applied when pre - or post cruise calibration is 
% missing for individual instruments

nnan   = find(~isnan(postC) & ~isnan(preC));%lines 1:5 are header
val    = find(abs(postC(nnan))<0.02 & abs(preC(nnan)) < 0.02);
Ctrend = -mean(preC(nnan(val)))+mean(postC(nnan(val)));  %allg. Vergleich pre - post cruise 

nnan   = find(~isnan(postT) & ~isnan(preT));%lines 1:5 are header
val    = find(abs(postT(nnan))<0.02 & abs(preT(nnan)) < 0.02);
Ttrend = -mean(preT(nnan(val)))+mean(postT(nnan(val)));

if isempty(val) 
  Ctrend = 0; Ttrend = 0;   
end

nnan   = find(~isnan(postP) & ~isnan(preP));%lines 1:5 are header
val    = find(abs(postP(nnan))<100 & abs(preP(nnan)) < 100);
Ptrend = -mean(preP(nnan(val)))+mean(postP(nnan(val)));

if isempty(val) 
  Ptrend = 0;    
end

% ---- 1) apply dedrift ---------

[typ,dep,serial] = rodbload([mooring_dir,mooring,'/',mooring,'info.dat'],...
                        ['instrument:z:serialnumber']);

mcI     = find(typ >= typ_id(1) & typ <= typ_id(end));
dep     = dep(mcI);
serial  = serial(mcI);
typ     = typ(mcI);

for mc = 1: length(serial)
  close all
  mcfile_out = [mooring_outdir,mooring,'/',sensortyp,'/',mooring,'_',sprintf('%3.3d',mc),ext]; 
  mcfig_out = [mooring_outdir,mooring,'/',sensortyp,'/',mooring,'_',sprintf('%3.3d',mc)]; 
  if exist(mcfile_out) == 2 
    over = input(['File ' [sprintf('%3.3d',mc) ext] ' for instrument ' num2str(serial(mc)) ' exists already: do you want to overwrite it? y/n '],'s');  
    if ~strcmp(over,'y')
         continue
    end
    [succ1,mess1] =copyfile(mcfile_out,[mcfile_out,date]);
    [succ2,mess2] =copyfile([mcfile_out,'.txt'],[mcfile_out,'.txt',date]);
  end  
     
  if exist('diag','var')
      fclose(diag);
      diag=fopen('CTD_near.txt','w+');
  end
  skipC  = 1;
  skipT  = 1;
  skipP  = 1;
  
  mcfile = [mooring_dir,mooring,'/',sensortyp,'/',mooring,'_',sprintf('%4.4d',serial(mc)),'.use'];

  if exist(mcfile)==0
      disp(['Problem! the file :\n' mcfile '\n didn t exist'])
      continue
  end
  
  % read .USE file
  if strcmp(sensortyp,'microcat') | strcmp(sensortyp,'rbr') | strcmp(sensortyp,'idr')  
    [yy,mm,dd,hh,t,c,p]= rodbload(mcfile,colsp);
  elseif strcmp(sensortyp,'arg')
    [yy,mm,dd,hh,t,tcat,p,pcat,c,u,v,w]= rodbload(mcfile,colsarg);
  end      
      
  [moo,serial0,wd,id,sd,st,ed,et,lt,ln,cl] = rodbload(mcfile,head);

  jd       = julian(yy,mm,dd,hh);
  jd0      = jd - jd(1);  
  xli      = [0 max(jd0)];

  spikeT   = find(t<-900);
  spikeC   = find(c<-900); 
  spikeP   = find(p<-0); 
  
  cinI     = find(cin == serial0);
  tinI     = find(tin == serial0); 
  p_exist  = 'n';
  sp       = 3;
  if ~isempty(find(~isnan(p)))
    p_exist = 'y';
    pinI     = find(pin == serial0); 
    sp       = 4;   
  end

%% CONDUCTIVITY
  warn = 0;
  if isempty(cinI)
     disp(['no conduct. calibration entry found for serial number ',num2str(serial0)])
     skipC    =  input('skip (0) or save data without calibration(1) ');
     corrC    =  0;
  else
    cal1     =  preC(cinI);
    cal1comment = preCcomment(cinI);
    cal2     =  postC(cinI);       
    cal2comment = postCcomment(cinI);
    
    
    disp('   ')
    disp(['Conductivity pre-cruise calibration: ',num2str(cal1)])
    disp(['Comment: ',cal1comment{1}])
    disp(['Conductivity post-cruise calibration: ',num2str(cal2)])
    disp(['Comment: ',cal2comment{1}])
    disp('  ')  
    keep = input('Type 0, 1 or 2 if you want to use both, the 1st or the 2nd coefficient ');
    disp('  ')
    if keep == 1
      postC(cinI) = NaN;
      disp('Only using the pre cruise calibration coefficient')
    elseif keep == 2
      preC(cinI) = NaN;
      disp('Only using the post cruise calibration coefficient')
    elseif keep == 0
      disp('Using both calibration coefficients')
    else
      disp(['Choice ',num2str(keep),' not defined'])
      break
    end
    
    % clibration offsets for conductivity
    cal1     =  preC(cinI);
    cal2     =  postC(cinI);       
      
    % IF one of the offsets is a nan, we can use a general trend
    if isnan(cal1)
      applytrendcg = input('SHOULD GENERAL TREND FOR CONDUCTIVITIES BE APPLIED? [y/n]','s');
      if strcmp(applytrendcg,'y')
        cal1   =  cal2-Ctrend;
      else
        cal1   =  cal2-Ctrend*0;
      end    
      warn = 1;  
    end
    
    if isnan(cal2)
      applytrendcg = input('SHOULD GENERAL TREND FOR CONDUCTIVITIES BE APPLIED? [y/n]','s');  
      if strcmp(applytrendcg,'y')
        cal2   =  cal1+Ctrend;
      else
        cal2   =  cal1+Ctrend*0;
      end  
      warn = 2;
    end  
    
    % Ask if want to use average 
    if keep == 0
        applytrendc = input('SHOULD AVERAGE OF CONDUCTIVITY COEFFICIENTS BE APPLIED? [y/n]','s');
      if strcmp(applytrendc,'y')
        cal1   =  (cal1+cal2)/2;
        cal2   =  cal1;
      else
        cal1   =  cal1;
        cal2   =  cal2;
      end
      warn = 3;
    end
      
    % calculate the correction to be applied to the data
    corrC    =  cal1 + (cal2 - cal1)/max(jd0)*jd0; 
    
    if isempty(find(~isnan(corrC)))
      disp(['no conduct. calibration found for serial number ',num2str(serial0)])
      skipC    =  input('skip (0) or save cond. data without calibration(1) ')
      corrC    =  0; 
    end    
  end  
  
%% TEMPERATURE    
  if isempty(tinI)
     disp(['no temp. calibration entry found for serial number ',num2str(serial0)])
     skipT    =  input('skip (0) or save data without calibration(1) ');
     corrT    =  0;
  else
    cal1comment = preTcomment(tinI); 
    cal2comment = postTcomment(tinI);
    cal1        =  preT(tinI);
    cal2        =  postT(tinI); 
    disp('   ')
    disp(['Temperature pre-cruise calibration: ',num2str(cal1)])
    disp(['Comment: ',cal1comment{1}])
    disp(['Temperature post-cruise calibration: ',num2str(cal2)])
    disp(['Comment: ',cal2comment{1}])
    disp('  ')  
    keep = input('Type 0, 1 or 2 if you want to use both, the 1st or the 2nd coefficient ');
    disp('  ')
    if keep == 1
      postT(tinI) = NaN;
      disp('Only using the pre cruise calibration coefficient')
    elseif keep == 2
      preT(tinI) = NaN;
      disp('Only using the post cruise calibration coefficient')
    elseif keep == 0
      disp('Using both calibration coefficients')
    else
      disp(['Choice ',num2str(keep),' not defined'])
      break
    end
     
    cal1     =  preT(tinI);
    cal2     =  postT(tinI); 
    if isnan(cal1)
      applytrendtg = input('SHOULD GENERAL TREND FOR TEMPERATURES BE APPLIED? [y/n]','s');
      if strcmp(applytrendtg,'y')
        cal1   =  cal2-Ttrend;
      else
        cal1   =  cal2-Ttrend*0;
      end    
      warn = 1;
    end
  
    if isnan(cal2)
      applytrendtg = input('SHOULD GENERAL TREND FOR TEMPERATURES BE APPLIED? [y/n]','s');  
      if strcmp(applytrendtg,'y')
        cal2   =  cal1+Ttrend;
      else
        cal2   =  cal1+Ttrend*0;
      end  
      warn = 2;
    end  
    
    if keep == 0
        applytrendt = input('SHOULD AVERAGE OF TEMPERATURE COEFFICIENTS BE APPLIED? [y/n]','s');
      if strcmp(applytrendt,'y')
        cal1   =  (cal1+cal2)/2;
        cal2   =  cal1;
      else
        cal1   =  cal1;
        cal2   =  cal2;
      end
      warn = 3;
    end

    corrT    =  cal1 + (cal2 - cal1)/max(jd0)*jd0; 
    if isempty(find(~isnan(corrT)))
      disp(['no temp. calibration found for serial number ',num2str(serial0)])
      skipT    =  input('skip (0) or save temp. data without calibration(1) ')
      corrT    =  0; 
    end    
  end

  %% PRESSURE
if strcmp(p_exist,'y')
    if isempty(pinI)
        disp(['no pressure calibration found for serial number ',num2str(serial0)])
        skipP    =  input('skip (0) or save data without calibration(1) ');
        corrP    =  0;
    else   
      cal1comment = prePcomment(pinI); 
      cal2comment = postPcomment(pinI);
      cal1        =  preP(pinI);
      cal2        =  postP(pinI); 
      
      disp(' ')
      disp(['Pressure precruise calibration: ',num2str(cal1)])
      disp(['Comment: ',cal1comment{1}])
      disp(['Pressure postcruise calibration: ',num2str(cal2)])
      disp(['Comment: ',cal2comment{1}])
      disp(' ')
      keep = input('Type 0, 1 or 2 if you want to use both, the 1st or the 2nd coefficient ');
      disp(' ')
      if keep == 1
          postP(pinI) = NaN;
          disp('Only using the pre cruise calibration coefficient')
      elseif keep == 2
          preP(pinI) = NaN;
          disp('Only using the post cruise calibration coefficient')
      elseif keep == 0
          disp('Using both calibration coefficients')
      else
          disp(['Choice ',num2str(keep),' not defined'])
          break
      end
      disp(' ')
      cal1     =  preP(pinI);
      cal2     =  postP(pinI); 
      
    if isnan(cal1)
      applytrendpg = input('SHOULD GENERAL TREND FOR PRESSURES BE APPLIED? [y/n]','s');
      if strcmp(applytrendpg,'y')
        cal1   =  cal2-Ptrend;
      else
        cal1   =  cal2-Ptrend*0;
      end    
      warn = 1;  
    end
    
    if isnan(cal2)
      applytrendpg = input('SHOULD GENERAL TREND FOR PRESSURES BE APPLIED? [y/n]','s');  
      if strcmp(applytrendpg,'y')
        cal2   =  cal1+Ptrend;
      else
        cal2   =  cal1+Ptrend*0;
      end  
      warn = 2;
    end  

    if keep == 0
        applytrendp = input('SHOULD AVERAGE OF PRESSURE COEFFICIENTS BE APPLIED? [y/n]','s');
      if strcmp(applytrendp,'y')
        cal1   =  (cal1+cal2)/2;
        cal2   =  cal1;
      else
        cal1   =  cal1;
        cal2   =  cal2;
      end
      warn = 3;
    end
      
    corrP    =  cal1 + (cal2 - cal1)/max(jd0)*jd0; 
        if isempty(find(~isnan(corrP)))
            disp(['no pressure calibration found for serial number ',num2str(serial0)])
            skipP    =  input('skip (0) or save pressure data without calibration(1) ');
            corrP    =  0; 
        end    

    end 
   
    pn       = p - corrP;
    
end

%% apply offsets 

  tn       = t - corrT;      
  cn       = c - corrC;

  tn(spikeT) = dum;
  cn(spikeC) = dum;
  pn(spikeP) = dum;
  
  figure(1)
  subplot(1,sp,1)
  hold off
  ii = find(tn>dum);
  plot(jd0(ii),t(ii),'k')
  hold on
  plot(jd0(ii),tn(ii),'r')
  grid on;  title(['temperature sn',num2str(serial0)]);xlim(xli);
  
  subplot(1,sp,2)
  hold off
  ii = find(cn>dum);
  plot(jd0(ii),c(ii),'k')
  hold on
  plot(jd0(ii),cn(ii),'r')
  grid on;  title('C raw=black plus offset=red');xlim(xli);

  if strcmp(p_exist,'y')
    subplot(1,sp,sp-1)
    hold off
    ii = find(pn>dum);
    plot(jd0(ii),p(ii),'k')
    hold on
    plot(jd0(ii),pn(ii),'r')
    grid on; title('pressure'); xlim(xli);
  end

  subplot(1,sp,sp)
  hold off
  if warn == 0
    plot(jd0([1 end]),corrT([1 end])*1000,'o-')
    hold on
    plot(jd0([1 end]),corrC([1 end])*1000,'ro-')
    legend('T drift','C drift','location',loclegend); xlim(xli);
    if strcmp(p_exist,'y')
      plot(jd0([1 end]),corrP([1 end]),'go-')
      legend('T drift','C drift','P drift','location',loclegend)
    end
  elseif warn ==1
    pl1=plot(jd0([1 end]),corrT([1 end])*1000,'--',jd0(end),corrT(end)*1000,'o');
    hold on    
    pl2=plot(jd0([1 end]),corrC([1 end])*1000,'r--',jd0(end),corrC(end)*1000,'ro');
    legend([pl1(1) pl2(1)],'T drift','C drift','location',loclegend)
    if strcmp(p_exist,'y')
      pl3=plot(jd0([1 end]),corrP([1 end]),'g-',jd0(end),corrP(end),'go');
      legend([pl1(1) pl2(1) pl3(1)],'T drift','C drift','P drift','location',loclegend)
    end
  elseif warn ==2
    pl1=plot(jd0([1 end]),corrT([1 end])*1000,'--',jd0(1),corrT(1)*1000,'o');
    hold on    
    pl2=plot(jd0([1 end]),corrC([1 end])*1000,'r--',jd0(1),corrC(1)*1000,'ro');
    legend([pl1(1) pl2(1)],'T drift','C drift','location',loclegend)
    if strcmp(p_exist,'y')
      pl3=plot(jd0([1 end]),corrP([1 end]),'g-',jd0(1),corrP(1),'go');
      legend([pl1(1) pl2(1) pl3(1)],'T drift','C drift','P drift','location',loclegend)
    end
  elseif warn ==3
      plot(jd0([1 end]),corrT([1 end])*1000,'o-')
    hold on
    plot(jd0([1 end]),corrC([1 end])*1000,'ro-')
    legend('T drift','C drift','location',loclegend); xlim(xli);
    if strcmp(p_exist,'y')
      plot(jd0([1 end]),corrP([1 end]),'go-')
      legend('T drift','C drift','P drift','location',loclegend)
    end
  end 
  grid on 
 
  % Conductivity at (35,15,0)
  c3515   = sw_c3515;
  
  ii = find(cn>dum & pn>dum); 
  
  % Salinity from cndr, T, P
  sn = sw_salt(cn(ii)/c3515,tn(ii)*t90_68,pn(ii));
 
  figure(6);clf; hold on
  
  plot(jd(ii)-jd(1),sn)
  ylabel('Salinity')
  title(['Salinity at ' num2str(nanmean(p),'%04.0f') 'm'])
  
  figure(6)
  
  print([mcfig_out,'_salinity', '.png'],'-dpng');
  sss=sw_salt(c/c3515,t*t90_68,p);
  figure(8);clf; hold on
  plot(sss,t,'.')
  xlabel('Salinity')
  ylabel('Temperature')
  
%%%%%%%%%%%%%% --- 2) conductivity pressure correction

acc_cpcor = input('Is a pressure dependant correction is required for the conductivity (see Fig. 1)? y/n ','s');

if strcmp(acc_cpcor,'y')
  pref = input('Insert reference pressure [dbar] ');
  invalcpcor = input('Insert time interval where correction is required, e.g. [3 56] denotes days 3-56 ');
  dumI = find(cn == dum);
  
  for i = 1 : length(invalcpcor)/2  
    corI = find(jd0>= invalcpcor(i*2-1) & jd0 <= invalcpcor(i*2));    
    cn(corI) = mc_concorr(cn(corI),pref,pn(corI));
  end
  cn(dumI) = dum;
  figure(6)
    sn = sw_salt(cn(ii)/c3515,tn(ii)*t90_68,pref*ones(length(ii),1));
    plot(jd(ii)-jd(1),sn,'r')
    grid on
    legend('raw','P corr.')

elseif strcmp(acc_cpcor,'n')
  disp('No pressure conductivity correction applied')
  
elseif ~strcmp(acc_cpcor,'y') & ~strcmp(acc_cpcor,'n')
  disp('Answer not defined - no pressure conductivity correction applied')  
  acc_cpcor = 'n';
  
end

%%%%%%%%% 3)  pressure drift removal ---------------

acc_drift = input('Does the pressure record require drift removal (see Fig. 1)?  y/n ','s');

if strcmp(acc_drift,'y')
  [coef,fit]= exp_lin_fit2(jd,pn,[1 1 1 1]);
  
  figure(3);clf;hold on

  plot(jd0,pn)

  plot(jd0,fit,'r','Linewidth',2)
  grid on
  
  if warn == 1
    pp = pn - fit + fit(end);
  elseif warn == 2
    pp = pn - fit + fit(1);
  else 
    pp = pn - fit + mean(fit);
  end
  
  so = sw_salt(cn(ii)/c3515,tn(ii)*t90_68,pn(ii));
  sn = sw_salt(cn(ii)/c3515,tn(ii)*t90_68,pp(ii));
  deltaS = std(sn)-std(so)
    fprintf(1,'\n\n Salinity standard deviation should decrease with dedrifted pressure \n')
  if deltaS > 0
    fprintf(1,' Salinity standard deviation has INCREASED instead  by %4.4f \n\n',deltaS)
  elseif deltaS < 0
    fprintf(1,' WARNING: Salinity standard deviation has DECREASED by %4.4f \n\n',-deltaS)
  end
  plot(jd0,pp,'g')
  
  figure(111);clf;hold on
    plot(so,tn(ii),'b.')
    plot(sn,tn(ii),'r.')
    grid on
    xlabel('S')
    ylabel('T')
    drawnow
    legend('P raw','P dedrifted')
    
  acc = input('Do you accept the the drift removal y/n ','s');

  if strcmp(acc,'y')
    disp('drift removal accepted')
    
    figure(1)
    subplot(1,sp,sp-1)
    plot(jd0(ii),pn(ii),'w')
    plot(jd0(ii),pp(ii),'r')
    pn =pp;
    figure(6)
    plot(jd(ii)-jd(1),sn,'r')
    legend('raw','P corr')
    pn(spikeP)=dum;
    
  else
    disp('drift removal discarded')
    acc_drift = 'n';
  end    
elseif strcmp(acc_drift,'n')
  disp('No drift correction applied')
  
elseif ~strcmp(acc_drift,'y') & ~strcmp(acc_drift,'n')
  disp('Answer not defined - no drift correction carried applied')  
  acc_drift = 'n';
end    

  %%%%%%%----- 4)  set suspicious data to dummies -------%%%%%%%%%%
  
  invalT = 0; invalC = 0; invalP = 0;
  
  inval= input('Should temperatures in a certain time interval be set to dummies y/n \n','s');

  if inval == 'y' % temp dummies loop
    invalT= input('Insert interval, e.g. [2 4 6 9] denotes days 2-4 and days 6-9 \n');      
    while 1
      if isodd(length(invalT))
        invalT= input('Repeat entry, vector must have even number of elements \n');  
      else
        break
      end
    end 

    for i = 1 : length(invalT)/2,  
      dumI = find(jd0>= invalT(i*2-1) & jd0 <= invalT(i*2));    
      tn(dumI) = dum;
    end
    figure(1)
    subplot(1,sp,1)
    hold off
    ii = find(tn>dum);
    plot(jd0(ii),t(ii),'k')
    hold on
    plot(jd0(ii),tn(ii),'b')
    title('temperature');xlim(xli);
  end %  temp. dummies loop 

  inval= input('Should conductivities in a certain time interval be set to dummies y/n \n','s');
  
  if inval == 'y' % cond dummies loop

    invalC= input('Enter interval, e.g. [2 4 6 9] denotes days 2-4 and days 6-9 \n');      
    while 1
      if isodd(length(invalC))
        invalC= input('Repeat entry, vector must have even number of elements \n')  
      else
        break
      end
    end 

    for i = 1 : length(invalC)/2,  
      dumI = find(jd0>= invalC(i*2-1) & jd0 <= invalC(i*2));    
      cn(dumI) = dum;
    end
    figure(1)
    subplot(1,sp,2)
    hold off
    ii = find(cn>dum);
    plot(jd0(ii),c(ii),'k')
    hold on
    plot(jd0(ii),cn(ii),'r')
    title('conductivity');xlim(xli);
  end %  cond. dummies loop 

  if ~isempty(find(~isnan(p)))
    inval= input('Should pressures in a certain time interval be set to dummies y/n \n','s');

    if inval == 'y' % cond dummies loop
      invalP= input('Enter interval, e.g. [2 4 6 9] denotes days 2-4 and days 6-9 \n');      

      while 1
        if isodd(length(invalP))
          inval= input('Repeat entry, interval must have even number of elements \n')      
        else
          break
        end
      end 

      for i = 1 : length(invalP)/2,  
        dumI = find(jd0>= invalP(i*2-1) & jd0 <= invalP(i*2));    
        pn(dumI) = dum;
      end
      figure(1)
      subplot(1,sp,3)
      hold off
      ii = find(pn>dum);
      plot(jd0(ii),p(ii),'k')
      hold on
      plot(jd0(ii),pn(ii),'r')
      title('pressure');xlim(xli);
    end %  cond. dummies loop 
  end

  %%%%%%%----- 4a)  apply constant offset to section of data -------%%%%%%%%%%
  
  invalTc = 0; invalCc = 0; invalPc = 0;
  invalofft = 0;  invaloffc = 0;  invaloffp = 0;
  
  inval = input('Does an offset need to be applied to part of the temperature series y/n \n','s');
    
    if strcmp(inval,'y')
        invalofft= input('Enter time interval where correction is required, e.g. [3 56] denotes days 3-56 ');
        invalTc= input('Insert offset: ');
  
      while 1
        if isodd(length(invalofft))
          inval= input('Repeat entry, interval must have even number of elements \n')      
        else
          break
        end
      end 
  
      for i = 1:length(invalofft)/2,
        dumI = find(jd0>= invalofft(i*2-1) & jd0 <= invalofft(i*2));    
        tn(dumI) = tn(dumI)+invalTc(i);
      end
      figure(1)
    subplot(1,sp,1)
    hold off
    ii = find(tn>dum);
    plot(jd0(ii),t(ii),'k')
    hold on
    plot(jd0(ii),tn(ii),'b')
    title('temperature');xlim(xli);
    end
    
  inval = input('Does an offset need to be applied to part of the conductivity series y/n \n','s');
   
    if strcmp(inval,'y')
      if strcmp(mooring,'eb1_9_201113')&&~isempty(findstr(mcfile,'5763'))
        cn(1387)=cn(1387)-0.26;
      end
        invaloffc= input('Enter time interval where correction is required, e.g. [3 56] denotes days 3-56 ');
        invalCc= input('Insert offset: ');
  
      while 1
        if isodd(length(invaloffc))
          inval= input('Repeat entry, interval must have even number of elements \n')      
        else
          break
        end
      end 
  
      for i = 1:length(invaloffc)/2,
        dumI = find(jd0>= invaloffc(i*2-1) & jd0 <= invaloffc(i*2));    
        cn(dumI) = cn(dumI)+invalCc(i);
      end
      
    figure(1)
    subplot(1,sp,2)
    hold off
    ii = find(cn>dum);
    plot(jd0(ii),c(ii),'k')
    hold on
    plot(jd0(ii),cn(ii),'r')
    title('conductivity');xlim(xli);
    
    
      figure(6)
      hold on
      if strcmp(acc_cpcor,'y')
       sn = sw_salt(cn(ii)/c3515,tn(ii)*t90_68,pref*ones(length(ii),1));
      else
       sn = sw_salt(cn(ii)/c3515,tn(ii)*t90_68,pn(ii));
      end
    plot(jd(ii)-jd(1),sn,'r')
    grid on
    end
  
  
    inval = input('Does an offset need to be applied to part of the pressure series y/n \n','s');
    
    if inval == 'y'
        invaloffp= input('Enter time interval where correction is required, e.g. [3 56] denotes days 3-56 ');
        invalPc= input('Insert offset: ');
  
      while 1
        if isodd(length(invaloffp))
          inval= input('Repeat entry, interval must have even number of elements \n')      
        else
          break
        end
      end 
  
      for i = 1:length(invaloffp)/2,
        dumI = find(jd0>= invaloffp(i*2-1) & jd0 <= invaloffp(i*2));    
        pn(dumI) = pn(dumI)+invalPc(i);
      end
      figure(1)
      subplot(1,sp,3)
      hold off
      ii = find(pn>dum);
      plot(jd0(ii),p(ii),'k')
      hold on
      plot(jd0(ii),pn(ii),'r')
      title('pressure');xlim(xli);
    end
 
%%%%%%%%%%%%%% 5) interactive despiking  %%%%%%%%%%%%%%%%%%
  c3515   = sw_c3515;

  val = find(cn>dum & pn>dum);  

  if isempty(val)
    val = find(cn>dum);
    sn = sw_salt(cn(val)/c3515,tn(val)*t90_68,pref*ones(length(val),1));
  else    
    sn = sw_salt(cn(val)/c3515,tn(val)*t90_68,pn(val));
  end 

  figure(111);hold off
  if ~isempty(find(pn>dum))
      plot(sn,sw_ptmp(sn,tn(val),pn(val),median(pn(val))),'.k')
  else       
      plot(sn,tn(val),'.r')
  end 
  grid on
  
for zz=1:size(external_ctd_file,1)
    if exist([external_ctd(zz,1:max(findstr(external_ctd(zz,:),'/'))),external_ctd_file(zz,(~isspace(external_ctd_file(zz,:))))]) == 2
      eval(['load ',external_ctd(zz,1:max(findstr(external_ctd(zz,:),'/'))),external_ctd_file(zz,:)])
  clear cnt dis
      for cnt = 1 :length(ctd_prof)
        dis(cnt) = dist2([ lt ctd_lat(cnt)],[ln ctd_lon(cnt)]);
      end
  clear near
  near = find(dis<distrangectd);
  distance=dis(near)/1000;
  diag=fopen('CTD_near.txt','a');      
  fprintf(diag,'%s',[external_ctd_file(zz,:)]);fprintf(diag,'\n');
      fprintf(diag,'%4.2f\t',near);
      fprintf(diag,'%.3f\n',distance);
  %    fclose(diag);
      for qq=1:length(near)
   if ~isempty(strfind(external_ctd(zz,:),'kn221'))
       ctd_t(ctd_t==0)=nan;
       ctd_s(ctd_s==0)=nan;   
       hold on
      xli = get(gca,'Xlim');
      yli = get(gca,'Ylim');

      if ~isempty(find(pn>dum))
        %%plot(ctd_s,theta(ctd_p,ctd_t,ctd_s,median(pn(val))),'g')
        plot(ctd_s,sw_ptmp(ctd_s,ctd_t,ctd_p,median(pn(val))),'g')          
      else     
        plot(ctd_s,ctd_t,'y')
      end 

      xlim(xli);ylim(yli);
      
   elseif ~isempty(strfind(external_ctd(zz,:),'pe399'))
   ctd_file = [external_ctd(zz,1:max(findstr(external_ctd(zz,:),'/'))),external_ctd_file(zz,1:5),'_',sprintf('%3.3d',ctd_prof(near(qq))),'.ctd'];
       [ctd_p,ctd_t,ctd_s] = rodbload(ctd_file,'p:t:s');
       ctd_t(ctd_t==0)=nan;
       ctd_s(ctd_s==0)=nan;       
      hold on
      xli = get(gca,'Xlim');
      yli = get(gca,'Ylim');
      xli = get(gca,'Xlim');
      yli = get(gca,'Ylim');
      if ~isempty(find(pn>dum))
        %%plot(ctd_s,theta(ctd_p,ctd_t,ctd_s,median(pn(val))),'g')
        plot(ctd_s,sw_ptmp(ctd_s,ctd_t,ctd_p,median(pn(val))),'b')
      else     
        plot(ctd_s,ctd_t,'y')
 
      end       
      
elseif ~isempty(strfind(external_ctd(zz,:),'pe400'))
       ctd_t(ctd_t==0)=nan;
       ctd_s(ctd_s==0)=nan;   
       hold on
      xli = get(gca,'Xlim');
      yli = get(gca,'Ylim');

      if ~isempty(find(pn>dum))
        %%plot(ctd_s,theta(ctd_p,ctd_t,ctd_s,median(pn(val))),'g')
        plot(ctd_s,sw_ptmp(ctd_s,ctd_t,ctd_p,median(pn(val))),'b')
              
      else     
        plot(ctd_s,ctd_t,'y')
 
      end 

      xlim(xli);ylim(yli);
else
   ctd_file = [external_ctd(zz,1:max(findstr(external_ctd(zz,:),'/'))),external_ctd_file(zz,1:4),'_',sprintf('%3.3d',ctd_prof(near(qq))),'.ctd'];
       [ctd_p,ctd_t,ctd_s] = rodbload(ctd_file,'p:t:s');
      hold on
      xli = get(gca,'Xlim');
      yli = get(gca,'Ylim');

      if ~isempty(find(pn>dum))
          
                if ~(isempty(strfind(external_ctd(zz,:),'d382')))   %'sj08')))                
                      plot(ctd_s,sw_ptmp(ctd_s,ctd_t,ctd_p,median(pn(val))),'r')
                elseif ~(isempty(strfind(external_ctd(zz,:),'jc103')))   
                      plot(ctd_s,sw_ptmp(ctd_s,ctd_t,ctd_p,median(pn(val))),'m')%color',[0.5 0.5 0.5])
                else
                     plot(ctd_s,sw_ptmp(ctd_s,ctd_t,ctd_p,median(pn(val))),'b')
                end
      else     
        plot(ctd_s,ctd_t,'y')
 
      end 
   end
      xlim(xli);ylim(yli);
   end
end
end

  ELIM = []; 
  if exist('diag')==1
  fclose(diag);
  open('CTD_near.txt')
  end
  while 1   
    cut = input('Want to exclude spikes from T/S plot interactively y/n ','s');
   
    if strcmp(cut,'y') 
      fprintf(1,'Define an area that covers spikes by clicking on graph \n') 
      [x,y] = ginput;
      x = [x ;x(1)];
      y = [y ;y(1)];
      elim = find(inpolygon(sn,tn(val),x,y) == 1);
     
      ELIM = [ELIM;elim];
      hold on
      if ~isempty(find(pn>dum))
         plot(sn(elim),sw_ptmp(sn(elim),tn(val(elim)),pn(val(elim)),median(pn(val))),'.w') 
      else       
         plot(sn(elim),tn(val(elim)),'.w')
      end 
      
     
    else
      figure(1)
      subplot(1,sp,2)
      
      cn(val(ELIM)) = dum;
      sn(ELIM)      = dum;
      ii = find(cn>dum);
      plot(jd0(ii),cn(ii),'g')
      fprintf(1,'%d data points eliminated \n\n\n',length(ELIM)) 
      fprintf(1,[' Name of the figure '  mcfig_out  '.png' '\n\n'])
      print([mcfig_out, '.png'],'-dpng');

      break
    end    
  end % end while
  
    % automatic de-spiking
    
    val   = find(cn> dum &pn>dum);
    if isempty(val)
      val = find(cn>dum);
      sn  = sw_salt(cn(val)/c3515,tn(val)*t90_68,pref*ones(length(val),1));
    else    
      sn =  sw_salt(cn(val)/c3515,tn(val)*t90_68,pn(val));
    end 
       ELIM2 = [];
    if ~isempty(val) 
      Tlim  = [min(tn(val)) max(tn(val))];
      dTlim = diff(Tlim); 
      Tstep = 15;
      Tgrid = linspace(Tlim(1),Tlim(2),Tstep);
   
   
    
      for i = 1 : Tstep -1
        ii     = find(tn(val)>=Tgrid(i) & tn(val)<=Tgrid(i+1));
        if ii < 3
          ssd(i) = ssd(i-1);  
        else
          ssd(i) = std(sn(ii));
        end  
        
        smd(i) = median(sn(ii));
        elim  =          find(sn(ii) > (smd(i)+6*ssd(i)) | sn(ii) < (smd(i) - 6*ssd(i)));
        ELIM2 = [ELIM2 ii(elim)'];
      end
      Tgrid = mean([Tgrid(1:end-1);Tgrid(2:end)]);
      figure(111) 
      hold on
      plot(smd+6*ssd,Tgrid,'m--','Linewidth',2)
      plot(smd-6*ssd,Tgrid,'m--','Linewidth',2)
    
      auto = input('Do you want to exclude all the values outside the 6 sigma area y/n ','s');
      if strcmp(auto,'y')
        disp(['Automatic despiking accepted - ',num2str(length(ELIM2)),' values discarded']);  
        cn(val(ELIM2)) = dum;
        plot(sn(ELIM2),tn(val(ELIM2)),'w.')
      elseif strcmp(auto,'n')
        disp('Automatic despiking rejected')    
        ELIM2 = [];
      else 
        disp('Answer not defined - Automatic despiking rejected')
        ELIM2 = [];
      end  
    end
    saveas(gcf,'111.fig')
 %% ---  replace NaN by dummy
 
 ii     = find(isnan(pn));
 pn(ii) = dum;
 ii     = find(isnan(cn));
 cn(ii) = dum;
 ii     = find(isnan(tn));
 tn(ii) = dum;
 
 
%%%%%%%%%%% ---- 6) save data ------%%%%%%%%%%% 

if skipT == 1 & skipC == 1 & skipP == 1    
    if isempty(find(~isnan(p)))
      dat = [yy mm dd hh tn cn];   
      rodbsave(mcfile_out,head,fort,moo,serial0,wd,id,sd,st,ed,et,lt,ln,cols,dat)
    else
      dat = [yy mm dd hh tn cn pn];
      rodbsave(mcfile_out,head,fortp,moo,serial0,wd,id,sd,st,ed,et,lt,ln,colsp,dat) % instr with pressure option
    end
    fprintf(1,[mcfile_out,' has been saved\n\n\n'])
else  
    fprintf(1,[mcfile_out,' has NOT been saved\n\n\n'])
end 
 
 
  
 %%% ---- Text Output -------------------------------------
 
 fidtxt  =  fopen([mcfile_out,'.txt'],'w');
 
 fprintf(fidtxt,'Microcat_apply_cal.m: \n');
 fprintf(fidtxt,'Date    : %s \n',date);
 fprintf(fidtxt,'Operator: %s \n',operator);
 fprintf(fidtxt,'Input file : %s \n',mcfile);
 fprintf(fidtxt,'Output file: %s \n',mcfile_out);
 fprintf(fidtxt,'Variable            | pre-cruise | post-cruise \n');
 fprintf(fidtxt,'Conductivity [mS/cm]:  %5.4f        %5.4f \n',preC(cinI),postC(cinI));
 fprintf(fidtxt,'Temperature  [degC]    :  %5.4f        %5.4f \n',preT(tinI),postT(tinI));
 fprintf(fidtxt,'Pressure     [dbar] :  %5.4f        %5.4f \n',preP(pinI),postP(pinI));
 if exist('applytrendc')
 fprintf(fidtxt,'Average conductivity applied? %s\n',applytrendc);
 end
 if exist('applytrendcg')
 fprintf(fidtxt,'General conductivity trend applied? %s\n',applytrendcg);
 end
 if exist('applytrendt')
 fprintf(fidtxt,'Average temperature applied? %s\n',applytrendt);
 end
 if exist('applytrendtg')
 fprintf(fidtxt,'General temperature trend applied? %s\n',applytrendtg);
 end
 if exist('applytrendp')
 fprintf(fidtxt,'Average pressure applied? %s\n',applytrendp);
 end
 if exist('applytrendpg')
 fprintf(fidtxt,'General pressure trend applied? %s\n',applytrendpg);
 end
 fprintf(fidtxt,'Pressure drift removal? %s\n',acc_drift);
 if strcmp(acc_drift,'y')
     fprintf(fidtxt,'Equation of pressure drift fit: yfit = a1*(1-exp(-a2*(x-x(1)))) + a3*(x-x(1)) + a4\n');
     fprintf(fidtxt,'Coefficients of pressure drift fit: a1 = %5.6f, a2 = %5.6f, a3 = %5.6f, a4 = %5.6f\n',coef(1),coef(2),coef(3),coef(4));
 end
 fprintf(fidtxt,'Conductivity pressure correction redone? %s\n',acc_cpcor);
 if strcmp(acc_cpcor,'y')
  fprintf(fidtxt,'Reference pressure for conductivity pressure correction: %s\n',num2str(pref));
  fprintf(fidtxt,'Time interval to which correction was applied: [%s]\n',num2str(invalcpcor));
 end
 fprintf(fidtxt,'Skipped conductivity intervals: [%s]\n',num2str(invalC));
 fprintf(fidtxt,'Skipped temperature intervals : [%s]\n',num2str(invalT));
 fprintf(fidtxt,'Skipped pressure intervals    : [%s]\n',num2str(invalP));
 if exist('invaloffc')
     for i=1:length(invaloffc)/2
 fprintf(fidtxt,'Offset of %5.4f applied to conductivity interval [%s]\n',invalCc(i),num2str(invaloffc(i*2-1:i*2)));
     end
 end
 if exist('invalofft')
     for i=1:length(invalofft)/2
 fprintf(fidtxt,'Offset of %5.4f applied to temperature interval [%s]\n',invalTc(i),num2str(invalofft(i*2-1:i*2)));
     end
 end
 if exist('invaloffp')
     for i=1:length(invaloffp)/2
 fprintf(fidtxt,'Offset of %5.4f applied to pressure interval [%s]\n',invalPc(i),num2str(invaloffp(i*2-1:i*2)));
     end
 end
 fprintf(fidtxt,'Number of additional C points skipped interactively: [%d]\n',length(ELIM));
 fprintf(fidtxt,'Number of additional C points skipped automatically: [%d]\n',length(ELIM2));

end

fclose all
