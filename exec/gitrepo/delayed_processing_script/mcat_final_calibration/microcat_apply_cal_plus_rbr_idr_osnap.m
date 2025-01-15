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
% edited by Tiago Dotto Jan 2025 - added the progdir in global, modified
%     the pathway of p_applycal.coef_dir, changed the save mode of .fig

close all; warning off

global basedir datadir execdir pathgit pathosnap progdir % 14/01/2025

% path of the mooring data define in the startup file under osnap/

% moor = 'rteb1_07_2022';
% moor = 'rtwb1_07_2022';
%moor = 'rtwb2_07_2022';
moor = 'ib5_03_2022';
%=========================================================================
% Apply calibration coefficients to series, removes bad data. If required, applies
% constant offsets, and conductivity pressure correction
p_applycal.operator  = 'TSD';
p_applycal.mooring  = moor;   
p_applycal.sensortyp = 'microcat';%'microcat';   % seaphox / arg / microcat / rbr / idr
p_applycal.delim = ',';
% input directories & files
p_applycal.mooring_dir         = [basedir '/osnap/data/moor/proc/'];
p_applycal.mooring_outdir      = [basedir '/osnap/data/moor/proc/'];
p_applycal.coef_dir            = [progdir '/m_moorproc_toolbox/metadata/cal_coef/osnap/']; % 14/01/2025
p_applycal.external_ctd_dir    = [basedir '/cruise_data/'];
p_applycal.ctd_ref_cruises     = {''};%{'pe400'}; %{'kn221-02';'pe399'}; % references cruises for the QC
p_applycal.distrangectd        = 100e3; % distance of the reference ctd from the mooring
p_applycal.strformat.mctemptxt = repmat('%s',1,38);
p_applycal.strformat.mctempnum = ['%f%f%s%s%f%s%s%f%f%s%s' repmat('%f%s%s',1,9)];
p_applycal.strformat.mcsaltxt  = repmat('%s',1,38);
p_applycal.strformat.mcsalnum  = ['%f%f%s%s%f%s%s%f%f%s%s' repmat('%f%s%s',1,9)];
p_applycal.strformat.mcprestxt = repmat('%s',1,81);
p_applycal.strformat.mcpresnum = ['%f' repmat('%s%f%s%s',1,20)];
loclegend = 'north';
opts.Interpreter = 'tex';
opts.Default = 'No';
opts.WindowStyle='normal';
% ---------------------------------------------------------------------------
operator  = p_applycal.operator; 
mooring   = p_applycal.mooring ;
sensortyp = p_applycal.sensortyp; % arg / microcat / rbr / idr
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
fclose(fid_cruise);

moorI                = find((strcmp(cruise{1},mooring))>0);
cruise               = [cruise{2} cruise{3}];
cruiseI              = cruise(moorI,:);

% ------ extract deployment can recovery cruise names from database--------
[num,txt,raw]=xlsread([coef_dir,'cruise_id.xls']);
deplcr = txt(find(num==cruiseI(1))+1);
reccr  = txt(find(num==cruiseI(2))+1);

% ------ output  data ---- 
dum      = -9999; 
pref     = 0;
t90_68   = 1.00024;  % convert its90 to its68 for cond. to sal.
c3515   = gsw_C3515 ;% Conductivity at (35,15,0)% GSW update

% ------ conversion ------------
cols      = 'YY:MM:DD:HH:T:C'; % column info for rodb header
colsp     = 'YY:MM:DD:HH:T:C:P'; % column info for rodb header (mc with pressure sensor)
colsphox  = 'YY:MM:DD:HH:T:S:P'; % column info for rodb header (mc with pressure sensor)
colsarg   = 'YY:MM:DD:HH:T:TCAT:P:PCAT:C:U:V:W';
fort      = '%4.4d  %2.2d  %2.2d  %7.5f   %6.4f  %6.4f'; %data output format
fortp     = '%4.4d  %2.2d  %2.2d  %7.5f   %6.4f  %6.4f  %5.1f'; %data output format(mc with pressure sensor)  

% ------ sensor specific settings ------

if strcmp('microcat',sensortyp)
    typ_id = [333 337];
elseif strcmp('arg',sensortyp)
    typ_id = 366;  
elseif strcmp('rbr',sensortyp)
    typ_id = 330;
elseif strcmp('idr',sensortyp)
    typ_id = 339;
elseif strcmp('seaphox',sensortyp)
    typ_id = 375;
end
ext = ['.',sensortyp];

% ---- microcat raw data ----
head    = ['Mooring:SerialNumber:WaterDepth:InstrDepth:Start_Date:Start_Time:End_Date:End_Time:Latitude:Longitude:Columns'];

% ---- load calib offsets (pre and post cruise) ---
if strcmp('microcat',sensortyp) | strcmp('idr',sensortyp) | strcmp('rbr',sensortyp)| strcmp('seaphox',sensortyp)

  % Temperature
  fidt  = fopen([coef_dir,'microcat_temp.csv'],'r');
  ttext = textscan(fidt,strformat.mctemptxt,'delimiter',delim);
  fseek(fidt,0,'bof');
  tnum  = textscan(fidt,strformat.mctempnum,'delimiter',delim,'HeaderLines',5);
  % Conductivity
  fidc  = fopen([coef_dir,'microcat_cond.csv'],'r');
  ctext = textscan(fidc,strformat.mcsaltxt ,'delimiter',delim);
  fseek(fidc,0,'bof');
  cnum  = textscan(fidc,strformat.mcsalnum,'delimiter',delim,'HeaderLines',5);
  % Pressure
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
    % ------ Conductivity ------
    for i = 1 :clen
        cvalI(i)=strcmp(ctext{i}{2},num2str(cruiseI(1)));
    end
    ccol = find(cvalI == 1);
    preC  =  cnum{ccol}; %pre cruise offsets
    preCcomment =  cnum{ccol+2};
      
    % ------ Temperature ------
    for i = 1 :tlen
        tvalI(i)=strcmp(ttext{i}{2},num2str(cruiseI(1)));
    end
    tcol = find(tvalI == 1);
    preT  =  tnum{tcol}; %pre cruise offsets
    preTcomment =  tnum{tcol+2};

  % ------ Pressure ------
  for i = 1 :plen
    pvalI(i)=strcmp(ptext{i}{2},num2str(cruiseI(1)));
  end
  pcol = find(pvalI == 1);
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

    % ------ Conductivity ------
    for i = 1 :clen
        cvalI(i)=strcmp(ctext{i}{2},num2str(cruiseI(2)));
    end
    ccol = find(cvalI == 1);
    postC  = cnum{ccol}; %post cruise offsets
    postCcomment =  cnum{ccol+2};  
    
    % ------ Temperature ------
    for i = 1 :tlen
        tvalI(i)=strcmp(ttext{i}{2},num2str(cruiseI(2)));
    end
    tcol = find(tvalI == 1);
    postT  = tnum{tcol}; %post cruise offsets
    postTcomment =  tnum{tcol+2};  


    for i = 1 :plen
        pvalI(i)=strcmp(ptext{i}{2},num2str(cruiseI(2)));
    end
    pcol = find(pvalI == 1);

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

%% 1 apply dedrift 
[typ,dep,serial] = rodbload([mooring_dir,mooring,'/',mooring,'info.dat'],...
                        ['instrument:z:serialnumber']);
mcI     = find(typ >= typ_id(1) & typ <= typ_id(end));
dep     = dep(mcI);
serial  = serial(mcI);
typ     = typ(mcI);

for mc = 1:length(serial)
  close all
  nextmcat='';
  mcfile_out = [mooring_outdir,mooring,'/',sensortyp,'/',mooring,'_',sprintf('%3.3d',mc),ext]; 
  mcfig_out = [mooring_outdir,mooring,'/',sensortyp,'/',mooring,'_',sprintf('%3.3d',mc)]; 
  if exist(mcfile_out,'file')  
    answer = questdlg(sprintf('%s',['\fontsize{16}File ' [sprintf('%3.3d',mc) ext] ' for instrument ' num2str(serial(mc)) ' exists already: do you want to overwrite it? ']), ...
        'Existing calibrated data file', ...
        'Yes','Yes, but make a dated copy','No, skip to next sensor',opts);
    switch answer
       case 'Yes'
            writefile=1;
       case 'Yes, but make a dated copy'
            [succ1,mess1] =copyfile(mcfile_out,[mcfile_out,date]);
            [succ2,mess2] =copyfile([mcfile_out,'.txt'],[mcfile_out,'.txt',date]);
            writefile=0;
        case 'No, skip to next sensor'
            nextmcat='skip';
    end
  end  

  if strcmpi(nextmcat,'skip')
      continue
  end
  
  % default is no calibration data 
  skipC  = 1;skipT  = 1;skipP  = 1;
  
  mcfile = [mooring_dir,mooring,'/',sensortyp,'/',mooring,'_',sprintf('%4.4d',serial(mc)),'.use'];

  if exist(mcfile)==0
      warndlg(sprintf('%s',['File: ' mcfile ' does not exist!']),...
          'WARNING',...
          opts)
  end
  
  % read .USE file
  if strcmp(sensortyp,'microcat') | strcmp(sensortyp,'rbr') | strcmp(sensortyp,'idr') 
    [yy,mm,dd,hh,t,c,p]= rodbload(mcfile,colsp);
  elseif strcmp(sensortyp,'seaphox')
    [yy,mm,dd,hh,t,s,p]= rodbload(mcfile,colsphox); 
    c = gsw_C_from_SP(s,t,p);
  elseif strcmp(sensortyp,'arg')
    [yy,mm,dd,hh,t,tcat,p,pcat,c,u,v,w]= rodbload(mcfile,colsarg);
  end      
      
  [moo,serial0,wd,id,sd,st,ed,et,lt,ln,cl] = rodbload(mcfile,head);

  jd       = julian(yy,mm,dd,hh);
  jd0      = jd - jd(1);  
  xli      = [0 max(jd0)];
    
  cinI     = find(cin == serial0);
  tinI     = find(tin == serial0); 
  p_exist  = 'n';

  if ~isempty(find(~isnan(p)))
    p_exist = 'y';
    pinI     = find(pin == serial0); 
  end

  %% 1.1 CONDUCTIVITY
  warn = 0;
  if isempty(cinI)           
     corrC    =  0;
     answer = questdlg(sprintf('%s',['\fontsize{16}No conductivity coeffecient found for serial number ',num2str(serial0)]), ...
	    [num2str(serial(mc)) ' conductivity coefficient'], ...
	    'Skip','Save with cal applied','Or, break here',opts);
        % Handle response
        switch answer
           case 'Skip'
                postC(cinI) = NaN;
                disp('Skipping')
                skipC = 0;
           case 'Save with no cal applied'
                preC(cinI) = NaN;
                disp('Saving without calibration!!')
                skipC = 1;
           case 'Or, break here'
                disp('Breaking')
                return
        end
  else
    
    % transfer coefs to temprary variables
    cal1        =  preC(cinI);
    cal1comment = preCcomment(cinI);
    cal2        =  postC(cinI);       
    cal2comment = postCcomment(cinI);

    answer = questdlg(sprintf('%s\n%s\n%s\n%s',['\fontsize{16}Conductivity pre-cruise calibration: ' num2str(cal1)], ...
                ['\fontsize{16}Comment: ' cal1comment{1}],...
                ['\fontsize{16}Conductivity post-cruise calibration: ' num2str(cal2)],...
                ['\fontsize{16}Comment: ' cal2comment{1}]), ...
	[num2str(serial(mc)) ' conductivity coefficient selection'], ...
	'Use both','Just pre-cruise','Or, just post-cruise',opts);
    % Handle response
    switch answer
        case 'Just pre-cruise'
            postC(cinI) = NaN;
            disp('Only using the pre cruise calibration coefficient')
            keep = 1;
       case 'Or, just post-cruise'
            preC(cinI) = NaN;
            disp('Only using the post cruise calibration coefficient')
            keep = 2;
       case 'Use both'
            disp('Using both calibration coefficients')
            keep = 0;
    end

    % update calibration offsets for conductivity
    cal1     =  preC(cinI);
    cal2     =  postC(cinI);       
      
    % IF one of the offsets is a nan, we can use a general trend
    if isnan(cal1)
        answer = questdlg('\fontsize{16}Should general trend generated from all sensors be applied?', ...
	    [num2str(serial(mc)) ' conductivity pre-cruise'], ...
	    'Yes','No',opts);
        % Handle response
        switch answer
            case 'Yes'
                cal1   =  cal2-Ctrend;
            case 'No'
                cal1   =  cal2-Ctrend*0;
        end
      warn = 1;  
    end

    if isnan(cal2)
        answer = questdlg('\fontsize{16}Should general trend generated from all sensors be applied?', ...
	    [num2str(serial(mc)) ' conductivity post-cruise'], ...
	    'Yes','No',opts);
        % Handle response
        switch answer
            case 'Yes'
                cal2   =  cal1+Ctrend;
            case 'No'
                cal2   =  cal1+Ctrend*0;
        end
      warn = 2;  
    end

    if keep == 0
        answer = questdlg('\fontsize{16}Should average of conductivity coeffecients be applied?', ...
	    [num2str(serial(mc)) ' conductivity'], ...
	    'Yes','No',opts);
        % Handle response
        switch answer
            case 'Yes'
                cal1   =  (cal1+cal2)/2;
                cal2   =  cal1;
            case 'No'
                cal1   =  cal1;
                cal2   =  cal2;
        end
      warn = 3;  
    end  

    % calculate the correction to be applied to the data
    corrC    =  cal1 + (cal2 - cal1)/max(jd0)*jd0; 
    
    if isempty(find(~isnan(corrC)))
        corrC=0;
        answer = questdlg(sprintf('%s',['\fontsize{16}Conductivity calibration for serial number ',num2str(serial0) ' is empty!']), ...
	        'Conductivity calibration', ...
	        'Skip','Save without cal applied','Or, break here',opts);
            % Handle response
            switch answer
               case 'Skip'
                    disp('Skipping')
                    skipC = 0;
               case 'Save with cal applied'
                    disp('Saving without calibration!!')
                    skipC = 1;
               case 'Or, break here'
                    disp('Breaking')
                    return
            end  
    end
  end
  
%% 1.2 TEMPERATURE    
  if isempty(tinI)
     corrT    =  0;
      answer = questdlg(sprintf('%s',['\fontsize{16}No temperature coeffecient found for serial number ',num2str(serial0)]), ...
	    'Conductivity coefficient entry', ...
	    'Skip','Save with cal applied','Or, break here',opts);
        % Handle response
        switch answer
           case 'Skip'
                disp('Skipping')
                skipC = 0;
           case 'Save with cal applied'
                disp('Saving without calibration!!')
                skipC = 1;
           case 'Or, break here'
                disp('Breaking')
                return
        end
  else
    
    cal1comment = preTcomment(tinI); 
    cal2comment = postTcomment(tinI);
    cal1        =  preT(tinI);
    cal2        =  postT(tinI); 
        
    answer = questdlg(sprintf('%s\n%s\n%s\n%s',['\fontsize{16}Temperature pre-cruise calibration: ' num2str(cal1)], ...
                ['\fontsize{16}Comment: ' cal1comment{1}],...
                ['\fontsize{16}Temperature post-cruise calibration: ' num2str(cal2)],...
                ['\fontsize{16}Comment: ' cal2comment{1}]), ...
	'Temperature coefficient selection', ...
	'Use both','Just pre-cruise','Or, just post-cruise',opts);
    % Handle response
    switch answer
        case 'Just pre-cruise'
            postT(tinI) = NaN;
            disp('Only using the pre cruise calibration coefficient')
            keep = 1;
       case 'Or, just post-cruise'
            preT(tinI) = NaN;
            disp('Only using the post cruise calibration coefficient')
            keep = 2;
       case 'Use both'
            disp('Using both calibration coefficients')
            keep = 0;
    end
     
    cal1     =  preT(tinI);
    cal2     =  postT(tinI); 

    % IF one of the offsets is a nan, we can use a general trend
    if isnan(cal1)
        answer = questdlg('\fontsize{16}Should general trend generated from all sensors be applied?', ...
	    'Temperature pre-cruise', ...
	    'Yes','No',opts);
        % Handle response
        switch answer
            case 'Yes'
                cal1   =  cal2-Ttrend;
            case 'No'
                cal1   =  cal2-Ttrend*0;
        end
      warn = 1;  
    end

    if isnan(cal2)
        answer = questdlg('\fontsize{16}Should general trend generated from all sensors be applied?', ...
	    'Temperature post-cruise', ...
	    'Yes','No',opts);
        % Handle response
        switch answer
            case 'Yes'
                cal2   =  cal1+Ttrend;
            case 'No'
                cal2   =  cal1+Ttrend*0;
        end
      warn = 2;  
    end

    if keep == 0
        answer = questdlg('\fontsize{16}Should average of Temperature coeffecients be applied?', ...
            'Temperature', ...
            'Yes','No',opts);
        switch answer
            case 'Yes'
                cal1   =  (cal1+cal2)/2;
                cal2   =  cal1;
            case 'No'
                cal1   =  cal1;
                cal2   =  cal2;
        end
        warn = 3;  
    end

    corrT    =  cal1 + (cal2 - cal1)/max(jd0)*jd0; 

    if isempty(find(~isnan(corrT)))
        corrT=0;
        answer = questdlg(sprintf('%s',['\fontsize{16}No Temperature calibration data found for serial number ',num2str(serial0)]), ...
	        'Temperature calibration', ...
	        'Skip','Save without cal applied','Or, break here',opts);
            % Handle response
            switch answer
               case 'Skip'
                    disp('Skipping')
                    skipT = 0;
               case 'Save with cal applied'
                    disp('Saving without calibration!!')
                    skipT = 1;
               case 'Or, break here'
                    disp('Breaking')
                    return
            end  
    end
  end
%% 1.3 PRESSURE
if strcmp(p_exist,'y')
    if isempty(pinI)
        corrP   =  0;
        answer = questdlg(sprintf('%s',['\fontsize{16}No Pressure coeffecient found for serial number ',num2str(serial0)]), ...
	    [num2str(serial(mc))  'pressure coefficient entry'], ...
	    'Skip','Save with cal applied','Or, break here',opts);
        % Handle response
        switch answer
           case 'Skip'
                disp('Skipping')
                skipP = 0;
           case 'Save with cal applied'
                disp('Saving without calibration!!')
                skipP = 1;
           case 'Or, break here'
                disp('Breaking')
                return
        end
    else   
      cal1comment = prePcomment(pinI); 
      cal2comment = postPcomment(pinI);
      cal1        =  preP(pinI);
      cal2        =  postP(pinI); 
      
        answer = questdlg(sprintf('%s\n%s\n%s\n%s',['\fontsize{16}Pressure pre-cruise calibration: ' num2str(cal1)], ...
                    ['\fontsize{16}Comment: ' cal1comment{1}],...
                    ['\fontsize{16}Pressure post-cruise calibration: ' num2str(cal2)],...
                    ['\fontsize{16}Comment: ' cal2comment{1}]), ...
	    [num2str(serial(mc)) ' pressure coefficient selection'], ...
	    'Use both','Just pre-cruise','Or, just post-cruise',opts);
        % Handle response
        switch answer
            case 'Just pre-cruise'
                postP(pinI) = NaN;
                disp('Only using the pre cruise calibration coefficient')
                keep = 1;
           case 'Or, just post-cruise'
                preP(pinI) = NaN;
                disp('Only using the post cruise calibration coefficient')
                keep = 2;
           case 'Use both'
                disp('Using both calibration coefficients')
                keep = 0;
        end
     
      disp(' ')
      cal1     =  preP(pinI);
      cal2     =  postP(pinI); 
      
    % IF one of the offsets is a nan, we can use a general trend
    if isnan(cal1)
        answer = questdlg('\fontsize{16}Should general trend generated from all sensors be applied?', ...
	    [num2str(serial(mc)) ' pressure pre-cruise'], ...
	    'Yes','No',opts);
        % Handle response
        switch answer
            case 'Yes'
                cal1   =  cal2-Ptrend;
            case 'No'
                cal1   =  cal2-Ptrend*0;
        end
      warn = 1;  
    end

    if isnan(cal2)
        answer = questdlg('\fontsize{16}Should general trend generated from all sensors be applied?', ...
	    [num2str(serial(mc)) ' pressure post-cruise'], ...
	    'Yes','No',opts);
        % Handle response
        switch answer
            case 'Yes'
                cal2   =  cal1+Ptrend;
            case 'No'
                cal2   =  cal1+Ptrend*0;
        end
      warn = 2;  
    end

    if keep == 0
        answer = questdlg('\fontsize{16}Should average of Pressure coeffecients be applied?', ...
            [num2str(serial(mc)) ' pressure'], ...
            'Yes','No',opts);
        switch answer
            case 'Yes'
                cal1   =  (cal1+cal2)/2;
                cal2   =  cal1;
            case 'No'
                cal1   =  cal1;
                cal2   =  cal2;
        end
        warn = 3;  
    end
    
    corrP    =  cal1 + (cal2 - cal1)/max(jd0)*jd0; 
        if isempty(find(~isnan(corrP)))
        corrP=0;
        answer = questdlg(sprintf('%s',['\fontsize{16}No Pressure calibration data found for serial number ',num2str(serial0)]), ...
            [num2str(serial(mc)) ' pressure calibration'], ...
            'Skip','Save without cal applied','Or, break here',opts);
            % Handle response
            switch answer
               case 'Skip'
                    disp('Skipping')
                    skipP = 0;
               case 'Save with cal applied'
                    disp('Saving without calibration!!')
                    skipP = 1;
               case 'Or, break here'
                    disp('Breaking')
                    return
            end    
        end    
    end 
   
    pn       = p - corrP;
    
end

%% 2 apply offsets to T and C and plot
    
    tn       = t - corrT;      
    cn       = c - corrC;
    
    %find dummy data and apply to corrected data
    spikeT   = find(t<-900);tn(spikeT) = dum;
    spikeC   = find(c<-900);cn(spikeC) = dum; 
    spikeP   = find(p<-0); pn(spikeP) = dum;
    
    %index of depsiked data
    cidx=cn>dum;
    tidx=tn>dum;
    pidx=pn>dum; 
    
    figure(1);sgt = sgtitle([mooring '  ' p_applycal.sensortyp ' sn',num2str(serial0) ' at ' num2str(nanmean(p),'%4.0f') ' m'],...
      'Color','black', 'Interpreter', 'none');
    sgt.FontSize = 20;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    ax1=subplot(4,3,1:2);
    ax2=subplot(4,3,4:5);
    ax3=subplot(4,3,7:8);
    ax4=subplot(4,3,3);
    ax5=subplot(4,3,[6 9]);
    ax6=subplot(4,3,10);
    ax7=subplot(4,3,11);
    ax8=subplot(4,3,12);

    axes(ax1);
    plot(jd0(tidx),t(tidx),'k')
    hold on
    plot(jd0(tidx),tn(tidx),'b')
    grid on;  xlim(xli);
    ylabel('Temperature [deg C]')
    legend('Temperature raw','Temperature plus offset','location','best');
    title('1. Temperature');
    ax1.TitleHorizontalAlignment = 'left';
    
    axes(ax2);
    plot(jd0(cidx),c(cidx),'k')
    hold on
    plot(jd0(cidx),cn(cidx),'r')
    grid on;xlim(xli);
    ylabel('Conductivity [mS/cm]')
    legend('Conductivity raw','Conductivity plus offset','location','best');
    title('2. Conductivity');
    ax2.TitleHorizontalAlignment = 'left';
    
    if strcmp(p_exist,'y')
        axes(ax3);
        plot(jd0(pidx),p(pidx),'k')
        hold on
        plot(jd0(pidx),pn(pidx),'green')
        grid on; xlim(xli);
        ylabel('Pressure [db]')
        legend('Pressure raw','Pressure plus offset','location','best');
        title('3. Pressure');
        ylabel('Number of days since deployment')
        ax3.TitleHorizontalAlignment = 'left';
    end
    
    axes(ax6);
    plot(jd0([1 end]),corrT([1 end]),'bo-')
    grid on;
    ylabel('Temperature [deg C]')
    title('6. Temperature correction');
    ax6.TitleHorizontalAlignment = 'left';
    axes(ax7);
    plot(jd0([1 end]),corrC([1 end]),'ro-')
    ylabel('Conductivity [mS/cm]');
    grid on
    title('7. Conductivity correction');
    ax7.TitleHorizontalAlignment = 'left';
    if strcmp(p_exist,'y')
        axes(ax8);
        plot(jd0([1 end]),corrP([1 end]),'go-')
        ylabel('Pressure [db]'); grid on;
        title('8. Pressure correction');
        ax8.TitleHorizontalAlignment = 'left';
    end 
    
    % index of cond and temp dummy values
    % need to meake ctpidx
    ctidx = cidx & tidx & pidx;

    % Salinity from cndr, T, P
    % GSW update
    %SP from conductivity ratio(PSS‐78),  sw_salt(R,t,p), gsw_SP_from_R(R,t,p) 
    sn(ctidx) = gsw_SP_from_R(cn(ctidx)/c3515,tn(ctidx)*t90_68,pn(ctidx));
    s = gsw_SP_from_R(c/c3515,t*t90_68,p);
    
    axes(ax4);
    plot(jd(ctidx)-jd(1),s(ctidx),'k');
    hold on
    plot(jd(ctidx)-jd(1),sn(ctidx),'r');
    legend('Pre-cal','post-cal','Location','Best')
    ylabel('Salinity')
    title('4. Post-calibration Salinity');
    grid on
    ax4.TitleHorizontalAlignment = 'left';
    
    axes(ax5);
    plot(s(ctidx),t(ctidx),'k.');
    hold on
    plot(sn(ctidx),t(ctidx),'r.');
    grid on
    title('5. Salinity offset applied')
    legend('pre-cal','post-cal','location','best');
    xlabel('Salinity')
    ylabel('Temperature')
    ax5.TitleHorizontalAlignment = 'left';
    hold off
%%  3 conductivity pressure correction (optional)
    % only required for instruments withour pressure sensor
      if strcmp(p_exist,'n')
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
          figure(2)
            sn = gsw_SP_from_R(cn(ctidx)/c3515,tn(ctidx)*t90_68,pref*ones(length(ctidx),1));
            plot(jd(ctidx)-jd(1),sn,'r')
            grid on
            legend('raw','P corr.')
        elseif strcmp(acc_cpcor,'n')
          disp('No pressure conductivity correction applied')
        elseif ~strcmp(acc_cpcor,'y') & ~strcmp(acc_cpcor,'n')
          disp('Answer not defined - no pressure conductivity correction applied')  
          acc_cpcor = 'n';      
        end
      end
    %% 4 pressure drift removal 
    % need to sort out the inclsion of dummy values in here
    % only deals with exponential drift and is potentially skewed by
    % knockdown. perhaps we need something that deals with linear drift
    group = "Updates";
    pref = "Conversion";
    boxtitle = "Pressure drift removal";
    quest = "Does the pressure record require drift removal?";
    pbtns = ["Yes","No"];
    [pval,~] = uigetpref(group,pref,boxtitle,quest,pbtns,opts);
    
    switch pval
        case 'yes'
            acc_drift   =  1;
        case 'no'
            acc_drift   =  0;
    end
    
    if acc_drift==1
        [coef,fit]= exp_lin_fit2(jd,pn,[1 1 1 1]);
        if warn == 1
            pp = pn - fit + fit(end);
        elseif warn == 2
            pp = pn - fit + fit(1);
        else 
            pp = pn - fit + mean(fit);
        end
        
        so = gsw_SP_from_R(cn(ctidx)/c3515,tn(ctidx)*t90_68,pn(ctidx));
        sn = gsw_SP_from_R(cn(ctidx)/c3515,tn(ctidx)*t90_68,pp(ctidx));
        
        deltaS = std(sn)-std(so);
        fprintf(1,'\n\n Salinity standard deviation should decrease with dedrifted pressure \n')
    
        if deltaS > 0
            fprintf(1,' Salinity standard deviation has INCREASED instead  by %4.4f \n\n',deltaS)
        elseif deltaS < 0
            fprintf(1,' WARNING: Salinity standard deviation has DECREASED by %4.4f \n\n',-deltaS)
        end
        
        figure(1);axes(ax3); hold off;
        plot(jd0,pn)
        hold on
        plot(jd0,pp,'g')
        plot(jd0,fit,'r','Linewidth',2);
        grid on
        xlim(xli);
        ylabel('Pressure [db]')
        title('3. Pressure');
        legend('depsiked + offset','exponential linear fit applied',...
          'data minus fit','location','best');
    
        figure(1);axes(ax5);hold on
        plot(so,tn(ctidx),'b.')
        plot(sn,tn(ctidx),'r.')
        grid on; xlabel('S'); ylabel('T')
        legend('P raw','P dedrifted')
        title('5. Pressure drift applid to salinity and pressure');
        hold off
        group = "Updates";
        pref = "Conversion";
        boxtitle = "Pressure drift removal";
        quest = ["Do you accept the the drift removal??"];
        pbtns = ["Yes","No"];
        [pval,tf] = uigetpref(group,pref,boxtitle,quest,pbtns,opts);
        
        switch pval
            case 'yes'
                acc   =  1;
            case 'no'
                acc   =  0;
    end
    
      if acc == 1
        disp('drift removal accepted')
        
        figure(1); axes(ax3); hold off
        plot(jd0(ctidx),pn(ctidx),'k')
        hold on
        plot(jd0(ctidx),pp(ctidx),'g')    
        legend('no fit applierd','P corrected with linear exp fit')
        grid on        
        title('3. Pressure (drift correction applied)');
        ylabel('Number of days since deployment')
        ax3.TitleHorizontalAlignment = 'left';
        pn =pp;
        axes(ax4)
        plot(jd(ctidx)-jd(1),sn,'r')
        
        % apply correction
        pn(spikeP)=dum;
      else
        disp('drift removal discarded')
        acc_drift = 'n';
        axes(ax3);
        plot(jd0(pidx),p(pidx),'k')
        hold on
        plot(jd0(pidx),pn(pidx),'green')
        grid on; xlim(xli);
        ylabel('Pressure [db]')
        legend('Pressure raw','Pressure plus offset','location','best');
        title('3. Pressure (no drift removal)');
        ylabel('Number of days since deployment')
        ax3.TitleHorizontalAlignment = 'left';        
      end    
    
    elseif acc_drift==0
        disp('No drift correction applied')
    end    

    %% 5  set suspicious data to dummies
    % set default 
    invalT = 0; invalC = 0; invalP = 0;

    % TEMPERATURE----------------------------------------------------------
    group = "Updates";
    pref = "Conversion";
    boxtitle = "Temperature removal";
    quest = ["Should temperatures in a certain time interval be set to dummies?"];
    pbtns = ["Yes","No"];
    [pval,tf] = uigetpref(group,pref,boxtitle,quest,pbtns,opts);

    if strcmpi(pval,'yes')
          answer = inputdlg('Enter interval as start and end days e.g. 200 230 400 480',...
                     'Sample', [1 50],{' '},opts);
          invalT=str2num(answer{1});
          while 1
              if isodd(length(invalT))
                  answer = inputdlg('Entry must have even numbr of elements e.g. 200 230 400 480',...
                             'Sample', [1 50],{' '},opts);
                  invalT=str2num(answer{1});
              else
                break
              end
           end 
    
        for i = 1 : length(invalT)/2,  
          dumI = find(jd0>= invalT(i*2-1) & jd0 <= invalT(i*2));    
          tn(dumI) = dum;
        end

        figure(1)
        axes(ax1); hold off
        ii = find(t>dum);
        plot(jd0(ii),t(ii),'k')
        hold on
        ii = find(tn>dum);
        plot(jd0(ii),tn(ii),'b')
        grid on
        xlim(xli);  ylabel('Temperature [deg C]')
        legend('Pre-cal','post-cal','Location','Best')
    end 

    % CONDUCTIVITY---------------------------------------------------------
    group = "Updates";
    pref = "Conversion";
    boxtitle = "Conductivity removal";
    quest = ["Should conductivity in a certain time interval be set to dummies?"];
    pbtns = ["Yes","No"];
    [pval,tf] = uigetpref(group,pref,boxtitle,quest,pbtns,opts);

    if strcmpi(pval,'yes') 
          answer = inputdlg('Enter interval as start and end days e.g. 200 230 400 480',...
                     'Sample', [1 50],{' '},opts);
          invalC=str2num(answer{1});
          while 1
              if isodd(length(invalC))
                  answer = inputdlg('Entry must have even numbr of elements e.g. 200 230 400 480',...
                             'Sample', [1 50],{' '},opts);
                  invalC=str2num(answer{1});
              else
                break
              end
           end 
    
        for i = 1 : length(invalC)/2,  
          dumI = find(jd0>= invalC(i*2-1) & jd0 <= invalC(i*2));    
          cn(dumI) = dum;
        end

        figure(1)
        axes(ax2); hold off
        ii = find(c>dum);
        plot(jd0(ii),t(ii),'k')
        hold on
        ii = find(cn>dum);
        plot(jd0(ii),cn(ii),'r')
        grid on
        xlim(xli);  ylabel('Conductivity [mS/cm]')
        legend('Pre-cal','post-cal','Location','Best')
    end 

    % PRESSURE------------------------------------------------------------
    group = "Updates";
    pref = "Conversion";
    boxtitle = "Pressure removal";
    quest = ["Should pressures in a certain time interval be set to dummies?"];
    pbtns = ["Yes","No"];
    [pval,tf] = uigetpref(group,pref,boxtitle,quest,pbtns,opts);

    if strcmpi(pval,'yes')
          answer = inputdlg('Enter interval as start and end days e.g. 200 230 400 480',...
                     'Sample', [1 50],{' '},opts);
          invalP=str2num(answer{1});
          while 1
              if isodd(length(invalP))
                  answer = inputdlg('Entry must have even numbr of elements e.g. 200 230 400 480',...
                             'Sample', [1 50],{' '},opts);
                  invalP=str2num(answer{1});
              else
                break
              end
           end 
    
        for i = 1 : length(invalP)/2,  
          dumI = find(jd0>= invalP(i*2-1) & jd0 <= invalP(i*2));    
          pn(dumI) = dum;
        end

        figure(1)
        axes(ax3); hold off
        ii = find(p>dum);
        plot(jd0(ii),t(ii),'k')
        hold on
        ii = find(pn>dum);
        plot(jd0(ii),tn(ii),'g')
        grid on
        xlim(xli);    ylabel('Pressure [db]')
        legend('Pre-cal','post-cal','Location','Best')
    end 

  %% 6 apply constant offset to section of data
    % set default values
    invalTc = 0; invalCc = 0; invalPc = 0;
    invalofft = 0;  invaloffc = 0;  invaloffp = 0;
    
    % 6.1 TEMPERATURE
    group = "Updates";
    pref = "Conversion";
    boxtitle = "Offset";
    quest = ["Does an offset need to be applied to part of the temperature series?"];
    pbtns = ["Yes","No"];
    [pval,~] = uigetpref(group,pref,boxtitle,quest,pbtns,opts);

    if strcmpi(pval,'yes')
          answer = inputdlg({'Enter interval where correction is required e.g. 200 230 400 480',...
                                'Enter offset value'},...
                     'Enter offset', [1 50; 1 50],{' ',' '},opts);
          invalofft=str2num(answer{1});
          invalTc=str2num(answer{2});
      while 1
        if isodd(length(invalofft))
          answer = inputdlg({'not even! - Enter interval where correction is required e.g. 200 230 400 480',...
                                'Enter offset value'},...
                     'Enter offset', [1 50; 1 50],{' ',' '},opts);
          invalofft=str2num(answer{1});
          invalTc=str2num(answer{2});
        else
          break
        end
      end 
  
    for i = 1:length(invalofft)/2
        dumI = find(jd0>= invalofft(i*2-1) & jd0 <= invalofft(i*2));    
        tn(dumI) = tn(dumI)+invalTc(i);
    end

        figure(1)
        axes(ax1)
        hold off
        ii = find(tn>dum);
        plot(jd0(ii),t(ii),'k')
        hold on
        plot(jd0(ii),tn(ii),'b')
        grid on
        xlim(xli);  ylabel('Temperature [deg C]')
        legend('Pre-cal','post-cal','Location','Best')
    end

    % 6.2 CONDUCTIVITY
    group = "Updates";
    pref = "Conversion";
    boxtitle = "Offset";
    quest = ["Does an offset need to be applied to part of the conductivity series?"];
    pbtns = ["Yes","No"];
    [pval,~] = uigetpref(group,pref,boxtitle,quest,pbtns,opts);


    if strcmpi(pval,'yes')
          answer = inputdlg({'Enter interval where correction is required e.g. 200 230 400 480',...
                                'Enter offset value'},...
                     'Enter offset', [1 50; 1 50],{' ',' '},opts);
          invaloffc=str2num(answer{1});
          invalCc=str2num(answer{2});
      while 1
        if isodd(length(invaloffc))
          answer = inputdlg({'not even! - Enter interval where correction is required e.g. 200 230 400 480',...
                                'Enter offset value'},...
                     'Enter offset', [1 50; 1 50],{' ',' '},opts);
          invaloffc=str2num(answer{1});
          invalCc=str2num(answer{2});
        else
          break
        end
      end 
  
        for i = 1:length(invaloffc)/2
            dumI = find(jd0>= invaloffc(i*2-1) & jd0 <= invaloffc(i*2));    
            tn(dumI) = tn(dumI)+invalCc(i);
        end
    
        figure(1)
        axes(ax4)
        hold off
        ii = find(cn>dum);
        plot(jd0(ii),c(ii),'k')
        hold on
        plot(jd0(ii),cn(ii),'r')
        grid on
        xlim(xli);    ylabel('Conductivity [mS/cm]')
        legend('Pre-cal','post-cal','Location','Best')
        title('Post-calibration Salinity');
        
        % update salinity tiume series
        subplot(4,3,3)
        hold on
          if strcmp(acc_cpcor,'y')
           sn = gsw_SP_from_R(cn(ii)/c3515,tn(ii)*t90_68,pref*ones(length(ii),1));
          else
           sn = gsw_SP_from_R(cn(ii)/c3515,tn(ii)*t90_68,pn(ii));
          end
         plot(jd(ii)-jd(1),sn,'r')
         grid on
    end
   
    % 6.3 PRESSURE
    group = "Updates";
    pref = "Conversion";
    boxtitle = "Offset";
    quest = ["Does an offset need to be applied to part of the pressure series?"];
    pbtns = ["Yes","No"];
    [pval,~] = uigetpref(group,pref,boxtitle,quest,pbtns,opts);

    if strcmpi(pval,'yes')
          answer = inputdlg({'Enter interval where correction is required e.g. 200 230 400 480',...
                                'Enter offset value'},...
                     'Enter offset', [1 50; 1 50],{' ',' '},opts);
          invaloffp=str2num(answer{1});
          invalPc=str2num(answer{2});
      while 1
        if isodd(length(invaloffp))
          answer = inputdlg({'not even! - Enter interval where correction is required e.g. 200 230 400 480',...
                                'Enter offset value'},...
                     'Enter offset', [1 50; 1 50],{' ',' '},opts);
          invaloffp=str2num(answer{1});
          invalPc=str2num(answer{2});
        else
          break
        end
      end 
  
        for i = 1:length(invaloffp)/2
            dumI = find(jd0>= invaloffp(i*2-1) & jd0 <= invaloffp(i*2));    
            tn(dumI) = tn(dumI)+invalPc(i);
        end
    
        figure(1)
        axes(ax3)
        hold off
        ii = find(pn>dum);
        plot(jd0(ii),p(ii),'k')
        hold on
        plot(jd0(ii),pn(ii),'g')
        grid on
        xlim(xli);    ylabel('Pressure [db]')
        legend('Pre-cal','post-cal','Location','Best')      
        

    end
  
%% 7 interactive despiking
  cn(cn==dum)=NaN;
  pn(pn==dum)=NaN;
  tn(tn==dum)=NaN;
  spidx = ~isnan(cn) & ~isnan(pn);
  % this was creating spurious data which may be wiped. Data are only
  % indexed if pressure and cond are dum so if pressure or cond alone are
  % dum they are included in the equation below, so the salinity become
  % zero. Or if any temperatures are dummy then the salinity also becomes
  % zero. Switch to nans?

  % convert to salinity for plotting purpose only
  if sum(spidx) == 0 % if not pressure sensor, use reference pressure
    idx = ~isnan(cn);
    sn = gsw_SP_from_R(cn(idx)/c3515,tn(idx)*t90_68,pref*ones(length(idx),1));
  else    
    sn = gsw_SP_from_R(cn(spidx)/c3515,tn(spidx)*t90_68,pn(spidx));
  end 
  
  val=find(tn(spidx));
  % plot T-S 
  figure(1); axes(ax5);
  plot(sn,tn(val),'.r');grid on; hold on
  title('5. De-spike: Non-Nan index of salinity and pressure');
  xlabel('Salinity'); ylabel('In-situ Temperature');
  
  % CTD comparison/calibration
  if exist('diag')==1
      fclose(diag);
      open('CTD_near.txt')
  end
  
  % infinite loop 
  ELIM = [];
  while 1   
    answer = questdlg('\fontsize{16}Do you wish to remove spikes interactively?', ...
    'Spike removal', ...
    'Yes','No',opts);
    % Handle response
    switch answer
        case 'Yes'
          [x,y] = ginput;
          x = [x ;x(1)];
          y = [y ;y(1)];
          elim = find(inpolygon(sn,tn(val),x,y) == 1);
          ELIM = [ELIM;elim];
          cn(val(ELIM)) = NaN;
          sn(ELIM)      = NaN;
          figure(1);axes(ax5) ; hold off;
          plot(sn,tn(val),'.r')
          grid on;
          title('5. De-spike: Non-Nan index of salinity and pressure');
          xlabel('Salinity'); ylabel('In-situ Temperature');
        case 'No'
          figure(1)
          axes(ax2)
          cn(val(ELIM)) = NaN;
          sn(ELIM)      = NaN;
          ii = find(~isnan(cn));
          plot(jd0(ii),cn(ii),'r')
          grid on
          xlim(xli);  ylabel('Conductivity [mS/cm]')
          legend('Pre-cal','post-cal','Location','Best')
          fprintf(1,'%d data points eliminated \n\n\n',length(ELIM)) 
          break
    end  
  end
  
%% 8 automatic de-spiking
    val   = find(cn> dum &pn>dum);
    if isempty(val)
        val = find(cn>dum);
        sn  = gsw_SP_from_R(cn(val)/c3515,tn(val)*t90_68,pref*ones(length(val),1));
    else    
        sn =  gsw_SP_from_R(cn(val)/c3515,tn(val)*t90_68,pn(val));
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
        axes(ax5)
        hold on
        plot(smd+6*ssd,Tgrid,'m--','Linewidth',2)
        plot(smd-6*ssd,Tgrid,'m--','Linewidth',2)
        answer = questdlg('\fontsize{16}Do you want to exclude all the values outside the 6 sigma area?', ...
        'Spike removal', ...
        'Yes','No',opts);    
        if strcmpi(answer,'yes')
            disp(['Automatic despiking accepted - ',num2str(length(ELIM2)),' values discarded']);  
            cn(val(ELIM2)) = NaN;
            sn(ELIM)      = NaN;
        elseif strcmpi(answer,'no')
            disp('Automatic despiking rejected')    
            ELIM2 = [];
        end  
        figure(1);axes(ax5) ; hold off;
        plot(sn,tn(val),'.r')
        grid on;
        text('Units', 'Normalized', 'Position', [0.05, 0.05],...
            'string',['Automatic despiking accepted - ',num2str(length(ELIM2)),' values discarded'])   ;     
        text('Units', 'Normalized', 'Position', [0.05, 0.1],...
            'string',['Manual despike ' num2str(length(ELIM)) ' data points eliminated']);
        title('5. De-spike: Non-Nan index of salinity and pressure');
        xlabel('Salinity'); ylabel('In-situ Temperature');
    end
    
    % replace any NaNs with dummy
    ii     = isnan(pn);pn(ii) = dum;
    ii     = isnan(cn);cn(ii) = dum;
    ii     = isnan(tn);tn(ii) = dum;

%% 10 save data and figure 
 print(figure(1),[mcfig_out, ''],'-dpng');
 % print(figure(1),[mcfig_out, '.fig']);
 saveas(figure(1),[mcfig_out, '.fig']); % 14/01/2025

if skipT == 1 && skipC == 1 && skipP == 1    
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

 if strcmp(p_exist,'n')
     fprintf(fidtxt,'Conductivity pressure correction redone? %s\n',acc_cpcor);
     if strcmp(acc_cpcor,'y')
      fprintf(fidtxt,'Reference pressure for conductivity pressure correction: %s\n',num2str(pref));
      fprintf(fidtxt,'Time interval to which correction was applied: [%s]\n',num2str(invalcpcor));
     end
 end

 fprintf(fidtxt,'Skipped conductivity intervals: [%s]\n',num2str(invalC));
 fprintf(fidtxt,'Skipped temperature intervals : [%s]\n',num2str(invalT));
 fprintf(fidtxt,'Skipped pressure intervals    : [%s]\n',num2str(invalP));

 if exist('invaloffc','var')
     for i=1:length(invaloffc)/2
        fprintf(fidtxt,'Offset of %5.4f applied to conductivity interval [%s]\n',invalCc(i),num2str(invaloffc(i*2-1:i*2)));
     end
 end

 if exist('invalofft','var')
     for i=1:length(invalofft)/2
        fprintf(fidtxt,'Offset of %5.4f applied to temperature interval [%s]\n',invalTc(i),num2str(invalofft(i*2-1:i*2)));
     end
 end
 if exist('invaloffp','var')
     for i=1:length(invaloffp)/2
        fprintf(fidtxt,'Offset of %5.4f applied to pressure interval [%s]\n',invalPc(i),num2str(invaloffp(i*2-1:i*2)));
     end
 end
 fprintf(fidtxt,'Number of additional C points skipped interactively: [%d]\n',length(ELIM));
 fprintf(fidtxt,'Number of additional C points skipped automatically: [%d]\n',length(ELIM2));
 if mc<length(serial)
    answer = questdlg(sprintf('%s',['\fontsize{16}Instrument ' num2str(serial(mc)) ' processed. Move to ' num2str(serial(mc+1)) '? ']), ...
        'Process', ...
        'Yes','No, stop here',opts);
    switch answer
       case 'Yes'
            continue
       case 'No, stop here'
            break
    end
 end
end

fclose all; close('all');