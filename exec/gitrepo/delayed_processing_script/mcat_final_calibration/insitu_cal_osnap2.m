% compare Mooring sensors to CTD at bottle-stops for post-calibration
%
% uses read_rosfile.m, ctd_impact.m, mc_concorr.m, num_legend.m, pload.m (+subroutines) 
%
% This version can be used for MicroCAT / Sontek Argonaut data 
%
%  Input 
%    - simultaneous Mooring sensor and CTD data (+ corresponding info.dat and mooring information) 
%
%
%  Output
%
%  dt,dc,dp                   -- T / C / P  diff. MicroCAT - CTD at the bottlestops
%  dt_mcdep,dc_mcdep,dp_mcdep -- same as above, but interp. onto sensor deployment depths
%  dt_av,dc_av,dp_av          -- same quantity, but aver. over pressure range (average_interval)
%

% Kanzow, adapted from 'microcat_insitu_cal.m'. 
%         'microcat_insitu_cal.m' is obsolete and will not be seviced any more
%         03.01.06  Sontek / Argonaut data option added
% Loic Houpert,  - adapted the function read_rosfile to be able to read time in julian days 
% 19.10.2015,    - adapted the function read_rosfile and ctd_impact to deal with cnv and ros files 
%                 with a problem in the  header variable #start_time
%                - add reference to structure p_insitucal in which the parameters
%                 are defined (see osnap/users/loh/mcatpostcruisecalib/postcalib_process)  
%               - add a parameter mc_cunit
close all
clearvars  -except pathgitrepo pathosnap p_insitucal pathgit
warning off
global cruise

% ---- parameters specified in the main script ----------------------------------------------------

sensor_id            = p_insitucal.sensor_id;
cast                 = p_insitucal.cast;
sensorselec          = p_insitucal.sensorselec;
average_interval     = p_insitucal.average_interval;
p_interval           = p_insitucal.p_interval;
cruise               = p_insitucal.cruise;
depl_period          = p_insitucal.depl_period;
basedir              = p_insitucal.basedir;
datadir              = p_insitucal.datadir;
coef_dir             = p_insitucal.coef_dir;
apply_offset         = p_insitucal.apply_offset;
interval_move        = p_insitucal.interval_move;  
ctd_latestart_offset = p_insitucal.ctd_latestart_offset;
dp_tol               = p_insitucal.dp_tol;
c_interval           = p_insitucal.c_interval; 
t_interval           = p_insitucal.t_interval; 
dp_interval          = p_insitucal.dp_interval;  

% --- futhers parameters (unlikely to be changed) -------------------

cond_threshold     = p_insitucal.cond_threshold;
impact_var         = p_insitucal.impact_var;
bottlestop_tmin    = p_insitucal.bottlestop_tmin;
bottlestop_dpmin   = p_insitucal.bottlestop_dpmin;
bottlestop_average = p_insitucal.bottlestop_average;                  
cnv_time_correction=  p_insitucal.cnv_time_correction;
lat                = p_insitucal.lat;
                          
                          
% ------ output / input files and directories ---------------

if sensor_id(1) >= 332 & sensor_id(end) <= 337
   sensor      = 'microcat';
elseif sensor_id(1) == 366
   sensor      = 'arg';
   cond_threshold = 10; % dbar
   impact_var     = 'p'; 
   c_interval     = [-2 2]; % cond. dev.plot range
   t_interval     = [-2 2]; % temp. dev.plot range 
   dp_interval    = [-35  35];   % pres. dev. plot range 
elseif sensor_id(1) == 330
    sensor   =   'rbr';
elseif sensor_id(1) == 339
    sensor   =   'idr';
else
   disp('Sensor type unknown - return')
   return
end   

if sensorselec == 1
    output_name    = [sensor,'_insitu_cal_',cruise,'_',p_insitucal.depl_period,'_','cast',num2str(cast),'_cond1'];
elseif sensorselec == 2
    output_name    = [sensor,'_insitu_cal_',cruise,'_',p_insitucal.depl_period,'_','cast',num2str(cast),'_cond2'];
end    

% exemple of other cruises data directory and structure in: insitu_cal_osnap_sub_list_cruises

%Data Directories
if strcmp('kn221-02',cruise)
  mc_ext      = '.raw';  
  mc_dir      = [basedir '/data/moor/proc_calib/',cruise,'/cal_dip/',sensor,'/'];
  info_dir    = [basedir '/data/moor/proc_calib/',cruise,'/cal_dip/'];
  ctd_dir     = [basedir '/cruise_data/',cruise,'/CTD_FINAL/CTD/binavg_1Hz_mat/'];
  ctdraw_dir  = [basedir '/cruise_data/',cruise,'/CTD_FINAL/CTD/os1407_ctd_binavg_1Hz/'];
  btldir  = [basedir '/cruise_data/',cruise,'/ctd_uncalib/datcnv/'];  
  ctdformat   = 'aoml'; % aoml format
  ctd_cunit   = 'S/m'; % cond. unit in .cnv file
  ctd_1hz     = 'S/m'; % cond. unit in 1 hz file
  mc_cunit    = 'S/m';
elseif strcmp('kn221-03',cruise)
  mc_ext      = '.raw';  
  mc_dir      = [basedir '/data/moor/proc_calib/',cruise,'/cal_dip/',sensor,'/'];
  info_dir    = [basedir '/data/moor/proc_calib/',cruise,'/cal_dip/'];
  ctd_dir     = [basedir '/cruise_data/',cruise,'/1hz_bin_averaged_mat/'];
  ctdraw_dir  = [basedir '/cruise_data/',cruise,'/1hz_bin_averaged/'];
  btldir      = [basedir '/cruise_data/',cruise,'/ros_files/'];  
  ctdformat   = 'aoml'; % aoml format
  ctd_cunit   = 'S/m'; % cond. unit in .cnv file
  ctd_1hz     = 'S/m'; % cond. unit in 1 hz file
  mc_cunit    = 'mS/cm';    
elseif strcmp('pe399',cruise)
  mc_ext      = '.raw';  
  mc_dir      = [basedir '/data/moor/proc_calib/',cruise,'/cal_dip/',sensor,'/'];
  info_dir    = [basedir '/data/moor/proc_calib/',cruise,'/cal_dip/'];
  ctd_dir     = [basedir '/../../Cruises/PE399_2015/mstar/',cruise,'/data/ctd/'];
  ctdraw_dir  = [basedir '/../../Cruises/PE399_2015/mstar/',cruise,'/data/ctd/ASCII_FILES/'];
  btldir      = [basedir '/../../Cruises/PE399_2015/mstar/',cruise,'/data/ctd/ASCII_FILES/CTD_Processed_by_Sven/CTD-' sprintf('%2.2d',cast) filesep ];
  ctdformat   = 'mstar'; % NOC format
  ctd_cunit   = 'S/m'; % cond. unit in .cnv file
  ctd_1hz     = 'mS/cm'; % cond. unit in 1 hz file
  mc_cunit    = 'mS/cm';
elseif strcmp('pe400',cruise)
  mc_ext      = '.raw';  
  mc_dir      = [basedir '/data/moor/proc_calib/',cruise,'/cal_dip/',sensor,'/'];
  info_dir    = [basedir '/data/moor/proc_calib/',cruise,'/cal_dip/'];
  ctd_dir     = [basedir '/cruise_data/',cruise,'/1hz_bin_averaged/'];
  ctdraw_dir  = [basedir '/cruise_data/',cruise,'/1hz_bin_averaged_uncal/'];
  btldir      = [basedir '/cruise_data/',cruise,'/ros_files/'];  
  ctdformat   = 'aoml'; % aoml format
  ctd_cunit   = 'S/m'; % cond. unit in .cnv file
  ctd_1hz     = 'S/m'; % cond. unit in 1 hz file
  mc_cunit    = 'mS/cm';    
elseif strcmp('dy053',cruise)
  mc_ext      = '.raw';  
  mc_dir      = [basedir,'/data/moor/proc_calib/',cruise,'/cal_dip/',sensor,'/'];
  info_dir    = [basedir,'/data/moor/proc_calib/',cruise,'/cal_dip/'];
  ctd_dir     = [basedir '/cruise_data/',cruise,'/ctd/'];
  ctdraw_dir  = [basedir '/cruise_data/',cruise,'/ctd/ASCII_FILES/'];
  btldir      = [basedir '/cruise_data/',cruise,'/ctd/ASCII_FILES/'];  
  ctdformat   = 'mstar'; % NOC format
  ctd_cunit   = 'S/m'; % cond. unit in .cnv file
  ctd_1hz     = 'mS/cm'; % cond. unit in 1 hz file
  mc_cunit    = 'mS/cm';     
elseif strcmp('dy078',cruise)
  mc_ext      = '.raw';  
  mc_dir      = [basedir,'/data/moor/proc_calib/',cruise,'/cal_dip/',sensor,'/'];
  info_dir    = [basedir,'/data/moor/proc_calib/',cruise,'/cal_dip/'];
  ctd_dir     = [basedir '/cruise_data/',cruise,'/ctd/'];
  ctdraw_dir  = [basedir '/cruise_data/',cruise,'/ctd/ASCII_FILES/'];
  btldir      = [basedir '/cruise_data/',cruise,'/ctd/ASCII_FILES/'];  
  ctdformat   = 'mstar'; % NOC format
  ctd_cunit   = 'S/m'; % cond. unit in .cnv file
  ctd_1hz     = 'mS/cm'; % cond. unit in 1 hz file
  mc_cunit    = 'mS/cm';
elseif strcmp('ar304',cruise)
  mc_ext      = '.raw';  
  mc_dir      = [basedir,'/data/moor/proc_calib/',cruise,'/cal_dip/',sensor,'/'];
  info_dir    = [basedir,'/data/moor/proc_calib/',cruise,'/cal_dip/'];
  ctd_dir     = [basedir '/cruise_data/',cruise,'/ctd/'];
  ctdraw_dir  = [basedir '/cruise_data/',cruise,'/ctd/ASCII_FILES/'];
  btldir      = [basedir '/cruise_data/',cruise,'/ctd/ASCII_FILES/'];  
  ctdformat   = 'mstar'; % NOC format
  ctd_cunit   = 'S/m'; % cond. unit in .cnv file
  ctd_1hz     = 'mS/cm'; % cond. unit in 1 hz file
  mc_cunit    = 'mS/cm';  
elseif strcmp('dy120',cruise)
  mc_ext      = '.raw';  
  mc_dir      = [basedir,'/data/moor/proc_calib/',cruise,'/cal_dip/',sensor,'/'];
  info_dir    = [basedir,'/data/moor/proc_calib/',cruise,'/cal_dip/'];
  ctd_dir     = [basedir '/cruise_data/',cruise,'/ctd/'];
  ctdraw_dir  = [basedir '/cruise_data/',cruise,'/ctd/ASCII_FILES/'];
  btldir      = [basedir '/cruise_data/',cruise,'/ctd/ASCII_FILES/'];  
  ctdformat   = 'mstar'; % NOC format
  ctd_cunit   = 'S/m'; % cond. unit in .cnv file
  ctd_1hz     = 'mS/cm'; % cond. unit in 1 hz file
  mc_cunit    = 'mS/cm';  
else
  error('Cruise unknown')  
end


% Deployment Depths File Directory 
ein_dir     = coef_dir ; 

%Data Files
if strcmp('kn221-02',cruise)
  ctd_file    = ['ctd_cal_kn221_allsta_1Hz.mat'] ; % ['ab1104_ctd_caldip_1Hz.mat'];
  bottle_file = ['OS1407_',sprintf('%3.3d',cast),'.ros'];
  ctd_cnvfile = ['OS1407_',sprintf('%3.3d',cast),'.cnv'];
  mc_file     = ['cast',num2str(cast),'/cast',num2str(cast),'_'];
  info_file   = ['cast',num2str(cast),'info.dat'];
elseif strcmp('pe399',cruise)
  ctd_file    = ['ctd_pe399_',sprintf('%3.3d',cast),'_1hz.nc'];
  dd=dir([btldir '/*.ros']);
  bottle_file = dd.name;%['PE399_',sprintf('%3.3d',cast),'.ros']; % .ros
  ctd_cnvfile = ['PE399_',sprintf('%3.3d',cast),'_align_ctm.cnv'];
  mc_file     = ['cast',num2str(cast),'/cast',num2str(cast),'_'];
  info_file   = ['cast',num2str(cast),'info.dat'];
elseif strcmp('pe400',cruise)
  ctd_file    = 'PE400_1Hz.mat' ; % ['ab1104_ctd_caldip_1Hz.mat'];
  bottle_file = ['PE400_',num2str(cast),'.ros'];
  ctd_cnvfile = ['PE400_',num2str(cast),'.cnv'];
  mc_file     = ['cast',num2str(cast),'/cast',num2str(cast),'_'];
  info_file   = ['cast',num2str(cast),'info.dat'];  
elseif strcmp('kn221-03',cruise)
  ctd_file    = 'ctd_cal_kn221-03_1Hz.mat' ; % ['ab1104_ctd_caldip_1Hz.mat'];
  bottle_file = ['kn221-03_',sprintf('%3.3d',cast),'.ros'];
  ctd_cnvfile = ['kn221-03_',sprintf('%3.3d',cast),'.cnv'];
  mc_file     = ['cast',num2str(cast),'/cast',num2str(cast),'_'];
  info_file   = ['cast',num2str(cast),'info.dat'];    
elseif strcmp('dy053',cruise)
  ctd_file    = ['ctd_dy053_',sprintf('%3.3d',cast),'_1hz.nc'];
  bottle_file = ['DY053_',sprintf('%3.3d',cast),'.ros']; % .ros
  ctd_cnvfile = ['DY053_',sprintf('%3.3d',cast),'_alignDO_actm.cnv'];
  mc_file     = ['cast',num2str(cast),'/cast',num2str(cast),'_'];
  info_file   = ['cast',num2str(cast),'info.dat'];  
elseif strcmp('dy078',cruise)
  ctd_file    = ['ctd_dy078_',sprintf('%3.3d',cast),'_1hz.nc'];
  bottle_file = ['CTD',sprintf('%3.3d',cast),'.ros']; % .ros
  ctd_cnvfile = ['CTD',sprintf('%3.3d',cast),'_align_actm.cnv'];
  mc_file     = ['cast',num2str(cast),'/cast',num2str(cast),'_'];
  info_file   = ['cast',num2str(cast),'info.dat'];    
elseif strcmp('ar304',cruise)
  ctd_file    = ['ctd_ar304_',sprintf('%3.3d',cast),'_raw.nc']; % post-calibrated ctd data generated by Adam and read in mstar (ar304new)
  bottle_file = ['ar30-04',sprintf('%3.3d',cast),'.ros']; % .ros
  ctd_cnvfile = ['ar30-04',sprintf('%3.3d',cast),'.cnv'];
  mc_file     = ['cast',num2str(cast),'/cast',num2str(cast),'_'];
  info_file   = ['cast',num2str(cast),'info.dat'];     
elseif strcmp('dy120',cruise)
  ctd_file    = ['ctd_dy120_',sprintf('%3.3d',cast),'_psal.nc']; % post-calibrated data from crusie dy120, processed by Kristin Burmeister
  bottle_file = ['DY120_',sprintf('%3.3d',cast),'.ros']; % .ros
  ctd_cnvfile = ['DY120_',sprintf('%3.3d',cast),'.cnv'];
  mc_file     = ['cast',num2str(cast),'/cast',num2str(cast),'_'];
  info_file   = ['cast',num2str(cast),'info.dat'];   
end  

ctd_file    = [ctd_dir,ctd_file];
bottle_file = [btldir,bottle_file];
ctd_cnvfile = [ctdraw_dir,ctd_cnvfile];
mc_file     = [mc_dir,mc_file];     
info_file   = [info_dir,info_file];
 
%-------------------------------------
% 1. --- load data  ------------------
%-------------------------------------

% ----MicroCAT ID -----------

[xxx,typ,instr,deploy] = rodbload(info_file,'z:instrument:serialnumber:deployment');
val           = find(typ>=sensor_id(1) & typ<=sensor_id(end));
instr         = instr(val);

						
% --- deployment depths --------

ein           = [ein_dir,depl_period,'_deploymentdepths.dat'];
[dep,typ,ssn] = rodbload(ein,'z:instrument:serialnumber');
if isempty(dep)
    error(['!! MAKE sure the file ' ein ' exists or edit the filepath !!' ])    
    return
end		   
% ---- CTD ----------------

fprintf(1,'\n Loading CTD data ...\n')
if strcmp(ctdformat,'aoml')
  eval(['load ',ctd_file])
  disp('aoml')
  if strcmp('sj0614',cruise)
      cnv_cor_save = cnv_cal;
  end
  if strcmp('rb0701',cruise)
      cnv_cor_save = cnv_cal;
  end
  if  strcmp('kn221-02',cruise) | strcmp('kn221-03',cruise)
    cnv_cor_save=sbe;
  end
   if  strcmp('pe400',cruise)
    cnv_cor_save=ctd_cal;
  end 
  
  for i = 1: length(cnv_cor_save)
    if strcmp('pe400',cruise)
    	strf = regexp(cnv_cor_save(i).file,'_','split');
    	if ~isempty(strfind(strf{2},num2str(cast)))
           d.cond  = cnv_cor_save(i).C; % better sensor
           d.temp  = cnv_cor_save(i).T;
           d.press = cnv_cor_save(i).P;   
           d.time =  cnv_cor_save(i).JD0 + julian([2015 1 1 0 0 0]) - 1;
    	end   
    else  
      if cnv_cor_save(i).station == cast         
          if strcmp('kn221-02',cruise) | strcmp('kn221-03',cruise)
          d.cond  = cnv_cor_save(i).c1S; % better sensor
          d.temp  = cnv_cor_save(i).t190C;
          d.press = cnv_cor_save(i).prDM;
          
          ctd_time_ori = julian([cnv_cor_save(i).gtime(1:3) 0]);
          jul_day_frac = cnv_cor_save(i).timeJ - floor(cnv_cor_save(i).timeJ);
          d.time = jul_day_frac + ctd_time_ori;                     
          else
           d.cond  = cnv_cor_save(i).conductivity;
          d.temp  = cnv_cor_save(i).temperature;
          d.press = cnv_cor_save(i).pressure;
          %if strcmp(ctd_cunit,'S/m')
          d.cond = d.cond;%*10;             
          
          ctd_time_ori = julian([cnv_cor_save(i).gtime(1:3) hms2h(cnv_cor_save(i).gtime(4:6))]);
          d.time = cnv_cor_save(i).elap_time_sec/86400 + ctd_time_ori;                           
          end
      end
    end
  end
  
  clear cnv_cor_save
  pack
elseif strcmp(ctdformat,'pstar')  
    ctd_file
  [d h]  = pload(ctd_file,'press temp cond time','silent');

  year   = floor(h.iymd/10000);
  month  = floor((h.iymd - year*10000)/100);
  day    = h.iymd -year*10000 - month*100;
  hour   = h.ihms;
  
  ctd_time_ori = julian(year+h.icent, month, day,hour);
  d.time       = d.time/86400 + ctd_time_ori; % ctd time in julian days
elseif strcmp(ctdformat,'mstar')  
    if strcmp(cruise,'d344') || strcmp(cruise,'d359')
        ctd_file
  dd  = netcdf(ctd_file); %,'press temp cond time','silent');
  d = struct('press',dd{'press'}(:),'temp',dd{'temp'}(:),'cond',dd{'cond'}(:),'time',dd{'time'}(:));
  year   = dd.data_time_origin(1); %floor(h.iymd/10000);
  month  = dd.data_time_origin(2); %floor((h.iymd - year*10000)/100);
  day    = dd.data_time_origin(3); %h.iymd -year*10000 - month*100;
  hour   = dd.data_time_origin(4) + (dd.data_time_origin(5)+(dd.data_time_origin(6)/60))/60 ; %h.ihms;
  
  ctd_time_ori = julian(year,month,day,hour); %julian(year+h.icent, month, day,hour);
  d.time       = d.time/86400 + ctd_time_ori; % ctd time in julian days 
% if strcmp(cruise,'jc064')
elseif strcmp(cruise,'dy120') | strcmp(cruise,'ar304') | strcmp(cruise,'dy078') | strcmp(cruise,'dy053') | strcmp(cruise,'pe399') 
    ctd_file
    if sensorselec==1    
        [d h]=mload(ctd_file,'time','press','temp1','cond1',' ','q');
        d.cond = d.cond1;
        d.temp = d.temp1;
    elseif sensorselec == 2
        [d h]=mload(ctd_file,'time','press','temp2','cond2',' ','q');
        d.cond = d.cond2;
        d.temp = d.temp2;        
    end
 HH = h.data_time_origin(4)+h.data_time_origin(5)/60+h.data_time_origin(6)/3600;
 jtime=h.data_time_origin(1:3);
 jtime=[jtime HH];
 d.time=julian(jtime)+d.time/86400;
else
    ctd_file
  dd  = netcdf(ctd_file); %,'press temp cond time','silent');
  d = struct('press',dd{'press'}(:),'temp',dd{'temp1'}(:),'cond',dd{'cond1'}(:),'time',dd{'time'}(:));
  year   = dd.data_time_origin(1); %floor(h.iymd/10000);
  month  = dd.data_time_origin(2); %floor((h.iymd - year*10000)/100);
  day    = dd.data_time_origin(3); %h.iymd -year*10000 - month*100;
  hour   = dd.data_time_origin(4) + (dd.data_time_origin(5)+(dd.data_time_origin(6)/60))/60 ; %h.ihms;
  
  ctd_time_ori = julian(year,month,day,hour); %julian(year+h.icent, month, day,hour);
  d.time       = d.time/86400 + ctd_time_ori; % ctd time in julian days
 
  end
end   
% --- CTD Bottle -----------
fprintf(1,'\n')

if findstr(bottle_file,'.btl')
    bottle = read_botfile(bottle_file);
    for jjj=1:length(bottle)
        bottle(jjj).p=bottle(jjj).pav;
    end
    bottle.jd = [julian(bottle.yy,bottle.mm,bottle.dd,bottle.hh)]';
elseif findstr(bottle_file,'.ros')
    bottle            = read_rosfile(bottle_file);
else
    bottle = [];
    disp('PATH TO BOTTLE FILES NOT DEFINED')
end

    bottle.jd         = bottle.jd - cnv_time_correction;
    bottle.start_time = bottle.start_time -cnv_time_correction;

% --- MicroCAT -------

fprintf(1,'\n loading MicroCAT/CTD data ... \n\n')

%instr=instr(13:15);
ninst = length(instr);

for mc = 1 : ninst

  [yy,mm,dd,hh,p,t,c]  = rodbload([mc_file,sprintf('%4.4d',instr(mc)),mc_ext],...
                                   'YY:MM:DD:HH:P:T:C'); 
  lyy                = length(yy);
  YY(1:lyy,mc)       = yy; 
  MM(1:lyy,mc)       = mm; 
  DD(1:lyy,mc)       = dd; 
  HH(1:lyy,mc)       = hh; 
  T(1:lyy,mc)        = t; 
  C(1:lyy,mc)        = c; 
  P(1:lyy,mc)        = p; 
end

ii    = find(YY == 0);
YY(ii)=NaN; MM(ii)=NaN; DD(ii)=NaN; HH(ii)=NaN; T(ii)=NaN; C(ii)=NaN;
JD    = julian(YY,MM,DD,HH);                    % julian day mc time 

% check which variables have been measured by sensor  

pstat = 1;
tstat = 1;
cstat = 1;

if isempty(find(~isnan(P)))
   pstat = 0;
end   
if isempty(find(~isnan(T)))
   tstat = 0;
end   
if isempty(find(~isnan(C)))
   cstat = 0;
end   

% --------------------------------------------------------------------------
% 2. ---- water impact times: determine time offsets between ctd and mc ----
% --------------------------------------------------------------------------
% CTD impact
ctd_cond_threshold = cond_threshold;
if strcmp(ctd_cunit,'S/m')
    ctd_cond_threshold = ctd_cond_threshold/10;
end    
[wit_ctd,dwit_ctd] = ctd_impact(ctd_cnvfile,impact_var,ctd_cond_threshold,cruise)% water impact time ctd
wit_ctd(4)         = wit_ctd(4) - cnv_time_correction*24;

if size(wit_ctd,2) == 6
  wit_ctd_jd = julian([wit_ctd(1:3) hms2h(wit_ctd(4:6))]);
elseif  size(wit_ctd,2) == 4
  wit_ctd_jd = julian(wit_ctd);
end

% MicroCAT impact
if strcmp(mc_cunit,'S/m')
    C = C*10;
end    

wit_mc_jd = nan(1,ninst);
for inst = 1 : ninst
  if strcmp(impact_var,'c') 
    ii = find(C(:,inst) > cond_threshold);           % water impact mc
  elseif strcmp(impact_var,'p')
     ii = find(P(:,inst) > cond_threshold);           % water impact mc    
  end 
  if isnan(nanmean(C(:,inst)))
      disp(' '); 
      error(['NO DATA FOR INSTRUMENT ' num2str(instr(inst)) ...
          '. MAKE SURE THAT A FILE cast' num2str(p_insitucal.cast) '_' num2str(instr(inst)) ...
          '.raw EXISTS IN cal_dip/microcat/cast' num2str(p_insitucal.cast) ...
          '/  OR REMOVE THE INSTRUMENT ' num2str(instr(inst))  ' FROM THE FILE ' ....
          'moor/proc_calib/' p_insitucal.cruise '/cal_dip/cast' num2str(p_insitucal.cast) 'info.dat' ]) 
  end
  ii = ii(1); 
  wit_mc_jd(inst)  = julian([YY(ii,inst),MM(ii,inst),DD(ii,inst),HH(ii,inst)]);
end

wit_mc_jd    =  wit_mc_jd'; 
wit_mc       = gregorian(wit_mc_jd);  % 
start_mc_jd  = julian([YY(1,:)' MM(1,:)' DD(1,:)' HH(1,:)']); 


% time difference between start of instrument and ocean surface impact [s]

dwit_mc   = (wit_mc_jd - start_mc_jd)*24*3600;

% time offset  between CTD and MC: needed to compare bottle stop values

ii        = find(dwit_mc > 60);   % only consider mc with impact time > 60 for others
                                  % may have started after surface impact

impact_offset =  wit_ctd_jd - wit_mc_jd; % individual impact time offsets ctd - mc


if ~isempty(ii)
  offset    = wit_ctd_jd -  median(wit_mc_jd(ii));         % odecimal days
  [ohms(1),ohms(2),ohms(3)] = s2hms(offset*24*3600);
else
  offset = 0;
  [ohms(1),ohms(2),ohms(3)] = s2hms(offset*24*3600);  
end

[oh,om,os] = s2hms(impact_offset*86400);
fprintf(1,'\n Time offset of MicroCATs rel. to CTD:\n\n') 
fprintf(1,'  ID    HH    MM   SS\n')
fprintf(1,'  %d: %d  %d  %d\n',[instr';oh';om';round(os')]);

%-------------------------------------------------------
% 3. ---- extract data from bottle stops ---------------
%-------------------------------------------------------


bot_start     = bottle.jd;
% % istop = find(diff(bot_start)
% bot_end       = bot_start + bottlestop_average/3600/24;

% correct ctd bottle stop time for offset

if apply_offset == 'y'
  bot_start = bot_start - offset;
  bot_end   = bot_end - offset;
  fprintf(1,' O F F S E T   H A S   B E E N   A P P L I E D  ! ! !\n')
elseif  apply_offset == 'n'
  fprintf(1,' O F F S E T   H A S   N O T   B E E N   A P P L I E D  ! ! !\n')
elseif  apply_offset == 'i'
  fprintf(1,' D A N G E R !!! INDIV.  OFFSETS  HAVE  BEEN   APPLIED  ! ! !\n')
  JD  = JD + ones(size(JD,1),1)*impact_offset' ;
end



bst=1:length(bot_start);
% to identify single value for each bottle stop
bst2     = find(diff(bottle.p(bst))<-bottlestop_dpmin);
bst2     = [1; bst2+1];
bst      = bst(bst2);


% % bot_end  = bot_end(bst); 
% % bot_start= bot_start(bst);
bot_start = bottle.jd(bst);
% bot_end = bottle.jd([bst-1 ;length(bottle.jd)]);

nstop = length(bot_start);     % number of bottle stops 
gcnt     = 0;

% ---- MicroCAT data ------------------------

%bot_start = bot_start + offset/86400;
 %  bot_end = bot_end + offset/86400;

interval_move1=interval_move/3600/24;
mcatbotstart0=nan(nstop,ninst);
mcatbotend0=nan(nstop,ninst);
mcatbotstart=nan(nstop,ninst);
mcatbotend=nan(nstop,ninst);
for stop = 1 : nstop,    % bottle_stops loop
   for mc = 1 : ninst,     % instrument loop
 
        dp=gradient(P(:,mc));
        dcond = gradient(C(:,mc));
            
            [~,indstop] = nearest(bot_start(stop),JD(:,mc));
            presstop = P(indstop,mc);
%            imcatbotok00 = find(abs(dp)<1 & JD(:,mc)>bot_start(stop)+interval_move1(1) & JD(:,mc)<bot_start(stop)+interval_move1(2));% & abs(dcond)<0.02 );    
            imcatbotok00 = find(P(:,mc)>presstop-3 & P(:,mc)<presstop+3 & JD(:,mc)>bot_start(stop)+interval_move1(1) & JD(:,mc)<bot_start(stop)+interval_move1(2));% & abs(dcond)<0.02 );               
%             %--------------------------------------------------------------
%             % find the longest consecutive index
%             diffimact = diff(imcatbotok00);
%             if length(imcatbotok00)<2 | isempty(find(diffimact==1)),continue;end           
%             cc=1;
%             iindm =cell(length(find(diffimact~=1))+1,1);
%             for ijk=1:length(diffimact)
%                 if diffimact(ijk)==1
%                    iindm{cc} = [  iindm{cc} imcatbotok00(ijk)];
%                 else
%                     cc = cc +1;
%                 end
%             end
% 
%             llmax = 0;
%             csel = 0;
%             for illm = 1:length(iindm)
%                 if llmax < length(iindm{illm})
%                     llmax = length(iindm{illm});   
%                     csel = illm;
%                 end
%             end
%             imcatbotok = iindm{csel};
%             %---------------------------------------------------------------

% Add a condition to remove the first  30sec of the bottle stop.
            jdtime0 = JD(imcatbotok00(1),mc);
            imcatbotok = imcatbotok00(find(JD(imcatbotok00,mc)>jdtime0 + 0.5/24/60));
% And check that the length of the bottlestop is at least == to bottlestop_tmin (in sec)

            if ~isempty(imcatbotok) & ((JD(imcatbotok(end),mc) - JD(imcatbotok(1),mc))*3600*24>bottlestop_tmin)

                mcatbotstart0(stop,mc)=JD(imcatbotok(1),mc); 
                mcatbotend0(stop,mc)=JD(imcatbotok(end),mc);    
            else
                 mcatbotstart0(stop,mc) = nan;
                 mcatbotend0(stop,mc) = nan;   
            end
   end
   
   diffmcatbot = mcatbotend0(stop,:)-mcatbotstart0(stop,:);
   imcsel = find(diffmcatbot==min(diffmcatbot),1);
   
   for mc = 1 : ninst,  
       if isempty(mcatbotstart0(stop,imcsel))
                  mcatbotstart(stop,mc) = nan;
                 mcatbotend(stop,mc) = nan;   
       else
            mcatbotstart(stop,mc) = mcatbotstart0(stop,imcsel);
            mcatbotend(stop,mc)   = mcatbotend0(stop,imcsel);
       end
   end
end

figure; plot((JD(:,:)-JD(1,1))*24*60,P(:,:))

for stop = 1 : nstop,    % bottle_stops loop
  
  figure(10+stop); clf; hold on; ooo = 200/60/60/24;

   for mc = 1 : ninst,   

       ii_move = find(JD(:,mc)<=mcatbotend(stop,mc) & ...
                      JD(:,mc)>=mcatbotstart(stop,mc));
       ii_move2 = find(JD(:,mc)<=(mcatbotend(stop,mc)+2*ooo) & ...
                      JD(:,mc)>=(mcatbotstart(stop,mc)-2*ooo));                 
       tav(stop,mc) = mean(T(ii_move,mc));
       cav(stop,mc) = mean(C(ii_move,mc));
       pav(stop,mc) = mean(P(ii_move,mc)); 
       if cstat == 1 
         plot((JD(ii_move2,mc)-JD(1,1))*24*60,C(ii_move2,mc),'k') 
         plot((JD(ii_move,mc)-JD(1,1))*24*60,C(ii_move,mc),'b') 
       else
         %plot((JD(ii_move2,mc)-JD(1,1))*24*60,T(ii_move2,mc),'k') 
         %plot((JD(ii_move,mc)-JD(1,1))*24*60,T(ii_move,mc),'b') 
         plot((JD(ii_move2,mc)-JD(1,1))*24*60,P(ii_move2,mc),'k') 
         plot((JD(ii_move,mc)-JD(1,1))*24*60,P(ii_move,mc),'b')          
       end    
   end
   title([cruise,'   cast',num2str(cast),'  depth: ',num2str(round(bottle.p(bst(stop))))])
end

[llllnn mmmmnn]=size(JD);
if mmmmnn>4
ii_move=find(JD(:,5)<=mcatbotend(stop,5) & ...
                      JD(:,5)>=mcatbotstart(stop,5));
Ttemp=T(ii_move,5);
Ctemp=C(ii_move,5);
else

end


% ------  extract bottlestop values from ctd 

if rms(interval_move) ~= 0
  if  apply_offset   == 'y'
    ctd_time =  - offset + d.time; 
  else 
    ctd_time =  d.time;
  end
if ctd_latestart_offset ~=0
    ctd_time = ctd_time + ctd_latestart_offset/86400;
    fprintf(1,'CTD LATESTART OFFSET HAS BEEN APPLIED !!!\n')
end

ctdbotstart = nan(nstop,1);
ctdbotend = nan(nstop,1);
  for stop = 1 : nstop, 
      
     ctdbotstart(stop) = min(mcatbotstart(stop,:));
     ctdbotend(stop) = max(mcatbotend(stop,:));     
      
     ii_move = find(ctd_time<=ctdbotend(stop)...
                    & ctd_time>=ctdbotstart(stop));
     ii_move2 = find(ctd_time<=ctdbotend(stop)+ooo...
                    & ctd_time>= ctdbotstart(stop)-2*ooo);
     figure(10+stop)
     if cstat == 1
         if strcmp('d334',cruise)
       plot((ctd_time(ii_move2)-JD(1,1))*24*60,d.cond(ii_move2)*10,'r') 
       plot((ctd_time(ii_move)-JD(1,1))*24*60,d.cond(ii_move)*10,'b') 
         elseif exist('ctd_1hz')==1 && strcmp(ctd_1hz,'S/m')
             plot((ctd_time(ii_move2)-JD(1,1))*24*60,d.cond(ii_move2)*10,'r') 
       plot((ctd_time(ii_move)-JD(1,1))*24*60,d.cond(ii_move)*10,'b')
         else
        plot((ctd_time(ii_move2)-JD(1,1))*24*60,d.cond(ii_move2),'r') 
       plot((ctd_time(ii_move)-JD(1,1))*24*60,d.cond(ii_move),'b') 
         end
     else
       plot((ctd_time(ii_move2)-JD(1,1))*24*60,d.temp(ii_move2),'r') 
       plot((ctd_time(ii_move)-JD(1,1))*24*60,d.temp(ii_move),'b') 
     end    
     ctd_pav(stop) = mean(d.press(ii_move));
     ctd_tav(stop) = mean(d.temp(ii_move));
     if strcmp('d334',cruise)
         ctd_cav(stop) = mean(d.cond(ii_move)*10); % pstar units are S/m and not mS/cm like normal
     elseif exist('ctd_1hz')==1 && strcmp(ctd_1hz,'S/m')
         ctd_cav(stop) = mean(d.cond(ii_move)*10); % cases where 1hz data in S/m and not mS/cm
     else
     ctd_cav(stop) = mean(d.cond(ii_move)); 
     end
     grid on 
  end
end
% 
% [iuna iunb iunc]=unique(ctd_pav);
% ctd_pav = ctd_pav(iunb);
% ctd_tav = ctd_tav(iunb);
% ctd_cav = ctd_cav(iunb);

MP = max(d.press) +100;

% --------------------------------------------------------------------
% ----------------- CALCULATE DEVIATIONS MC - CTD  ------------------- 
% --------------------------------------------------------------------

% -------- allocate sensors to deployment depths ----

for i = 1 : ninst,
   IN      = find((typ>=sensor_id(1) & typ<=sensor_id(end)) & ssn == instr(i));
   if isempty(IN) 
     mcdep2(i) = NaN;
     mcdep(i) = NaN;
   else 
     mcdep2(i) = dep(IN);
     mcdep(i) = sw_pres(mcdep2(i),lat);
   end
end


% --- compute dc, dt, dp

  bot_p0av    =  ctd_pav'*ones(1,ninst);
  dp          = pav - (bot_p0av) ;   
                                                      % dp = P_mc - P_ctd  
  pproblem    = find(rms(dp) > dp_tol);
                                                      
  bot_t0av    = ctd_tav'*ones(1,ninst);
  dt          = tav - (bot_t0av);   
                                                    % dt = T_mc - T_ctd 
  kompx       = find(isnan(P(1,:)));   % index of sensors without own pressure sensor
  cavc        = mc_concorr(cav(:,kompx),bot_p0av(:,kompx),0); 
                                                    % compressibility corrected C
  cav(:,kompx)= cavc;      % insert corrected C    

  cav_problem = mc_concorr(cav(:,pproblem),bot_p0av(:,pproblem),pav(:,pproblem));
  
  bot_c0av    = ctd_cav'*ones(1,ninst);
  dc          = cav - (bot_c0av );  
  dc_pproblem  = cav_problem - bot_c0av(:,pproblem);
                                                      % dc = C_mc - C_ctd
  average_interval = average_interval(1):20:average_interval(2);
 
  % compute interpolated and averaged versions of dt,dc,dp 
  
  if max(average_interval) > max(d.press)
    fprintf(1,'WARNING: average_interval exceeds max. press. of CTD cast, please fix!!!\n')  
  end
  dp_mcdep_ext(1:ninst) = NaN; % dp extrap. beyond  max. cast pressure to deploym. depth
  dc_av_pproblem(1:ninst)= NaN; % dc_av for MicroCATs with problematic pressure records 

  pp_count   = 1; % pressure problem index

  for i = 1:ninst,
     oknan = find(~isnan(bot_p0av(:,1)));
    dt_mcdep(i) =  interp1([MP bot_p0av(oknan,1)'],...
                                [dt(1,i) dt(oknan,i)'],mcdep(i));
    dc_mcdep(i) =  interp1([MP bot_p0av(oknan,1)'],...
                                [dc(1,i) dc(oknan,i)'],mcdep(i));
    dp_mcdep(i) =  interp1([MP bot_p0av(oknan,1)'],...
                                [dp(1,i) dp(oknan,i)'],mcdep(i));
                            

    if isnan(dp_mcdep(i)) % extrapolate if deployment pressure > max. pressure of cast
        pol             = polyfit([bot_p0av(:,1)'],[dp(:,i)'],3);  
        dp_mcdep_ext(i) = polyval(pol,mcdep(i));
    end    
    
    dt_av(i)  =  mean(interp1([MP bot_p0av(oknan,1)'],...
                              [dt(1,i) dt(oknan,i)'],average_interval));
    dc_av(i)  =  mean(interp1([MP bot_p0av(oknan,1)'],...
                              [dc(1,i) dc(oknan,i)'],average_interval));
    dp_av(i)  =  mean(interp1([MP bot_p0av(oknan,1)'],...
                              [dp(1,i) dp(oknan,i)'],average_interval));

    dt_sd(i)  =  std(interp1([MP bot_p0av(oknan,1)'],...
                              [dt(1,i) dt(oknan,i)'],average_interval));
    dc_sd(i)  =  std(interp1([MP bot_p0av(oknan,1)'],...
                              [dc(1,i) dc(oknan,i)'],average_interval)); 
                          
    if ~isempty(find(i==pproblem)) % if MC press. bad, use ctd press. to correct conduct.
       dc_av_pproblem(i) = mean(interp1([MP bot_p0av(oknan,1)'],...
                    [dc_pproblem(1,pp_count) dc_pproblem(oknan,pp_count)'],average_interval));
       pp_count = pp_count + 1; 
    end    

  end

% ---- INTERNAL SUBROUTINES ----------------------

% ------ save ----------------------------------

instr_id = ones(nstop,1) * instr';
%%%eval(['save ',mc_dir,'cast',num2str(cast),'/cp_mc_bot.mat dt dc dp bot_pav pav bot_c0av cav bot_t0av tav dtemi dconi dprei apply_offset offset cast instr_id'])

% -------------------------------------------------
% ------ GRAPHICS --------------------------------
% -------------------------------------------------
 
% plot temp_deviations

col  = ['brgkmcybrgkmcybrgkmcybrgkmcy'];
mrk  = ['ddddddd+++++++xxxxxxxsssssss'];
lin  = ['----------------------------'];
mrkd = [col;mrk];
mrk  = [col;mrk;lin];

% ----------  Graphics  ----------------------
 
if ~isempty(find(~isnan(dp)))
    sub = 3; 
    figure(2);subplot(sub,1,1);hold off;subplot(sub,1,2);hold off;
             subplot(sub,1,3);hold off;hold off; 
else 
    sub = 2; 
    figure(2);subplot(sub,1,1);hold off;subplot(sub,1,2);hold off;hold off
end


for i = 1 : ninst,

  subplot(1,sub,1)
    plot(bot_p0av(:,i),dt(:,i),mrk(:,i))
    hold on

  subplot(1,sub,2)
    plot(bot_p0av(:,i),dc(:,i),mrk(:,i))
    hold on

  if ~isempty(find(~isnan(dp)))
      subplot(1,sub,3)
      plot(bot_p0av(:,i),dp(:,i),mrk(:,i))
      hold on
  end 
end

 %plot values at MicroCAT deployment depths

for i = 1 : ninst
    subplot(1,sub,1)
    plot(mcdep(i),dt_mcdep(i),[mrkd(:,i)],'Linewidth',2,'MarkerSize',11)
    subplot(1,sub,2)
    plot(mcdep(i),dc_mcdep(i),[mrkd(:,i)],'Linewidth',2,'MarkerSize',11)
   if ~isempty(find(~isnan(dp)))
      subplot(1,sub,3)
      plot(mcdep(i),dp_mcdep(i),[mrkd(:,i)],'Linewidth',2,'MarkerSize',11)  
      if isnan(dp_mcdep(i))
        plot(mcdep(i),dp_mcdep_ext(i),[mrkd(:,i)],'Linewidth',2,'MarkerSize',11)    
      end    
   end

end

subplot(1,sub,1)
   grid on   
   set(gca,'Fontsize',12,'xlim',p_interval,'ylim',t_interval)
   xlabel('pressure [dbar]')
   ylabel(['temp. diff. [C]'])
   title([sprintf('Deviations %s-CTD //%s',upper(sensor),date) ])
   num_legend(instr','''best''',5)  %legend 
   %num_legend(instr','southeast',5)  %legend
subplot(1,sub,2)
   grid on 
   set(gca,'Fontsize',12,'xlim',p_interval,'ylim',c_interval)
   xlabel('pressure [dbar]')
   ylabel('cond. diff. [mS/cm]') 
   title(['Cruise: ',cruise,'  Cast: ',num2str(cast)]) 

   if ~isempty(find(~isnan(dp)))
   subplot(1,sub,3)
   grid on 
   set(gca,'Fontsize',12,'xlim',p_interval,'ylim',dp_interval)
   xlabel('pressure [dbar]')
   ylabel('pres. diff. [dbar]') 
   end
  


set(figure(2),'Paperunits','centimeters','Paperposition',[0 0 29 21])
set(0,'currentfigure',figure(2))     

print([mc_dir,'cast',num2str(cast),'/',output_name,'.png'],'-dpng');

set(0,'currentfigure',figure(2))  
orient landscape

print([mc_dir,'cast',num2str(cast),'/',output_name,'.ps'],'-dpsc');
    

if ~isempty(pproblem) 
  figure(4);clf;hold on
  
  for i = 1 : length(pproblem),

    plot(bot_p0av(:,i),dc_pproblem(:,i),mrk(:,pproblem(i)))
    hold on
    %ylim(c_interval)
    xlabel('pressure [dbar]')  
    ylabel('cond. dev. [mS/cm]') 
    xlim(p_interval)
    num_legend(instr(pproblem),'''best''',5)
    grid on
    title('MicroCATs pressure problems: cond. diff. MC-CTD (corrected by CTD pressure)') 
    set(gca,'FontSize',12)
  end
    
end   
% -- TEXT OUTPUT ---------
if strcmp('rb0701',cruise)
    if cast == 27
fid = fopen([mc_dir,'cast',num2str(cast+1),'/',output_name,'.txt'],'w');
fprintf(fid,['%s',cruise,'  cast:',num2str(cast+1),' MicroCAT - CTD // processing date: ',date,' \n'],'%');
    elseif cast == 28
fid = fopen([mc_dir,'cast',num2str(cast-1),'/',output_name,'.txt'],'w');
fprintf(fid,['%s',cruise,'  cast:',num2str(cast-1),' MicroCAT - CTD // processing date: ',date,' \n'],'%');
    else
fid = fopen([mc_dir,'cast',num2str(cast),'/',output_name,'.txt'],'w');
fprintf(fid,['%s',cruise,'  cast:',num2str(cast),' MicroCAT - CTD // processing date: ',date,' \n'],'%');
    end
else
fid = fopen([mc_dir,'cast',num2str(cast),'/',output_name,'.txt'],'w');
fprintf(fid,['%s',cruise,'  cast:',num2str(cast),' MicroCAT - CTD // processing date: ',date,' \n'],'%');
end

fprintf(fid,'%sID   CAST   DEPTH RANGE  DEPTH   dt  dt_sd     dc  dc_sd       dp    |  dp_ext | dc_pproblem\n','%'); 
fprintf(fid,'%sID   CAST    T/C [m]     P [m]   [K]   [K] [mS/cm]   [ms/cm]  [dbar] |  [dbar] |   [mS/cm] \n','%'); 

caldata     = ...
              [instr cast*ones(ninst,1) ones(ninst,1)*average_interval([1 end])...
               mcdep2' dt_av' dt_sd' dc_av' dc_sd' dp_mcdep' dp_mcdep_ext' dc_av_pproblem'];

fprintf(fid,' %4.4d  %2.2d    %4.4d-%4.4d  %4.4d    %5.4f %5.4f   %5.4f  %5.4f   %3.1f  |  %3.1f  | %5.4f\n',caldata');   

fprintf(fid,['\n At nominal depth of instrument \n' ],'%');
fprintf(fid,'%sID   CAST   DEPTH   dt   dc     dp    |  dp_ext \n','%'); 
fprintf(fid,'%sID   CAST   P [m]   [K]  [mS/cm]  [dbar] |  [dbar] \n','%'); 

caldata2     = ...
              [instr cast*ones(ninst,1) ...
               mcdep2' dt_mcdep' dc_mcdep' dp_mcdep' dp_mcdep_ext' ];

fprintf(fid,' %4.4d  %2.2d   %4.4d    %5.4f    %5.4f  %3.1f  |  %3.1f  \n',caldata2');   
  
