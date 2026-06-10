% Final stage of bottom pressure processing to be run as a script.
% Multiple moorings can be run by using (e.g.) bp_stage3_2011_08_31.m.
%
% Data is 2 day low pass filtered and interpolated onto 12 hour grid
% A fit for removing fortnightly and monthly tides is calculated.
%
% INPUT: 'moor' and 'type' will be used if they already exist in
% the workspace, otherwise set them in this script.
%

% modified by ZBS on 2011/09/01:
%  - converted to a script, but still maintained so that all bottom
%  pressure records that need processing can be run with
%  bp_stage3_2011_08_03.m
%  - now generates a log-file
%
% modified by ZBS on 2009/01/30:
%  - recognition of PIES data was done by including the 'bp_id =
%  316' option and direction to the proper subdirectory for stage2 PIES
%  BP data
%  - an exception for stage-2 input files was made for an
%  irregularity in naming seagauge files for the mochabl or mochael
%  moorings
  

% tags for operator
cruise        = 'rb1201';
operator      = 'skye';

% assign which mooring to do, if not already specified by a calling script
if ~exist('moor','var') || isempty(moor)

  % the official mooring name/directory (e.g. 'eb1l5_5_200824')
  moor = input('enter a mooring name: ','s');

  % the numerical instrument code for the bottom pressure sensors.
  % for seagauge, ixsbpr, wlr, seacat, bourdon, and pies, these are:
  bpr_ids = [465 470 460 332 480 316]; 

  while 1==1
    bp_id = input('enter the instrment code (e.g. 465 for seagauge): ');
  
    if ~any(bp_id==bpr_ids)
      error('incompatible instrument code, please enter again')
    else
      break
    end
  end
  
end


basedir = '/Users/hydrosea5/Desktop/RB1201/rapid/data/moor/proc/';
infofile = [basedir moor '/' moor 'info.dat'];
outpath    = [basedir moor  '/seagauge/'];

programme     = 'bottom_pressure_grid_rb1201.m';


cutoff  = 1/2;  % cut off frequency [1/day]     

[sn,ins,sdate,stime,edate,etime] = rodbload(infofile,'serialnumber:instrument:StartDate:StartTime:EndDate:EndTime');    

insI = find(ins == bp_id);
sn   = sn(insI);

if bp_id == 465
  ext = '.seagauge';
elseif bp_id == 470
  ext = '.ixsbpr';
elseif bp_id == 460
  ext = '.wlr';
elseif bp_id == 332
  ext = '.seacat';
elseif bp_id == 480
  ext = '.bourdon';   
elseif bp_id == 316
  ext = '.pies';
end  

if isempty(sn)
  disp('No bottom pressure sensor exists in this mooring')
end  


for sensor = 1 :length(sn)

  SN = sn(sensor);

  % find the infile and define the outfile
  
  infile  = [basedir moor '/' ext(2:end) '/' moor '_' sprintf('%4.4d',SN) '.use'];

  if ins(sensor) == 465
    % deal with perhaps improperly named mochabl_1_369 or
    % mochael_1_370 data files
    if ( ~isempty(strfind(moor,'mochabl')) | ~isempty(strfind(moor,'mochael'))) & ...
          ~exist(infile,'file')
      infile  = [basedir moor '/' ext(2:end) '/' moor sprintf('%4.4d',SN) '.use'];
    end
    infile  = [basedir moor '/' ext(2:end) '/' moor '_' sprintf('%4.4d',SN) '.use'];
  elseif ins(sensor) == 316
    % full resolution pressure PIES data after stage2 processing
    infile = [basedir moor '/' ext(2:end) '/' moor '_' sprintf('%3.3d',SN) '_p.use'];
  end

  if ~exist(infile)
    %warning('using *.raw data file because *.use does not exist!')
    %infile  = [basedir moor '/' ext(2:end) '/' moor '_' sprintf('%3.3d',sensor) '.raw']; 
  
    %warning('*.use file does not exist - breaking bottom_pressure_grid2.m!')
    error('*.use file does not exist - breaking bottom_pressure_grid_rb1201.m!')
    break
    
  end   

  outfile = [basedir moor '/' ext(2:end) '/' moor '_' sprintf('%3.5d',SN) ext];


  % initialize the log file
  
  logfile = [outpath,'stage3_log'];
  disp(['instrument SN=' num2str(SN) ' of type=' ext(2:end) ' has been found'])

  disp(logfile);
  log  = fopen(logfile,'w');
  fprintf(log,'Stage 3 processing of BPR data from mooring %s \n',moor');
  fprintf(log,'Date: %s \n',datestr(clock)');
  fprintf(log,'Operator: %s \n',operator);
  fprintf(log,'Processing performed on/at: %s \n',cruise);
  fprintf(log,'Programme: %s \n',programme);
 
  fprintf(log,'Serialnumber %d \n',SN);
  fprintf(log,'Instrument code %d \n',bp_id);
  fprintf(log,'Extension %s \n',ext);
  
  fprintf(log,'Infile %s \n',infile);
  fprintf(log,'Outfile %s \n',outfile);


  
  % start the stage3 processing
  % this remove the trend and tide prediction; the residual is an anomaly
  [yy,mm,dd,hh,P,pfit,ptide1,ptide2,T] = rodbload(infile,'yy:mm:dd:hh:p:pfit:ptide:ptidem:t');
  P                      = P - pfit - ptide1;
  jd                     = julian(yy,mm,dd,hh);

  sf                     = round(1/median(diff(jd)));
  
  % not sure the following is very good if there are NaNs in the middle
  tnnan                  = find(~isnan(T));
  pnnan                  = find(~isnan(P));
  
  Pfilt                  = auto_filt(P(pnnan),sf,cutoff,'low',4);
  Tfilt                  = auto_filt(T(tnnan),sf,cutoff,'low',4);
    
  jd_grid                = round(jd(pnnan(1))+0.5/cutoff):1/(4*cutoff):round(jd(pnnan(end))-0.5/cutoff); 

  p_grid                 = interp1(jd(pnnan),Pfilt,jd_grid);
  t_grid                 = interp1(jd(tnnan),Tfilt,jd_grid);


  % plot  and save to file
  jd_start = julian([sdate' hms2h([stime;0]')']);
  jd_end   = julian([edate' hms2h([etime;0]')']);
  
  figure(100+sensor)
  subplot(2,1,1)
  hold off
  plot(datenum(gregorian(jd_grid)),p_grid,'k');
  hold on
  plot(datenum(gregorian(jd)),ptide2,'r');
  foo = char(moor);
  foo(foo=='_')='  ';
  title([foo,'   sn: ',num2str(SN)]);
  grid on
  ylabel({'pressure anomaly' ;'[dbar]'});
  set(gca,'xlim',[ datenum(gregorian(jd_start))-10 datenum(gregorian(jd_end))+10])
  datetick('x','mmmyy');
  set(gca,'ylim',max(abs(p_grid(:)'))*[-1.1 1.1]);
  set(gca,'xlim',[ datenum(gregorian(jd_start))-10 datenum(gregorian(jd_end))+10])
    
  subplot(2,1,2)
  hold off
  plot(datenum(gregorian(jd_grid)),t_grid,'k');
  grid on
  ylabel('temperature [C]');
  set(gca,'xlim',[ datenum(gregorian(jd_start))-10 datenum(gregorian(jd_end))+10]);
  datetick('x','mmmyy');
  set(gca,'xlim',[ datenum(gregorian(jd_start))-10 datenum(gregorian(jd_end))+10])
      
  orient tall
  print(gcf,'-depsc',[outfile '.eps']);

  

  % ---- save data to rodb format ---------------
  
  header ='Mooring:SerialNumber:WaterDepth:InstrDepth:Start_Date:Start_Time:End_Date:End_Time:Latitude:Longitude:Columns';  


  fort      = '%4.4d  %2.2d  %2.2d  %2.2f  %8.4f  %7.4f';
  time      = gregorian(jd_grid);
  data      = [time(:,1:3) hms2h(time(:,4:6)) p_grid(:) t_grid(:)];
  cols      = 'YY:MM:DD:HH:P:T'; 
  
  [a,b,c,d,e,f,g,h,k,l,m] = rodbload(infile,header);

  e = time(1,1:3);
  f = time(1,4:5);
  g = time(end,1:3);
  h = time(end,4:5);
  
  if exist(outfile)
    warning(['overwriting existing outfile: ' outfile])
  end
  rodbsave(outfile,header,fort,...
           moor,b,c,d,e,f,g,h,k,l,cols,data); 

  
  fprintf(log,'\n Saved to outfile');
  fclose(log);

end  
