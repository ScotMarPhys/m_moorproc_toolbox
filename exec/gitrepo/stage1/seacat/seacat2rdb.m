% Save raw ascii data from Seabird SBE16 'Seagauge' to rdb format
%
% 
%
% Kanzow, 16.04.05
% jym 22 April 2005: Start with clean slate


 clear all
 close all

% input - likely to be changed 
 mooring       = ['mochael_4_393'];       %
 cruise        = 'en517';
 operator      = 'dr400';
  
 % The path of the offsetfile is still preliminary

 inpath        = ['/noc/mpoc/rpdmoc/rapid/data/moor/raw/',cruise,'/seagauge/'];
 outpath       = ['/noc/mpoc/rpdmoc/rapid/data/moor/proc/',mooring,'/seacat/'];
 offsetfile    = ['/noc/mpoc/rpdmoc/rapid/data/moor/raw/',cruise,'/clock_offset.dat'];
 info  = ['/noc/mpoc/rpdmoc/rapid/data/moor/proc/',mooring,'/',mooring,'info.dat'];
 

 
 %%wrapped_sensor = [388   389  390   391   394   395   396   397   414  ];
 %%wrap_correction= [47.7  47.7 47.7  47.7  47.7  79.4  79.4  79.4  47.7  ]; 

 sg_code        = 332;
 toffset        = 0;  % time offset

% load infodat
 [lat,lon,wd,sdate,stime,edate,etime,z,type,serial,mocha]= ...
 rodbload(info,...
'Latitude:Longitude:WaterDepth:StartDate:StartTime:EndDate:EndTime:z:instrument:serialnumber:mochaname');


% ----------------------------------------

 sgI          = find(type==sg_code);
 sn           = serial(sgI);

 % NB The seagauge dir in outpath must be created first

 logfile       = [outpath,'stage1_log'];
 log           = fopen(logfile,'w');
 fprintf(log,'Processing of SBE16s  from mooring %s \n Date: %s \n',mooring,datestr(clock)');
 fprintf(log,'Operator: %s\n',operator');



for ins = 1 : length(sgI)

  serialnumber = sn(ins);

  instrdepth   = z(sgI(ins));

  outfile      = [outpath,mooring,'_',sprintf('%4.5d',serialnumber),'.raw'];

% ------ data input file and log file

  %file         = [inpath,sprintf('%4.4d',serialnumber),'_data.tid'];
  if ~isempty(mocha)
    file         = [inpath,sprintf('%d',mocha(sgI(ins))),'J.CNV'];
  else
    file         = [inpath,sprintf('m%s-01',mooring(end-2:end)),'.cnv']; % setup specifically for BPR/Seacat files received from Miami in July 2013. DR.
  end  
  disp(['SBE-16 instrument ',num2str(serialnumber),' has been found'])

  fprintf(log,'Serial number: %4.4d\n',serialnumber);
  fprintf(log,'Target Depth [m]: %4.4d\n',instrdepth);
  fprintf(log,'Input file: %s\n',file);
  fprintf(log,'Output file: %s\n',outfile);


% load data and carry out basic checks

  [P,T,C,S,jd,meas,sampling_rate,qflag] = load_seacat(file,log);

% ----------------------------------------------------------------
% -------check data ----------------
% ----------------------------------------------------------------


  % ---- check for potential time offset entry in the 'offsetfile' -----


   try 
     fidoff = fopen(offsetfile,'r');
   catch
   
     fprintf(fidlog,['Time offset file NOT found, \n NO check for ',...
                  'potential offset in recorded time applied \n'],offsetfile) 
     fprintf(1,['Time offset file  NOT found, \n NO check for ',...
              potential offset in recorded time applied \n'],offsetfile) 
   end   
  
   if fidoff > 0
           
     while 1 
       zeile = fgetl(fidoff);

       if ~ischar(zeile) 
          break 
       end 

       if ~isempty(findstr(zeile,mooring))   
          val = findstr(zeile,' ');
	  val = str2num(zeile(val:end));

          if val(1)== sg_code & val(2)== serialnumber 
	    toffset = val(3) + val(4)/24;      
            fprintf(1,['Instrument has time offset of %8.4f days ',...
                       '(rel. GMT) \n'],toffset); 
            fprintf(log,['Instrument has time offset of %8.4f days',...
                       ' (rel. GMT) \n'],toffset);
            break 
          end
       end
     end   % while

   end  % fidlog

  jd = jd - toffset; %apply potential correction of recorded time

  bottomstart     = julian([sdate(:)' hms2h([ stime(:)' 0])]);
  bottomstop      = julian([edate(:)' hms2h([ etime(:)' 0])]);

  if jd(end) < bottomstop
    fprintf(1,'\n\n W A R N I N G: Record already ends before end of mooring deployment period: \n Recorded time could be wrong!! \n')
  end 
  if jd(1) > bottomstart
    fprintf(1,'\n\n W A R N I N G: Record only starts after beginning of mooring deployment period: \n Recorded time could be wrong!! \n')
  end 


  % --- check if data if wrapped ----------------------------------------
if 0 
  wrap =  find(serialnumber == wrapped_sensor);
  figure(99);clf
  plot(P)
  grid on 
  title(['Serial number: ',num2str(serialnumber),'; target depth:',num2str(instrdepth)])
  hold on
  disp(['The bottom pressure may contain "wrapped" data - inspect figure 99'])
  corr= input(['Do you want to correct for "wrapped" data? y/n '],'s');

  %%  if ~isempty(wrap)& 
  if strcmp(corr,'y') 
    if ~isempty(wrap) 
      offset = 2^19 / wrap_correction(wrap)/1.4503774;
    else
      wrap_correction = input(['Insert pressure coefficient M from the ',...
		     'hex data file of sensor #',num2str(serialnumber),': ']);
      offset = 2^19 / wrap_correction/1.4503774;
    end 

    ind = input(['From figure 99 select index range of wrapped',...
                ' data [ind1 ind2]: ']); 
    Pneu = P;
    Pneu(ind(1):ind(end)) = Pneu(ind(1):ind(end)) + offset;
% jym 22 April 2005: Change plottin slightly to ensure both curves are 
% easily visible
    clf,     
    plot(Pneu,'r')
    hold on
    plot(P,'b')
    grid on
    legend('P_{wrapped}','P_{corrected}')
    success = input('Was correction successful? y/n ','s');

    if strcmp(success,'n')
        disp('Break')
        break
    else 
       P = Pneu; 
       fprintf(log,[['Correction for "wrapped" data applied for elements'], ...
            [' : [%d - %d] \n']],ind);
    end      
  end 
  
%% end
end % if 0

% ----- plot ---------------

 jd_start = julian([sdate' hms2h([stime;0]')']);
 jd_end   = julian([edate' hms2h([etime;0]')']);


 ii = find(jd<jd_end & jd>jd_start+.3);
 figure(1)
 plot(jd(ii)-jd(1),P(ii))

% --------- make logfile entries ----------------

   sz = size(P);
   fprintf(log,'Median Pressure [dbar]: %5.1f\n',median(P));
   fprintf(log,'Median Temperature [deg C]: %5.1f\n',median(T));
   fprintf(log,'Start date and time: %s \n',datestr(gregorian(jd(1))));
   fprintf(log,'End date and time:   %s \n',datestr(gregorian(jd(end))));
   sampling_rate = round(1./median(diff(jd)));
   ex_samples = round((jd(end)-jd(1))*sampling_rate+1);
   fprintf(log,'Sampling Frequency [per day]: %d \n',sampling_rate);
   fprintf(log,'Number of samples: %d; expected: %d \n\n\n',sz(1),ex_samples);
   if toffset ~= 0
     fprintf(log,'Offset of %8.4f days has been subtracted from recored time \n',toffset);
     fprintf(1,'Offset of %8.4f days has been subtracted from recored time \n',toffset);
   end
  if jd(end) < bottomstop
  fprintf(log,['\n\n W A R N I N G: Record already ends before end',...
  'of mooring deployment period: \n Recorded time could be wrong!! \n'])
  end 
  if jd(1) > bottomstart
  fprintf(log,['\n\n W A R N I N G: Record only starts after beginning',...
  'of mooring deployment period: \n Recorded time could be wrong!! \n'])
  end 

% ------save data (in rodb) ---------


  out = save_seacat(outfile,P,T,C,jd,mooring,lat,lon,wd,sdate,stime,edate,etime,instrdepth,serialnumber);

 fclose(log);



end
