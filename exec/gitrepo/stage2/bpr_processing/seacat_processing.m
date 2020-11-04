% Process data from Seabird SBE16 'Seagauge'
%
% uses  purge_bp.m 
%

% Kanzow, 24.04.06

clear 

% input - likely to be changed 

 mooring       = ['mochabl_4_394'];       %
 cruise        = 'en517';
 basepath      = '/Volumes/rpdmoc/rapid/data/moor/';
 operator      = 'dr400';
 plot_interval = [2009 12 01 00;   % start time of time axis on plot
	       	 2011 05 01 00];  % end time of time axis on plot

 
 inpath        = [basepath '/proc/',mooring,'/seacat/'];
 outpath       = [basepath '/proc/',mooring,'/seacat/'];
 clock_file    = [basepath '/raw/',cruise, '/seagauge/bpr_clock_offset.dat'];
 info          = [basepath '/proc/',mooring,'/',mooring,'info.dat'];

%% inpath        = ['/local/users/pstar/data/moor/proc/',mooring,'/seagauge/'];
%% outpath       = ['/local/users/pstar/data/moor/proc/',mooring,'/seagauge/'];
%% clock_file    = ['/local/users/pstar/data/moor/raw/seagauge/bpr_clock_offset.dat'];
%% info  = ['/local/users/pstar/data/moor/proc/',mooring,'/',mooring,'info.dat'];


 programme     = 'seacat_processing.m';
  
 
 sg_code       = 332;  % seacat identifier
 dum           = -9999; % dummy value
 
% --- load infodat -----------

 [lat,lon,wd,sdate,stime,edate,etime,z,type,serial]= ...
  rodbload(info,...
  'Latitude:Longitude:WaterDepth:StartDate:StartTime:EndDate:EndTime:z:instrument:serialnumber');


% ----------------------------------------

 sgI           = find(type==sg_code);
 sn            = serial(sgI);

 logfile       = [outpath,'stage2_log'];
 log           = fopen(logfile,'w');
 fprintf(log,'Stage 2 processing of SBE16 \n data  from mooring %s \n Date: %s \n',mooring,datestr(clock)');
 fprintf(log,'Operator: %s\n',operator');
 fprintf(log,'Programme: %s\n',programme);


for ins = 1 : length(sgI)

  serialnumber = sn(ins);
  instrdepth   = z(sgI(ins));
  outfile      = [outpath,mooring,'_',sprintf('%4.4d',serialnumber),'.use'];

 
% ------ data input file and log file

  file          = [inpath,mooring,'_',sprintf('%5.5d',serialnumber),'.raw'];

 disp(['SBE-16 instrument ',num2str(serialnumber),' has been found'])

 fprintf(log,'Serialnumber %d \n',serialnumber);
 fprintf(log,'Infile %s \n',file);
 fprintf(log,'Outfile %s \n',outfile);


% load data and carry out basic checks

 [yy,mm,dd,hh,P,T,C] = rodbload(file,'YY:MM:DD:HH:P:T:C');

 jd = julian(yy,mm,dd,hh);

% -------check data ----------------

  sampling_rate = round(1./median(diff(jd)));%nominal sampling rate [per day] 


 time_corr = input(['Shall time be corrected for clock offset? y/n '],'s');

  if strcmp(time_corr,'y')
    eval(['load ',clock_file])
    ins = find(serialnumber == bpr_clock_offset(:,2) & ...
               sg_code      == bpr_clock_offset(:,1));

    if isempty(ins)
      fprintf(1,['Clock offset not found in %s'],clock_file);
      fprintf(log,['Clock offset not found in %s'],clock_file);
      coff = input('Enter post deployment clock offset [s] from log sheet ');
    else 
       coff= bpr_clock_offset(ins,3);
    end
       
    disp(['Applying post recovery clock offset correction of ',...
           num2str(coff),' seconds (rel. UTC)']); 
    fprintf(log,'Post recovery clock offset of %d s subtracted \n',coff);
 

    len  = length(jd);
    coff = linspace(0,coff/86400,len); % assume linear clock drift

    jd   = jd - coff';  % subtract drift
  end

  jd_start = julian([sdate' hms2h([stime;0]')']);
  jd_end   = julian([edate' hms2h([etime;0]')']);


  [pr_int,pr_exfit,t_int,jd_grid,pr_linfit] =  ... 
          purge_bp(sampling_rate,jd,P,T,[jd_start jd_end],log);

   c_int = interp1(jd,C,jd_grid);
 
%------------ graphics ------------------------

  jd0 = julian(-1,1,1,24);
  jd1 = julian(plot_interval(1,:))-jd0;
  jd2 = julian(plot_interval(2,:))-jd0; 

  figure(21);clf

   subplot(3,1,1); ii = find(~isnan(pr_int) & pr_int~=dum);

   plot(jd_grid(ii)-jd0,pr_int(ii))
   title(['Seagauge s/n: ',num2str(serialnumber),...
          '; Target Depth: ',num2str(instrdepth)])
   ylabel('Pressure [dbar]')
   grid on
   xlim([jd1 jd2])
   %timeaxis(plot_interval(1,1:3));   
   datetick('x',12)

   subplot(3,1,2); ii = find(~isnan(t_int) & t_int~=dum);

   plot(jd_grid(ii)-jd0,t_int(ii))
   ylabel('Temperature [deg C]')
   grid on
   xlim([jd1 jd2])
%   timeaxis(plot_interval(1,1:3));   
   orient tall
   datetick('x',12)
  subplot(3,1,3); ii = find(~isnan(c_int) & c_int~=dum);

   plot(jd_grid(ii)-jd0,c_int(ii))
   ylabel('Cond [mS/cm]')
   grid on
   xlim([jd1 jd2])
%   timeaxis(plot_interval(1,1:3));   
   orient tall
   datetick('x',12)
   eval(['print -depsc ',outfile,'.1.eps']) 

  figure(22);clf     % 2 day low pass

   subplot(3,1,1); ii = find(~isnan(pr_int) & pr_int~=dum);
 
    plot(jd_grid(ii)-jd0,auto_filt(pr_int(ii),sampling_rate,1/2,'low',4))
    hold on
    plot(jd_grid-jd0,pr_exfit,'r')
    plot(jd_grid-jd0,pr_linfit,'g')
    legend('data','exp.-lin. fit','lin fit')
    title(['Seagauge s/n: ',num2str(serialnumber),...
	   '; Target Depth: ',num2str(instrdepth),' 2-day low pass'])
    ylabel('Pressure [dbar]')
    grid on
    xlim([jd1 jd2])
%    timeaxis(plot_interval(1,1:3));   
   datetick('x',12)


   %----- decide whether exponential or linear fit is more apropriate

   fit_select = input('Judge from Fig. 3 (or 22) which empircal fit should be stored: exponential-linear (1) or linear (2): ');

   if fit_select == 1
     fprintf(log,'Expontential-linear pressure fit to be stored with data \n');
     fprintf(1,'Expontential-linear pressure fit to be stored with data \n');
     pfit = pr_exfit;

   elseif fit_select == 2
     fprintf(log,'Linear pressure fit to be stored with data \n');
     fprintf(1,'Linear pressure fit to be stored with data \n');
     pfit = pr_linfit;

   end 
  % - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - 

   subplot(3,1,2); 
    
    plot(jd_grid(ii)-jd0,auto_filt(pr_int(ii),sampling_rate,1/2,'low',4)-pfit)
    title('Empirical drift estimate subtracted') 
    ylabel('Pressure [dbar]')
    grid on
    xlim([jd1 jd2])
%    timeaxis(plot_interval(1,1:3));   
   datetick('x',12)

   subplot(3,1,3); ii = find(~isnan(t_int) & t_int~=dum);

    plot(jd_grid(ii)-jd0,auto_filt(t_int(ii),sampling_rate,1/2,'low',4))
     
    ylabel('Temperature [deg C]')
    grid on
    xlim([jd1 jd2])
%    timeaxis(plot_interval(1,1:3));   
   datetick('x',12)

    orient tall

   eval(['print -depsc ',outfile,'.2.eps']) 



%----- save data (in rodb)------------

 cols      = 'YY:MM:DD:HH:P:T:C:PFIT';


 fort      = '%4.4d  %2.2d  %2.2d  %7.5f   %8.4f  %7.4f  %7.4f   %8.4f';
 time      = gregorian(jd_grid);
 data      = [time(:,1:3) hms2h(time(:,4:6)) pr_int(:) t_int(:) c_int(:) pfit(:)];


 rodbsave(outfile,...
       'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
         fort,...
         lat,lon,cols,sdate,stime(1:2),serialnumber,mooring,wd,instrdepth,edate,etime(1:2),...
         data);



 fclose(log);





end % end ins
