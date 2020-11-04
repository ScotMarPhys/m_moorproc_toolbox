% basic preprocessing for combined argonaut sd / microcat data
%   
% features
%      1. eliminate lauching and recovery period
%      2. save data to rodb file
%      3. create data overview sheet
%
% uses timeaxis.m, auto_filt.m, julian.m 

% 14.04.05 Kanzow 
%
% 25.04.05 Kanzow; debugged     

function argocat_raw2use
  

% --- get moring information from infofile 

moor          = 'eb3_1_200407';    % Mooring name
operator      = 'Kanzow';

plot_interval = [2004 02 01 00;   % start time of time axis on plot
		 2005 06 01 00];  % end time of time axis on plot

in_code       = [366337] ;         % argocat id


% -- set path for data input and output -----------
inpath   = ['/local/users/pstar/data/moor/proc/',moor,'/arg/'];
outpath  = ['/local/users/pstar/data/moor/proc/',moor,'/arg/'];
infofile = ['/local/users/pstar/data/moor/proc/',moor,'/',moor,'info.dat'];

%%inpath   = ['/data/rapid/cd170/moorings/',moor,'/arg/'];
%%outpath  = ['/data/rapid/cd170/moorings/',moor,'/arg/'];
%%infofile = ['/data/rapid/cd170/moorings/',moor,'/',moor,'info.dat'];

% --- load infodat ----------------------------------

[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');


combo    =  find(id == 366337);
combo_sn =  sn(combo);
combo_z  =  z(combo);


fid_stat= fopen([outpath,'stage2_log'],'w');
fprintf(fid_stat,'Processing steps taken by argocat_raw2use.m:\n');
fprintf(fid_stat,'  1. eliminate lauching and recovery period\n');

fprintf(fid_stat,'  2. save data to rdb file\n');
fprintf(fid_stat,'\n Operated by:%s on %s\n',operator,datestr(clock)); 
fprintf(fid_stat,['        Argonaut/MicroCATs in Mooring ',moor,'\n\n\n']);

dummy    = -9999;

columns = 'YY:MM:DD:HH:T:TCAT:P:PCAT:C:U:V:W:HDG:PIT:ROL:USD:VSD:WSD:USS:VSS:WSS:HDGSD:PITSD:ROLSD:IPOW';


%-----------------------------------------
% --- preprocessing loop -------------------
% ----------------------------------------


jd_s  = julian(s_d(1),s_d(2),s_d(3),s_t(1)+s_t(2)/60);  % start time
jd_e  = julian(e_d(1),e_d(2),e_d(3),e_t(1)+e_t(2)/60);  % end time

for proc = 1 : length(combo_sn),


  indep  = combo_z(proc);

  infile  = [inpath,moor,'_',sprintf('%4.4d',combo_sn(proc)),'.raw'];
  if exist(infile)   > 0 
 
    rodbfile= [moor,'_',sprintf('%4.4d',combo_sn(proc)),'.use']; 
    outfile = [outpath,rodbfile];

    fprintf(fid_stat,'Serialnumber %d \n',combo_sn(proc));
    fprintf(fid_stat,'Infile %s \n',infile);
    fprintf(fid_stat,'Outfile %s \n',outfile);

  
    [YY,MM,DD,HH,t,tc,p,pc,c,u,v,w,...
     hdg,pit,rol,usd,vsd,wsd,uss,vss,wss,hdgsd,pitsd,rolsd,ipow,...
     serial_mc,serial_arg] = ...
     rodbload(infile,[columns,':MicrocatSN:ArgonautSN']);

    %------------------------------------------ 
    %----- cut off launching and recovery period
    %------------------------------------------
    disp('cut off launching and recovery period')
 
    jd               = julian(YY,MM,DD,HH);
    ii               = find(jd <= jd_e & jd >= jd_s );

    YY=YY(ii);MM=MM(ii);DD=DD(ii);HH=HH(ii);
    c=c(ii);t=t(ii);tc=tc(ii);p=p(ii);pc=pc(ii);
    u=u(ii);v=v(ii);w=w(ii);

    hdg=hdg(ii);pit=pit(ii);rol=rol(ii);usd=usd(ii);vsd=vsd(ii);wsd=wsd(ii);
    uss=uss(ii);vss=vss(ii);wss=wss(ii);hdgsd=hdgsd(ii);pitsd=pitsd(ii);
    rolsd=rolsd(ii);ipow=ipow(ii);

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

if 0

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
         
    fort =['%4.4d %2.2d %2.2d  %6.4f    %7.4f %7.4f %6.2f %6.2f %7.4f  ',...
          ' %6.2f %6.2f %6.2f   %4.1f %4.1f %4.1f   %6.2f %6.2f %6.2f  ',...
          ' %3.3d %3.3d %3.3d   %4.1f %4.1f %4.1f   %4.2f'];

      
    rodbsave(outfile,...
          'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time:MicrocatSN:ArgonautSN',...
           fort,...
          lat,lon,columns,Start_Date,Start_Time,combo_sn(proc),mr,wd,indep,End_Date,End_Time,serial_mc,serial_arg,...
         [ YY,MM,DD,HH,t,tc,p,pc,c,u,v,w,...
     hdg,pit,rol,usd,vsd,wsd,uss,vss,wss,hdgsd,pitsd,rolsd,ipow]);
   

  %%%%%%%%%% Graphics %%%%%%%%%%%%%%%%

    jd0 = julian(-1,1,1,24);
    jd1 = julian(plot_interval(1,:))-jd0;
    jd2 = julian(plot_interval(2,:))-jd0; 

      
    sampling_rate = 1/median(diff(jd));


    STR = ['Temperature [deg C]   ';
         'Pressure [dbar]       ';  
         'Conductivity [mS/cm]  ';
         'Zonal Velocity [cm/s] ';
         'Merid. Velocity [cm/s]'];
    VAR1= ['t';'p';'c';'u';'v']; 
    VAR2= ['tc';'pc'];        

    figure(1);clf

    for sub = 1 : 5
      eval(['var1 = ',VAR1(sub),';'])
      if sub < 3
        eval(['var2 = ',VAR2(sub,:),';'])
      else
        var2 = [];
      end  
    ok = plot_timeseries(jd,var1,var2,sampling_rate,STR(sub,:),sub,[jd1 jd2],'n');
    end

    subplot(5,1,1)
       title(['ArgoCAT s/n: ',num2str(combo_sn(proc)), ...
             '; Target Depth: ',num2str(indep)])
  
       orient tall
    eval(['print -depsc ',outfile,'.1.eps']) 


    figure(2);clf

    for sub = 1 : 5
      eval(['var1 = ',VAR1(sub),';'])
      if sub < 3
        eval(['var2 = ',VAR2(sub,:),';'])
      else
        var2 = [];
      end  
    ok = plot_timeseries(jd,var1,var2,sampling_rate,STR(sub,:),sub,[jd1 jd2],'y');
    end

    subplot(5,1,1)
      title(['ArgoCAT s/n: ',num2str(combo_sn(proc)), ...
             '; Target Depth: ',num2str(indep)])
  
      orient tall
    eval(['print -depsc ',outfile,'.2.eps']) 

  end  %proc
end 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function ok = plot_timeseries(jd,var1,var2,sr,str,sub,jdlim,filt)
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

  subplot(5,1,sub);
         
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

