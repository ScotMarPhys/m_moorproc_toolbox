function [p_grid,jd_grid,tide_fit]=bottom_pressure_grid(moor,infopath,bp_id)
%function [p_grid,jd_grid,tide_fit]=bottom_pressure_grid(moor,infopath,bp_id)
%
% final stage of bottom pressure processing
% Data is 2 day low pass filtered and interpolated onto 12 hour grid
% A fit for removing fortnightly and monthly tides is  
%
% INPUT:
% moor      -- mooring name
% infopath  -- mooring directory, where infodat is located 
%              (default /rpdmoc/rapid/data/moor/proc/)
% bp_id     -- identifier for pressure sensor (default 465 = seagauge) 
%              470 = ixsbpr; 460 = wlr; 480 = Bourdon; 316 = PIES;
%              332 = seacat;
%
% OUTPUT:
% jd_grid   -- 
% p_grid    --
% tide_fit  --
%

% modified by ZBS on 30/1/2009:
%  - recognition of PIES data was done by including the 'bp_id =
%  316' option and direction to the proper subdirectory for stage2 PIES
%  BP data
%  - an exception for stage-2 input files was made for an
%  irregularity in naming seagauge files for the mochabl or mochael
%  moorings
  
  

if nargin <2
  infopath = '/rpdmoc/rapid/data/moor/proc/';
end

if nargin < 3
  bp_id = 465;
end  

if isempty(infopath)
  infopath = '/rpdmoc/rapid/data/moor/proc/';
end 

cutoff  = 1/2;  % cut off frequency [1/day]     

infofile = [infopath moor '/' moor 'info.dat'];

[sn,ins] = rodbload(infofile,'serialnumber:instrument');    


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
  jd_grid  = [];
  p_grid   = [];
  tide_fit = []; 
  return
end  


for sensor = 1 :length(sn)

  SN = sn(sensor);

  infile  = [infopath moor '/' ext(2:end) '/' moor '_' sprintf('%4.4d',SN) '.use'];

  if ins(sensor) == 465
    % deal with perhaps improperly named mochabl_1_369 or
    % mochael_1_370 data files
    if ( ~isempty(strfind(moor,'mochabl')) | ~isempty(strfind(moor,'mochael'))) & ...
          ~exist(infile,'file')
      infile  = [infopath moor '/' ext(2:end) '/' moor sprintf('%4.4d',SN) '.use'];
    end
    infile  = [infopath moor '/' ext(2:end) '/' moor '_' sprintf('%4.4d',SN) '.use'];
  elseif ins(sensor) == 316
    % full resolution pressure PIES data after stage2 processing
    infile = [infopath moor '/' ext(2:end) '/' moor '_' sprintf('%3.3d',SN) '_p.use'];
  end

  if ~exist(infile)
    %warning('using *.raw data file because *.use does not exist!')
    %infile  = [infopath moor '/' ext(2:end) '/' moor '_' sprintf('%3.3d',sensor) '.raw']; 
  
    warning('*.use file does not exist - breaking bottom_pressure_grid.m!')
    break
    
  end   
  
  outfile = [infopath moor '/' ext(2:end) '/' moor '_' sprintf('%3.3d',sensor) ext];

  [yy,mm,dd,hh,P,pfit,T] = rodbload(infile,'yy:mm:dd:hh:p:pfit:t');
  P                      = P - pfit;
  jd                     = julian(yy,mm,dd,hh);


  sf                     = round(1/median(diff(jd)));

  tnnan                  = find(~isnan(T));
  pnnan                  = find(~isnan(P));
  Pfilt                  = auto_filt(P(pnnan),sf,cutoff,'low',4);
  Tfilt                  = auto_filt(T(tnnan),sf,cutoff,'low',4);

  jd_grid                = round(jd(pnnan(1))+2):1/(4*cutoff):round(jd(pnnan(end))-2); 

  p_grid                 = interp1(jd(pnnan),Pfilt,jd_grid);
  p_grid                 = p_grid - mean(p_grid);
  t_grid                 = interp1(jd(tnnan),Tfilt,jd_grid);

  jd_tide                = round(jd(1)+1):1/12:round(jd(end)-1);  

  p_tide                 = interp1(jd(pnnan),Pfilt,jd_tide);
  p_tide                 = p_tide - mean(p_tide);  

  % -- tide fit -----------------------  
  
  tide_freq           = [0.471521 0.544374 1.015895 1.098033]; %[degree/hour]
  
  [XX,XX,XX,fit]      = tid_anal(jd_tide(:),p_tide(:),tide_freq(:),0.001);

  tide_fit            = interp1(jd_tide(:),fit,jd_grid);                  

  
  % plot correction and save to file

  figure(100+sensor)
  subplot(2,1,1)
  hold off
  plot(jd_grid-jd_grid(1),p_grid)
  hold on
  plot(jd_grid-jd_grid(1),tide_fit,'r')
  title([moor,'   sn: ',num2str(SN)])
  grid on
  ylabel('pressure [dbar]')
  
  subplot(2,1,2)
  hold off
  plot(jd_grid-jd_grid(1),t_grid)
  grid on
  ylabel('temperature [C]')
  
  
  orient tall
  print(gcf,'-depsc',[outfile '.eps']) 

  
  %keyboard



  % ---- save data to rodb format ---------------
  
  header ='Mooring:SerialNumber:WaterDepth:InstrDepth:Start_Date:Start_Time:End_Date:End_Time:Latitude:Longitude:Columns';  


  fort      = '%4.4d  %2.2d  %2.2d  %3.1f   %8.4f  %7.4f   %8.4f';
  time      = gregorian(jd_grid);
  data      = [time(:,1:3) hms2h(time(:,4:6)) p_grid(:) t_grid(:) tide_fit(:)];
  cols      = 'YY:MM:DD:HH:P:T:PFIT'; 
  
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

end  