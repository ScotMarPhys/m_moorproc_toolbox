% function rodb_to_oceansitesnetcdf(moor,procpath,moorinfo)
%
% Function for converting RODB mooring data to OceanSites netcdf format
%
% required inputs:-
%   moor: complete mooring name as string. e.g. 'wb1_1_200420'
%   procpath: can specify exact procpath if not using standard data paths. 
%           e.g. pressure_overlay('wb1_1_200420','inpath','/Volumes/noc/mpoc/hydro/rpdmoc/rapid/data/moor/proc/')
%   moorinfo: structure with metadata for the different type of instruments
%
% functions called:-
%   rodbload, julian, auto_filt
%   from .../exec/moor/tools and .../exec/moor/rodb paths
% 
% 18/11/16 - Loic Houpert
%
function rodb_to_oceansitesnetcdf(moor,procpath,moorinfo)
if nargin <1
    help rodb_to_oceansitesnetcdf
    return
end

if isunix
    infofile=[procpath,moor,'/',moor,'info.dat'];
elseif ispc
    infofile=[procpath,moor,'\',moor,'info.dat'];
end

ncdir = 'oceansites_format'; %outpudirectory for oceansites netcdf


% Load vectors of mooring information
% id instrument id, sn serial number, z nominal depth of each instrument
% s_t, e_t, s_d, e_d start and end times and dates
% lat lon mooring position, wd corrected water depth (m)
% mr mooring name
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,...
    'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');


% JULIAN Convert Gregorian date to Julian day.
% JD = JULIAN(YY,MM,DD,HH) or JD = JULIAN([YY,MM,DD,HH]) returns the Julian
% day number of the calendar date specified by year, month, day, and decimal
% hour.
% JD = JULIAN(YY,MM,DD) or JD = JULIAN([YY,MM,DD]) if decimal hour is absent,
% it is assumed to be zero.
% Although the formal definition holds that Julian days start and end at
% noon, here Julian days start and end at midnight. In this convention,
% Julian day 2440000 began at 00:00 hours, May 23, 1968.
jd_start = julian([s_d' hms2h([s_t;0]')]);
jd_end   = julian([e_d' hms2h([e_t;0]')]);

% find the index number of Microcats
iiMC = find(id == 337 | id == 335);
vecMC = sn(iiMC);
zMC = z(iiMC);
% find the index number of RBRs
iiRBR = find(id == 330);
vecRBR = sn(iiRBR);
% find the index number of Idronauts
iiIDR = find(id == 339);
vecIDR = sn(iiIDR);
% find the index number of S4s
iiS4 = find(id == 302);
vecS4 = sn(iiS4);
% and find index number of RCM11s
iiRCM11 = find(id == 310);
vecRCM11 = sn(iiRCM11);
% and find index number of Sontek Argonauts
iiARG = find(id == 366);
vecARG = sn(iiARG);
% and find index number of Nortek Aquadopps
iiNOR = find((id == 368|id==370));
vecNOR = sn(iiNOR);
zNOR = z(iiNOR);
% and find index number of AADI Seaguards
iiSG = find(id == 301);
vecSG = sn(iiSG);
% Possible ADCP codes - taken from IMP moorings package
iiADCP = find((id>=319) & (id <=328));
vecADCP = sn(iiADCP);


depths(:,1) = id([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG;iiADCP]);
depths(:,2) = z([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG;iiADCP]);
depths=sortrows(depths,2);
iiiMC=find(depths(:,1)==337 | depths(:,1)== 335);
iiiRBR=find(depths(:,1)==330);
iiiIDR=find(depths(:,1)==339);
iiiS4=find(depths(:,1)==302);  
iiiRCM11=find(depths(:,1)==310);  
iiiARG=find(depths(:,1)==366 | depths(:,1)==366337); 
iiiNOR=find(depths(:,1)==368|depths(:,1)==370);
iiiSG=find(depths(:,1)==301);
iiiADCP=find(depths(:,1)>=319 & depths(:,1)<=328);
iii=[iiiS4;iiiRCM11;iiiARG;iiiNOR;iiiMC;iiiRBR;iiiIDR;iiiSG;iiiADCP];

% -----------------------------------
% START OF READING IN INSTRUMENT DATA
% -----------------------------------
%--------------------------------------
% Now read in Microcat data if required
%--------------------------------------
if iiMC>0
    j=1;
    
    ibad = [];
    jbad = [];
    % loop to read one file at a time
    netcdffilepath = [ncdir '/' moor '_MCAT'];
    for i=1:length(vecMC);
        serialno = vecMC(i);
        disp('*************************************************************')
        disp(['Reading MICROCAT - ',num2str(serialno)])
        disp('*************************************************************')
        
        
        if isunix
            infile = [procpath,moor,'/microcat/',moor,'_',sprintf('%0.3d',i),'.microcat'];
        elseif ispc
            infile = [procpath,moor,'\microcat\',moor,'_',sprintf('%0.3d',i),'.microcat'];
        end
        
        
        % check if file exists
        fileopen=fopen(infile,'r');
        
        if fileopen>0
            % read data into vectors and then into structure array          
            [yy,mm,dd,hh,t,c,p,instrdpth] = ...
                rodbload(infile,'yy:mm:dd:hh:t:c:p:InstrDepth');
            jd=julian(yy,mm,dd,hh);
            
            % set to nan data outside time range:
            badtime = (jd<jd_start | jd>jd_end);
            p(badtime)=-9999;
            t(badtime)=-9999;
            c(badtime)=-9999;
            
            bad_data=find(p==-9999); p(bad_data)=nan;
            bad_data=find(t==-9999); t(bad_data)=nan;
            bad_data=find(c==-9999); c(bad_data)=nan;
            
            % interpolation of all the instrument on the same time_axis (the first instrument)
            if i==1
                timeref =  datenum(yy,mm,dd,hh,0*hh,0*hh);
                mcatdata.time(j,:) = datenum(yy,mm,dd,hh,0*hh,0*hh); %julian(YY,MM,DD,HH);
                goodtimeindex = j; % SJ: rarely needed except first instrument missing
                
                p(find(isnan(p))) = 99999;
                t(find(isnan(t))) = 99999;
                c(find(isnan(c))) = 99999;
                
                mcatdata.p(j,:) = p;
                mcatdata.t(j,:) = t;
                mcatdata.c(j,:) = c;
                
            else
                if exist('needtime'); % SJ: in the unlikely case that the first instrument is missing so we still need a base timeframe
                    timeref =  datenum(yy,mm,dd,hh,0*hh,0*hh);
                    %mcatdata.time(j,:) = datenum(yy,mm,dd,hh,0*hh,0*hh); %julian(YY,MM,DD,HH);
                    goodtimeindex = j; % SJ: rarely needed except first instrument missing
                    clear needtime
                    p(find(isnan(p))) = 99999;
                    t(find(isnan(t))) = 99999;
                    c(find(isnan(c))) = 99999;
                    
                    mcatdata.p(j,:) = p;
                    mcatdata.t(j,:) = t;
                    mcatdata.c(j,:) = c;
                    
                end
                
                timeinst =  datenum(yy,mm,dd,hh,0*hh,0*hh);
                mcatdata.time(j,:) = timeref; %julian(YY,MM,DD,HH);
                
                if length(p(~isnan(p)))>2
                    pint = interp1(timeinst(~isnan(p)),p(~isnan(p)),timeref);
                else
                    pint = timeref*0 + 99999;
                end
                if length(t(~isnan(t)))>2
                    tint = interp1(timeinst(~isnan(t)),t(~isnan(t)),timeref);
                else
                    tint = timeref*0 + 99999;
                end
                if length(c(~isnan(c)))>2
                    cint = interp1(timeinst(~isnan(c)),c(~isnan(c)),timeref);
                else
                    cint = timeref*0 + 99999;
                end
                
                baddatamaskp = interp1(timeinst,single(isnan(p)),timeref,'nearest');
                baddatamaskt = interp1(timeinst,single(isnan(t)),timeref,'nearest');
                baddatamaskc = interp1(timeinst,single(isnan(c)),timeref,'nearest');
                baddatamaskp(isnan(baddatamaskp)) = 1;
                baddatamaskt(isnan(baddatamaskt)) = 1;
                baddatamaskc(isnan(baddatamaskc)) = 1;
                
                pint(logical(baddatamaskp))= nan;
                tint(logical(baddatamaskt))= nan;
                cint(logical(baddatamaskc))= nan;
                
                pint(find(isnan(pint))) = 99999;
                tint(find(isnan(tint))) = 99999;
                cint(find(isnan(cint))) = 99999;
                
                mcatdata.p(j,:)    = pint;
                mcatdata.t(j,:)    = tint;
                mcatdata.c(j,:)    = cint;
            end
            mcatdata.serialnum(j) = serialno;
            mcatdata.instrdpth(j) = instrdpth;
        else
            
            ibad = [ibad i];
            jbad = [jbad j];
            needtime = 1;  % if the first instrument is missing, we still need
            % a timeref to interpolate the other instruments onto. take the next avbailable
        end
        
        j=j+1;
    end
    
    
    % for instrument without data
    for i=1:length(ibad);
        serialno = vecMC(ibad(i));
        mcatdata.serialnum(jbad(i)) = serialno;
        mcatdata.instrdpth(jbad(i)) = zMC(ibad(i));
        mcatdata.p(jbad(i),:)    = 99999*ones(size(mcatdata.p(jbad(i),:)));
        mcatdata.t(jbad(i),:)    = 99999*ones(size(mcatdata.p(jbad(i),:)));
        mcatdata.c(jbad(i),:)    = 99999*ones(size(mcatdata.p(jbad(i),:)));
    end
    
    % export oceansite format
    mcatinfo = moorinfo.mcat;
    mcatinfo.serial_num = mcatdata.serialnum;
    instrdepth = mcatdata.instrdpth; % nominal depth of the bin
    mcatdata.s = sw_salt( mcatdata.c/ sw_c3515,  mcatdata.t,  mcatdata.p);
    mcatdata.s(mcatdata.c==99999) = 99999;
    
    %cd export_Oceansites % SJ
    write_MCTD_to_NetCDF(netcdffilepath, moor, lat, lon, mcatinfo, instrdepth, mcatdata.time(goodtimeindex,:),  mcatdata.p, mcatdata.s, mcatdata.t, mcatdata.c)
    
    info=ncinfo([netcdffilepath '.nc']);
    ncdisp([netcdffilepath '.nc']);
    for ikk = 1:length(info.Variables)
        mcatncdata.(info.Variables(ikk).Name) = ncread([netcdffilepath '.nc'],info.Variables(ikk).Name);
    end
    
    rodb_to_oceansitesnetcdf_testplotncmcat
    
end

%--------------------------------------
% Now read in Aquadopp data if required
%--------------------------------------
if iiNOR>0
mintime = jd_start;
maxtime = jd_end;
%---------------------------------------------------------
% find longest non-nan timelimits   
    j=1;
    
    ibad = [];
    jbad = [];
    % loop to read one file at a time
        netcdffilepath = [ncdir '/' moor '_Nortek'];
    for i=1:length(vecNOR);
       serialno = vecNOR(i);
       disp('*************************************************************')
       disp(['Reading AQUADOPP - ',num2str(serialno)])
       disp('*************************************************************')


	       if isunix
        	   infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
       		elseif ispc
        	   infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
       		end	
    
       
       % check if file exists
       fileopen=fopen(infile,'r');
      
       if fileopen>0
           % read data into vectors and then into structure array

           [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd,instrdpth] = ...
               rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd:InstrDepth');
           jd=julian(yy,mm,dd,hh);           
           u = u/100;
           v = v/100;
           w = w/100;
           
           timeok = (jd>jd_start & jd<jd_end);
 
           bad_data=find(t==-9999); t(bad_data)=nan;
           bad_data=find(p==-9999); p(bad_data)=nan;
           bad_data=find(u==-9999/100); u(bad_data)=nan;
           bad_data=find(v==-9999/100); v(bad_data)=nan;
           bad_data=find(w==-9999/100); w(bad_data)=nan;

           if min(timeok) < mintime
	   	mintime = min(jd(timeok));
	   end
	   if  max(timeok) < maxtime
	   	maxtime = max(jd(timeok));
	   end
	     
       else
           
          ibad = [ibad i];
          jbad = [jbad j];
       end
       
       j=j+1;
    end
    

%---------------------------------------------------------
% Load data
    j=1;
    
    ibad = [];
    jbad = [];
    % loop to read one file at a time
        netcdffilepath = [ncdir '/' moor '_Nortek'];
    for i=1:length(vecNOR);
       serialno = vecNOR(i);
       disp('*************************************************************')
       disp(['Reading AQUADOPP - ',num2str(serialno)])
       disp('*************************************************************')


	       if isunix
        	   infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
       		elseif ispc
        	   infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
       		end	
    
       
       % check if file exists
       fileopen=fopen(infile,'r');
      
       if fileopen>0
           % read data into vectors and then into structure array

           [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd,instrdpth] = ...
               rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd:InstrDepth');
           jd=julian(yy,mm,dd,hh);           
           u = u/100;
           v = v/100;
           w = w/100;
           
           timeok = (jd>jd_start & jd<jd_end);
 
           bad_data=find(t==-9999); t(bad_data)=nan;
           bad_data=find(p==-9999); p(bad_data)=nan;
           bad_data=find(u==-9999/100); u(bad_data)=nan;
           bad_data=find(v==-9999/100); v(bad_data)=nan;
           bad_data=find(w==-9999/100); w(bad_data)=nan;

      % % % interpolation of all the instrument on the same time_axis (the first instrument) 
      % %    timeref =  datenum(yy(timeok),mm(timeok),dd(timeok),hh(timeok),0*hh(timeok),0*hh(timeok));            


        % interpolation on the longer instrument time
            if j==1
                timerefjd =  mintime:nanmean(diff(jd)):maxtime;
                timeref = datenum(gregorian(timerefjd));
            end
            
            timeinst =  datenum(yy,mm,dd,hh,0*hh,0*hh);
            nortekdata.time(j,:) = timeref; %julian(YY,MM,DD,HH);
            
            pint = interp1(timeinst(~isnan(p)),p(~isnan(p)),timeref); 
            tint = interp1(timeinst(~isnan(t)),t(~isnan(t)),timeref); 
            uint = interp1(timeinst(~isnan(u)),u(~isnan(u)),timeref); 
            vint = interp1(timeinst(~isnan(v)),v(~isnan(v)),timeref); 
            wint = interp1(timeinst(~isnan(w)),w(~isnan(w)),timeref); 
             
            baddatamaskp = interp1(timeinst,single(isnan(p)),timeref,'nearest');
            baddatamaskt = interp1(timeinst,single(isnan(t)),timeref,'nearest');      
            baddatamasku = interp1(timeinst,single(isnan(u)),timeref,'nearest');
            baddatamaskv = interp1(timeinst,single(isnan(v)),timeref,'nearest');   
            baddatamaskw = interp1(timeinst,single(isnan(w)),timeref,'nearest');

            
            baddatamaskp(isnan(baddatamaskp)) = 1;
            baddatamaskt(isnan(baddatamaskt)) = 1;            
            baddatamasku(isnan(baddatamasku)) = 1;        
            baddatamaskv(isnan(baddatamaskv)) = 1;        
            baddatamaskw(isnan(baddatamaskw)) = 1;        
            
            pint(logical(baddatamaskp))= nan;
            tint(logical(baddatamaskt))= nan;          
            uint(logical(baddatamasku))= nan;
            vint(logical(baddatamaskv))= nan;  
            wint(logical(baddatamaskw))= nan;             
            
            pint(find(isnan(pint))) = 99999; 
            tint(find(isnan(tint))) = 99999;
            uint(find(isnan(uint))) = 99999;   
            vint(find(isnan(vint))) = 99999;
            wint(find(isnan(wint))) = 99999;             
                        
            nortekdata.p(j,:)    = pint;
            nortekdata.t(j,:)    = tint;
            nortekdata.u(j,:) = uint;
            nortekdata.v(j,:) = vint;       
            nortekdata.w(j,:) = wint;               

         
            nortekdata.serialnum(j) = serialno;
            nortekdata.instrdpth(j) = instrdpth;     
       else
           
          ibad = [ibad i];
          jbad = [jbad j];
       end
       
       j=j+1;
    end
    


    % for instrument without data
    for i=1:length(ibad);
       serialno = vecNOR(ibad(i));
            nortekdata.serialnum(jbad(i)) = serialno;
            nortekdata.instrdpth(jbad(i)) = zNOR(ibad(i));          
            nortekdata.p(jbad(i),:)    = 99999*ones(size(nortekdata.p(jbad(i),:)));
            nortekdata.t(jbad(i),:)    = 99999*ones(size(nortekdata.p(jbad(i),:)));         
            nortekdata.u(jbad(i),:)    = 99999*ones(size(nortekdata.p(jbad(i),:)));      
            nortekdata.v(jbad(i),:)    = 99999*ones(size(nortekdata.p(jbad(i),:)));    
            nortekdata.w(jbad(i),:)    = 99999*ones(size(nortekdata.p(jbad(i),:)));      
    end

   % export oceansite format
   nortekinfo = moorinfo.nortek;
   nortekinfo.serial_num = nortekdata.serialnum;
   instrdepth = nortekdata.instrdpth; % nominal depth of the bin

   write_CM_to_NetCDF(netcdffilepath, moor, lat, lon, nortekinfo, instrdepth, nortekdata.time(1,:), nortekdata.p,  nortekdata.u, nortekdata.v,  nortekdata.w)

   info=ncinfo([netcdffilepath '.nc']);
   ncdisp([netcdffilepath '.nc']);
   for ikk = 1:length(info.Variables)
    nortekncdata.(info.Variables(ikk).Name) = ncread([netcdffilepath '.nc'],info.Variables(ikk).Name);
   end

   rodb_to_oceansitesnetcdf_testplotncnortek
 
end


%--------------------------------------
% Now read in ADCP data if required
%--------------------------------------
if iiADCP>0
    j=1;
    
    % loop to read one file at a time

    for i=1:length(vecADCP);
        
        netcdffilepath = [ncdir '/' moor '_ADCP' num2str(vecADCP(i))];
        
       serialno = vecADCP(i);
       disp('*************************************************************')
       disp(['Reading ADCP - ',num2str(serialno)])
       disp('*************************************************************')

    
    % first determine how many bin files are to be processed
    % trying to do automatically
    num_bins=dir([procpath,moor,'/adp/',moor,'_',num2str(vecADCP(i)),'_bin*' '.edt'])
    num_bins=length(num_bins)
    num_bins

    for j=1:num_bins % loop for total number of bins
        
        columns = ['YY:MM:DD:HH:Z:T:U:V:W:HDG:PIT:ROL:CS:CD:BEAM1SS:BEAM2SS:BEAM3SS'...
            ':BEAM4SS:BEAM1COR:BEAM2COR:BEAM3COR:BEAM4COR:EV:BEAM1PGP:BEAM2PGP:BEAM3PGP:BEAM4PGP:InstrDepth'];
        
        indep  = z(i);
        
        if j<=9
            infile  = [procpath,moor,'/adp/',moor,'_',num2str(vecADCP(i)),'_bin0',num2str(j),'.edt'];
        else
            infile  = [procpath,moor,'/adp/',moor,'_',num2str(vecADCP(i)),'_bin',num2str(j),'.edt'];
        end
                
        if exist(infile,'file')==0
            
            disp(['infile: ' infile ' does not exist.'])

        elseif exist(infile,'file')   > 0 

            [YY,MM,DD,HH,z,t,u,v,w,hdg,pit,rol,spd,direction,Amp1,Amp2,Amp3,Amp4,...
                Beam1Cor,Beam2Cor,Beam3Cor,Beam4Cor,err,PG1,PG2,PG3,PG4,instrument_depth] = ...
                rodbload(infile,[columns]);
           jd=julian(YY,MM,DD,HH);    
           
           u = u/100;
           v = v/100;
           w = w/100;
           err = err/100;
           
           timeok = (jd>jd_start & jd<jd_end);
           
           bad_data=find(t==-9999); t(bad_data)=99999;
           bad_data=find(z==-9999); z(bad_data)=99999;
           bad_data=find(u==-9999/100); u(bad_data)=99999;
           bad_data=find(v==-9999/100); v(bad_data)=99999;
           bad_data=find(w==-9999/100); w(bad_data)=99999;
           bad_data=find(hdg==-9999); hdg(bad_data)=99999;
           bad_data=find(pit==-9999); pit(bad_data)=99999;
           bad_data=find(rol==-9999); rol(bad_data)=99999;          
           bad_data=find(err==-9999/100); err(bad_data)=99999;       
           bad_data=find(PG1==-9999); PG1(bad_data)=99999;        
           bad_data=find(PG2==-9999); PG2(bad_data)=99999;   
           bad_data=find(PG3==-9999); PG3(bad_data)=99999;        
           bad_data=find(PG4==-9999); PG4(bad_data)=99999;   
           
            ADCPdata(i).time(j,:) = datenum(YY(timeok),MM(timeok),DD(timeok),HH(timeok),0*HH(timeok),0*HH(timeok)); %julian(YY,MM,DD,HH);
            ADCPdata(i).z(j,:) = z(timeok);
            ADCPdata(i).t(j,:) = t(timeok);          
            ADCPdata(i).u(j,:) = u(timeok);
            ADCPdata(i).v(j,:) = v(timeok);       
            ADCPdata(i).w(j,:) = w(timeok);   
            ADCPdata(i).err(j,:) = err(timeok);    
            ADCPdata(i).hdg(j,:) = hdg(timeok);       
            ADCPdata(i).pit(j,:) = pit(timeok);   
            ADCPdata(i).rol(j,:) = rol(timeok);                
            ADCPdata(i).PG1(j,:) = PG1(timeok);       
            ADCPdata(i).PG2(j,:) = PG2(timeok);    
            ADCPdata(i).PG3(j,:) = PG3(timeok);       
            ADCPdata(i).PG4(j,:) = PG4(timeok);       

        end
        
        
    end

        if exist('ADCPdata','var')~=1; continue;end
        ADCPdata(i).name = ['ADCP_' num2str(serialno)];
        ADCPdata(i).moor = moor;            


       j=j+1;

       % export oceansite format
       adcpinfo = moorinfo.adcp;
       adcpinfo.serial_num = serialno;
       bin_depth = round(nanmean(ADCPdata(i).z,2)); % nominal depth of the bin

       pres = sw_pres(ADCPdata(i).z,lat);

       write_ADCP_to_NetCDF(netcdffilepath, moor, lat, lon, adcpinfo, bin_depth, instrument_depth, pres, ADCPdata(i).time(1,:), ADCPdata(i).u, ADCPdata(i).v,  ADCPdata(i).w, ADCPdata(i).err)

       info=ncinfo([netcdffilepath '.nc']);  
       ncdisp([netcdffilepath '.nc']);
       for ikk = 1:length(info.Variables)
        ADCPncdata.(info.Variables(ikk).Name) = ncread([netcdffilepath '.nc'],info.Variables(ikk).Name);
       end
       
       rodb_to_oceansitesnetcdf_testplotncadcp
    end
end

% put the Oceansites conversion function here

