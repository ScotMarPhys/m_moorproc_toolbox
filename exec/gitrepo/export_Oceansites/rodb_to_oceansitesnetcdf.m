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
function rodb_to_oceansitesnetcdf(moor,moorinfo)

if nargin <1
    help rodb_to_oceansitesnetcdf
    return
end

global MOORPROC_G

pdo = moor_inoutpaths('oceansites',moor);
if ~exist(fileparts(pdo.ncpre),'dir')
    mkdir(fileparts(pdo.ncpre))
end

% Load vectors of mooring information
% id instrument id, sn serial number, z nominal depth of each instrument
% s_t, e_t, s_d, e_d start and end times and dates
% lat lon mooring position, wd corrected water depth (m)
% mr mooring name
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(pdo.infofile,...
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
    isgood = true(size(vecMC));

    ncfilep = [pdo.ncpre '_MCAT'];
    pd = moor_inoutpaths('microcat',moor);

    % loop to read one file at a time
    for j=1:length(vecMC)
        serialno = vecMC(j);
        instrdpth = zMC(j);
        mcatdata.serialnum(j) = serialno;
        mcatdata.instrdpth(j) = instrdpth;        
        disp('*************************************************************')
        disp(['Reading MICROCAT - ',num2str(serialno)])
        disp('*************************************************************')
        infile = fullfile(pd.stage3path,sprintf(pd.stage3form,j));
        
        % check if file exists
        fileopen=fopen(infile,'r');
        if fileopen>0
            % read data into vectors and then into structure array          
            [yy,mm,dd,hh,t,c,p,instrdpth] = ...
                rodbload(infile,'yy:mm:dd:hh:t:c:p:InstrDepth');
            jd=julian(yy,mm,dd,hh);
            % set to nan data outside time range:
            badtime = (jd<jd_start | jd>jd_end);
            p(badtime | p==-9999) = nan;
            t(badtime | t==-9999) = nan;
            c(badtime | c==-9999) = nan;
            
            % interpolation of all the instrument on the same time_axis (the first instrument)
            if ~exist('mctimeref','var')
                mctimeref =  datenum(yy,mm,dd,hh,0*hh,0*hh);
                mcatdata.time(j,:) = datenum(yy,mm,dd,hh,0*hh,0*hh); %julian(YY,MM,DD,HH);
                goodtimeindex = j; % SJ: rarely needed except first instrument missing
                mcatdata.p(j,:) = p;
                mcatdata.t(j,:) = t;
                mcatdata.c(j,:) = c;
            else
                timeinst =  datenum(yy,mm,dd,hh,0*hh,0*hh);
                mcatdata.time(j,:) = mctimeref; %julian(YY,MM,DD,HH);
                mcatdata.p(j,:) = baddatamask(timeinst,p,mctimeref,nan);
                mcatdata.t(j,:) = baddatamask(timeinst,t,mctimeref,nan);
                mcatdata.c(j,:) = baddatamask(timeinst,c,mctimeref,nan);
            end
        else
            isgood(j) = false;
        end

            ibad = ~isgood;            
        mcatdata.p(ibad,:)    = 99999;
        mcatdata.t(ibad,:)    = 99999;
        mcatdata.c(ibad,:)    = 99999;
    end
    
    % export oceansite format
    mcatinfo = moorinfo.mcat;
    mcatinfo.serial_num = mcatdata.serialnum;
    instrdepth = mcatdata.instrdpth; % nominal depth of the bin
    mcatdata.s = gsw_SP_from_R( mcatdata.c/ gsw_C3515,  mcatdata.t,  mcatdata.p);
    mcatdata.s(isnan(mcatdata.s))=99999;
    mcatdata.p(isnan(mcatdata.p))=99999;
    mcatdata.t(isnan(mcatdata.t))=99999;

    %cd export_Oceansites % SJ
    write_MCTD_to_NetCDF(ncfilep, moor, lat, lon, mcatinfo, instrdepth, mcatdata.time(goodtimeindex,:),  mcatdata.p, mcatdata.s, mcatdata.t, mcatdata.c)

    info=ncinfo([ncfilep '.nc']);
    ncdisp([ncfilep '.nc']);
    for ikk = 1:length(info.Variables)
        mcatncdata.(info.Variables(ikk).Name) = ncread([ncfilep '.nc'],info.Variables(ikk).Name);
    end

    rodb_to_oceansitesnetcdf_testplotncmcat

end

%--------------------------------------
% Now read in Aquadopp data if required
%--------------------------------------
if iiNOR>0
    mintime = jd_start;
    maxtime = jd_end;
    isgood = true(size(vecNOR));
    ncfilep = [pdo.ncpre '_Nortek'];
    pd = moor_inoutpaths('nortek',moor);

    %---------------------------------------------------------
    % find longest non-nan timelimits
    for j=1:length(vecNOR)
        serialno = vecNOR(j);
        instrdpth = zNOR(j);
        nortekdata.serialnum(j) = serialno;
        nortekdata.instrdpth(j) = instrdpth;
        disp('*************************************************************')
        disp(['Reading AQUADOPP - ',num2str(serialno)])
        disp('*************************************************************')
        infile = fullfile(pd.stage3path,sprintf(pd.stage3form,serialno));
        % check if file exists
        fileopen=fopen(infile,'r');
        if fileopen>0
            % read data into vectors and then into structure array
            [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd,instrdpth] = ...
                rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd:InstrDepth');
            jd=julian(yy,mm,dd,hh);
            timeok = (jd>jd_start & jd<jd_end);
            mintime = min(mintime,jd(timeok));
            maxtime = max(maxtime,jd(timeok));
        else
            isgood(j) = false;
        end
    end

    %---------------------------------------------------------
    % Load data
    for j=1:length(vecNOR)
        if isgood(j)
        serialno = vecNOR(j);
        disp('*************************************************************')
        disp(['Reading AQUADOPP - ',num2str(serialno)])
        disp('*************************************************************')
        infile = fullfile(pd.stage3path,sprintf(pd.stage3form,serialno));
        % check if file exists
        fileopen=fopen(infile,'r');
        % read data into vectors and then into structure array
        [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd,instrdpth] = ...
            rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd:InstrDepth');
        jd=julian(yy,mm,dd,hh);
        t(t==-9999) = nan;
        p(p==-9999) = nan;
        u(u==-9999) = nan;
        v(v==-9999) = nan;
        w(w==-9999) = nan;
        u = u/100;
        v = v/100;
        w = w/100;

        % % % interpolation of all the instrument on the same time_axis (the first instrument)
        % %    timeref =  datenum(yy(timeok),mm(timeok),dd(timeok),hh(timeok),0*hh(timeok),0*hh(timeok));
        % interpolation on the longer instrument time
        if ~exist('nortimeref','var')
            timerefjd =  mintime:nanmean(diff(jd)):maxtime;
            nortimeref = datenum(gregorian(timerefjd));
        end
        timeinst =  datenum(yy,mm,dd,hh,0*hh,0*hh);
        nortekdata.time(j,:) = nortimeref; %julian(YY,MM,DD,HH);
        nortekdata.p(j,:) = baddatamask(timeinst,p,nortimeref,99999);
        nortekdata.t(j,:) = baddatamask(timeinst,t,nortimeref,99999);
        nortekdata.u(j,:) = baddatamask(timeinst,u,nortimeref,99999);
        nortekdata.v(j,:) = baddatamask(timeinst,v,nortimeref,99999);
        nortekdata.w(j,:) = baddatamask(timeinst,w,nortimeref,99999);
        end
    end

    % export oceansite format
    nortekinfo = moorinfo.nortek;
    nortekinfo.serial_num = nortekdata.serialnum;
    instrdepth = nortekdata.instrdpth; % nominal depth of the bin

    write_CM_to_NetCDF(ncfilep, moor, lat, lon, nortekinfo, instrdepth, nortekdata.time(1,:), nortekdata.p,  nortekdata.u, nortekdata.v,  nortekdata.w)

    info=ncinfo([ncfilep '.nc']);
    ncdisp([ncfilep '.nc']);
    for ikk = 1:length(info.Variables)
        nortekncdata.(info.Variables(ikk).Name) = ncread([ncfilep '.nc'],info.Variables(ikk).Name);
    end

    rodb_to_oceansitesnetcdf_testplotncnortek

end


%--------------------------------------
% Now read in ADCP data if required
%--------------------------------------
if iiADCP>0
    isgood = true(size(vecADCP));
    pd = moor_inoutpaths('nortek',moor);

    for j=1:length(vecADCP)
        ncfilep = [pdo.ncpre '_ADCP' num2str(vecADCP(j))];
        serialno = vecADCP(j);
        disp('*************************************************************')
        disp(['Reading ADCP - ',num2str(serialno)])
        disp('*************************************************************')
        % first determine how many bin files are to be processed
    % trying to do automatically
    num_bins=dir(fullfile(pd.stage3path,sprintf([pd.stage3form(1:end-9) '*.edt'],j)));
    num_bins=length(num_bins);
    num_bins

    for jz=1:num_bins % loop for total number of bins  
        columns = ['YY:MM:DD:HH:Z:T:U:V:W:HDG:PIT:ROL:CS:CD:BEAM1SS:BEAM2SS:BEAM3SS'...
            ':BEAM4SS:BEAM1COR:BEAM2COR:BEAM3COR:BEAM4COR:EV:BEAM1PGP:BEAM2PGP:BEAM3PGP:BEAM4PGP:InstrDepth'];        
        indep  = z(jz);
        infile = fullfile(pd.stage3path,sprintf(pd.stage3form,j,jz));
        if exist(infile,'file')==0            
            disp(['infile: ' infile ' does not exist.'])
        elseif exist(infile,'file')   > 0 
            [YY,MM,DD,HH,z,t,u,v,w,hdg,pit,rol,spd,direction,Amp1,Amp2,Amp3,Amp4,...
                Beam1Cor,Beam2Cor,Beam3Cor,Beam4Cor,err,PG1,PG2,PG3,PG4,instrument_depth] = ...
                rodbload(infile,[columns]);
           jd=julian(YY,MM,DD,HH);    
           t(t==-9999) = 99999;
           z(z==-9999) = 99999;
           u(u==-9999) = 99999;
           v(v==-9999) = 99999;
           w(w==-9999) = 99999;
           hdg(hdg==-9999) = 99999;
           pit(pit==-9999) = 99999;
           rol(rol==-9999) = 99999;
           err(err==-9999) = 99999;
           PG1(PG1==-9999) = 99999;
           PG2(PG2==-9999) = 99999;
           PG3(PG3==-9999) = 99999;
           PG4(PG4==-9999) = 99999;
           u = u/100;
           v = v/100;
           w = w/100;
           err = err/100;
           
           timeok = (jd>jd_start & jd<jd_end);
           
            ADCPdata(j).time(jz,:) = datenum(YY(timeok),MM(timeok),DD(timeok),HH(timeok),0*HH(timeok),0*HH(timeok)); %julian(YY,MM,DD,HH);
            ADCPdata(j).z(jz,:) = z(timeok);
            ADCPdata(j).t(jz,:) = t(timeok);          
            ADCPdata(j).u(jz,:) = u(timeok);
            ADCPdata(j).v(jz,:) = v(timeok);       
            ADCPdata(j).w(jz,:) = w(timeok);   
            ADCPdata(j).err(jz,:) = err(timeok);    
            ADCPdata(j).hdg(jz,:) = hdg(timeok);       
            ADCPdata(j).pit(jz,:) = pit(timeok);   
            ADCPdata(j).rol(jz,:) = rol(timeok);                
            ADCPdata(j).PG1(jz,:) = PG1(timeok);       
            ADCPdata(j).PG2(jz,:) = PG2(timeok);    
            ADCPdata(j).PG3(jz,:) = PG3(timeok);       
            ADCPdata(j).PG4(jz,:) = PG4(timeok);       
        end
        
    end

        if exist('ADCPdata','var')~=1; continue;end
        ADCPdata(j).name = ['ADCP_' num2str(serialno)];
        ADCPdata(j).moor = moor;            

        % export oceansite format
        adcpinfo = moorinfo.adcp;
        adcpinfo.serial_num = serialno;
        bin_depth = round(nanmean(ADCPdata(j).z,2)); % nominal depth of the bin

        pres = sw_pres(ADCPdata(j).z,lat);
        write_ADCP_to_NetCDF(ncfilep, moor, lat, lon, adcpinfo, bin_depth, instrument_depth, pres, ADCPdata(j).time(1,:), ADCPdata(j).u, ADCPdata(j).v,  ADCPdata(j).w, ADCPdata(j).err)

        info=ncinfo([ncfilep '.nc']);
        ncdisp([ncfilep '.nc']);
        for ikk = 1:length(info.Variables)
            ADCPncdata.(info.Variables(ikk).Name) = ncread([ncfilep '.nc'],info.Variables(ikk).Name);
        end

        rodb_to_oceansitesnetcdf_testplotncadcp
    end
end

% put the Oceansites conversion function here

function datai = baddatamask(timeinst,data,timeref,badval)
if length(data(~isnan(data)))>2
    datai = interp1(timeinst(~isnan(data)),data(~isnan(data)),timeref);
else
    datai = NaN;
end
baddatamask = interp1(timeinst,single(isnan(data)),timeref,'nearest');
baddatamask(isnan(baddatamask)) = 1;
datai(logical(baddatamask))=nan;
if ~isnan(badval)
    datai(isnan(datai)) = badval;
end

