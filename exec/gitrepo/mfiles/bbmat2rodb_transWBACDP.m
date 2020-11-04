% PROGRAM - bbmat2rodb.m
%
% This code reads the bin mapped data from the .mat file and
% writes RODB format files for each bin
% .cor files are corrected for magnetic variation and speed of sound
% .edt iles have been de-spiked, filtered and interpolated onto 12hr grid
%
% ----------------------------------------------------------
% MODIFICATIONS FOR RB0901

% 1  -  lat is multiplied by ones(size(variable)) for functionality
% 2  -  header_file manually entered as the function to read it does not work
% 3  -  gregorian routine name changed to gregorian_decimal_hr to avoid
%       clashing with other gregorian functions - FEEL FREE TO CHANGE!!!!
% 4  -  file names code changed
% Suggestion: REWRITE THIS DAMN CODE TO USE rodbload INSTEAD!!!  
% (in order to enter the headers....)
% ----------------------------------------------------------
% Molina 200x  (modified 24 April 2009 Paul Wright)

clear all;
close all;

% PATHS FOR CRUISE RB0901
% -------------------------

addpath /Users/hydrosea3/Documents/RB0901/rapid/data/exec/rb0901/ADCP_processing/;
addpath /Users/hydrosea3/Documents/RB0901/rapid/data/moor/raw/rb0901/adp/;
addpath(genpath('/Users/hydrosea3/Documents/RB0901/rapid/data/exec/mfiles/sea/'));
inpath = ('/Users/hydrosea3/Documents/RB0901/rapid/data/moor/proc/wbadcp_5_200805/');
outpath = ('/Users/hydrosea3/Documents/RB0901/rapid/data/moor/proc/wbadcp_5_200805/');

%%%%%% load data here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  YOU WILL NEED TO CHECK INPUTS UP TO LINE 90 

% Set up mooring info here..
moor ='wbadcp_5_200805';


load wbadcp_5_200805_binmap.mat   % full record, use AA and BB to cut unwanted ensembles

% clear RDI* An* Ser*

%moor='36604001.edt'; % just use the Sontek CM data to grab lat and magdec info
% [lat,magdec]= jm_rodbload(moor,'Latitude:Mag_Variation'); % lat for sw_dpth or sw_pres routine

lat= 26.32;  magdec= -7.97; %  ALWAYS CHECK THIS!!!!!!!!



svel=cn;
% when you change nbins also change zd, uprof and vprof below and also znew
% (for akima)
nbins=40; % use up to bin ? for good data: **** USED ALL AVAILABLE BINS ******
% need to add 1.5 bins below the shadow zone (see Rainer's report)
%>> mean(zd(31,102:end-4))-mean(zs(102:end-4))
%ans =  39.9555  not within 24m (or 1.5 *  binsize)



% -------------------------
    % call function wbadcp_header - does not work :(
    % header_file=wbadcp_header(confile);%
    
    header_file = [ 'Mooring              = wbadcp_5_200805';
                    'Cruise               = RB0901         ';
                    'Instrument           = RDI 75kHz ADCP ';
                    'SerialNumber         = 1767           ';      
                    'Latitude             = 26 31.52 N     ';
                    'Longitude            = 76 52.12 W     ';
                    'WaterDepth           = 598            ';
                    'MagDeviation         = -7.97          ';
                    'StartDate            = 2008/04/23     ';
                    'StartTime            = 23:00          ';
                    'EndDate              = 2009/04/18     ';
                    'EndTime              = 14:00          '];  
% this is a bodge job :( (PW)% ----------------------------


AA=102;
BB=length(SerEnsembles)-4;
clear yy mm dd hh mn ss
jday=jday(AA:BB);
t=bb_t(AA:BB);         % Temperature


%%


confile= [inpath,moor,'.conf'];
ofile_cor=[outpath,moor,'.cor']; % this is hourly non filtered but corrected for magdec and sound speed
ofile_edt=[outpath,moor,'.edt']; % this is 12h res filtered
ofile_vel=[outpath,moor,'.vel'];


%%%%%% for acm_edit routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dt=(jday(2)-jday(1))*24    % sampling interval in hrs
nt=round(40/dt);            % 40 hr lowpass
filfac=10;                  % filter factor for despiking

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clf
set(gcf,'PaperPosition',[0 0 11 8.5])
set(gcf,'color',[1 1 1])

orient tall;


for i=1:nbins
    ofile_bin=[outpath,moor,'_bin',num2str(i),'.cor'] % this is hourly non filtered but corrected for magdec and sound speed
    ofile_12h=[outpath,moor,'_bin',num2str(i),'.edt'] % this is 12h res filtered
    uraw=bb_u(i,AA:BB);
    vraw=bb_v(i,AA:BB);
    p=sw_pres(zd(i,AA:BB),lat*ones(size(t))); % p=p'; edited: 'ones(size(t))'

    %         % correct for velocities here...
    %         % QUESTION: When correcting the velocities for sound speed, Do I assume that the
    %         % temperature and salinity at the transducer is constant above the
    %         % transducers? Or can I use a CTD cast to give me a sound speed profile
    %
    %
    %         % CORRECT FOR SOUND VELOCITY OF SEAWATER
    %         % FIRST CALC Sound Velocity as fxn of T P S


    % You use the sound speed measured at the transducers to
    % correct the velocities for all the bins, but you use the mean sound speed
    % over the profiling range to map the bins into depth.

    % FOR US MOOR (OLD RDI BB150) use sound speed at the transducer with S=35 (see RDI Practical Primer on ADCPs)
    % use sw_svel.m by Phil Morgan 1993 (CSIRO) for new sound speed
    %     if i==1
    %
    %      C=1449.2 + (4.6.*t) - (0.055.*t.*t) + (0.00029.*t.*t.*t) + ((0.0134.*t).*(35-35)) + (0.016.*(zd(1,AA:BB)));
    %     svel_old = mean(C); % was 35 used here as salinity? did they use 1500 for sound vel? see deployment records of ADCP
    %         S(1:length(t),1)=35; S=S';% estimate new sound speed at 35
    %         svel = sw_svel(S,t,p);
    %     end
    %
    %     S(1:length(t))=36;  % they used a fixed S=36 and T=14 for year 3
    %     svel = sw_svel(S,t,p);
    %     svel_old=sw_svel(36,14,597.5);
    %
    %     corfac= (svel/svel_old);
    %     us = uraw .* corfac;
    %     vs = vraw .* corfac;

    % UK YEAR 4: THEY USED measured T AND FIXED S=36
    us=uraw; % used S=36 and measured Tempreture therefore measured Sound speed
    vs=vraw;

    % CORRECT FOR MAGNETIC DEVIATION GIVEN BY MAGDEC
    % use uvrot by Visbeck
    % [u,v]=uvrot_visbeck(us,vs,magdec);
    % find magnetic declination (get median of deployment time)

    [ur,vr]=uvrot_visbeck(us,vs,magdec);
    vrotd(i,:)=vr; % for debugging

    %   smooth those that have erroneous surface vel data (< -1000 cm/s) by
    %   interpolation (interp1)
    jdate=jday;
    v_spike=find(vr < -1000);
    u_spike=find(ur < -1000);
    v_whatis=vr(v_spike);
    u_whatis=ur(u_spike);
    vr(v_spike)=NaN;
    ur(u_spike)=NaN;
    real_idx=isnan(vr)==0;
    jtemp=jdate(real_idx);
    v_temp=vr(real_idx);
    u_temp=ur(real_idx);
    %     [b,k,m]=unique(jtemp); % Remove duplicates from x
    %     v_new=interp1(b,v_temp(k),jdate);
    %     u_new=interp1(b,u_temp(k),jdate);
    v_new=interp1(jtemp,v_temp,jdate);
    u_new=interp1(jtemp,u_temp,jdate);
    clear b k m
    u=u_new;
    v=v_new;
    %     u=ur;
    %     v=vr;
    z=sw_dpth(p,lat*ones(size(t))); % edited 24Apr09 PW to make lat a matrix = p
    vsmth(i,:)=v;  % for debugging

    % now convert non filtered data into rodb
    [oyy,omm,odd,ohh]=gregorian_decimal_hr(jday);
    if i==1 % temp data for first bin only (closest to transducer)
        otab = [oyy;omm;odd;ohh;t;z;u;v];
        fmt = '%4d %2d %2d %8.1f %7.2f %7.2f %7.2f %7.2f\n';

    else
        otab = [oyy;omm;odd;ohh;z;u;v];
        fmt = '%4d %2d %2d %8.1f %7.2f %7.2f %7.2f\n';
    end


   fido = fopen(ofile_bin,'w');
   nn = 1;
   for nn = 1:12
    
      fprintf(fido,'%s',header_file(nn,1:38));
      fprintf(fido,'\n')
      nn = nn+1
 
   end
   

    binnum=['Bin                  = ',num2str(i)];
    mndp=num2str(mean(z));
    meandepth=['Depth                = ',mndp]; %mean depth
    fprintf(fido,binnum);
    fprintf(fido,'\n')
    fprintf(fido,meandepth);
    fprintf(fido,'\n')
    if i==1
        fprintf(fido,'Columns              = yy:mm:dd:hh:t:z:u:v\n');
    else
        fprintf(fido,'Columns              = yy:mm:dd:hh:z:u:v\n');
    end
    fprintf(fido,fmt,otab);
    fclose(fido);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Now filter data
    close all
    % this routine is from jm_acm_edit.m used for the CMs
    % edit data for dummies before filtering
    iu=isnan(u);
    u(iu)=meanmiss(u);
    iv=isnan(v);
    v(iv)=meanmiss(v);
    ip=isnan(p);
    p(ip)=meanmiss(p);
    if i==1
        it=isnan(t);
        t(it)=meanmiss(t);
    end

    % high pass data
    uhf=sqrt(mfilter(u,1,0,1/nt).^2);
    vhf=sqrt(mfilter(v,1,0,1/nt).^2);
    ii=find(uhf>filfac*mean(uhf) | vhf>filfac*mean(vhf));
    %             figure(1)        % spikes identified
    %             plot(jday,uhf,jday(ii),uhf(ii),'*r'); gregtick;
    %             figure(2)        % spikes identified
    %             plot(jday,vhf,jday(ii),vhf(ii),'or'); gregtick;


    % determine likely value at spike position first iteration
    uf=mfilter(u,1,1/nt,0);
    vf=mfilter(v,1,1/nt,0);
    pf=mfilter(p,1,1/nt,0);
    if i==1
        tf=mfilter(t,1,1/nt,0);
    end

    % put value
    u(ii)=uf(ii); v(ii)=vf(ii); p(ii)=pf(ii);
    if i==1
        t(ii)=tf(ii);
    end

    %despiked

    % low pass product again
    uf=mfilter(u,1,1/nt,0);
    vf=mfilter(v,1,1/nt,0);
    pf=mfilter(p,1,1/nt,0);
    df=sw_dpth(pf,lat*ones(size(t))); % calc for depth
    if i==1
        tf=mfilter(t,1,1/nt,0);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% make 12h resolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     tim=[fix(jday(1)+0.5):0.5:fix(jday(end))];
    tim=[fix(jday(1)+0.5):0.5:fix(jday(end))];

    %     [b,k,m]=unique(jday); % Remove duplicates from x
    %     x55=interp1(b, y(i), 5.5) % Obtain the interpolated value at x=5.5 - this succeeds
    %    ufhh=interp1(b,uf(k),tim);
    %    vfhh=interp1(b,vf(k),tim);
    %    pfhh=interp1(b,pf(k),tim);
    ufhh=interp1(jday,uf,tim);
    vfhh=interp1(jday,vf,tim);
    pfhh=interp1(jday,pf,tim);
    zfhh=sw_dpth(pfhh,lat*ones(size(pfhh))); %
    if i==1
        %     tfhh=interp1(b,tf(k),tim);
        tfhh=interp1(jday,tf,tim);
    end


    figure(i)
    clf
    plot(tim,ufhh); hold on;
    ylabel('velocity')
    plot(tim,vfhh,'-r');
    legend('u','v')
    timeaxis([2008,1,1,0])%gregtick;

    ufhhdb(i,:)=ufhh;
    vfhhdb(i,:)=vfhh;
    pfhhdb(i,:)=pfhh;
    zfhhdb(i,:)=sw_dpth(pfhhdb(i,:),lat*ones(size(pfhhdb(i,:)))); % calc for depth again but from filtered press data

    [yy,mm,dd,hh]=gregorian_decimal_hr(tim);

    % now convert FILTERED 12H RES data into rodb
    if i==1 % temp data for first bin only (closest to transducer)
        tab = [yy;mm;dd;hh;tfhh;zfhh;ufhh;vfhh];
        fmt = '%4d %2d %2d %8.1f %7.2f %7.2f %7.2f %7.2f\n';
    else
        tab = [yy;mm;dd;hh;zfhh;ufhh;vfhh];
        fmt = '%4d %2d %2d %8.1f %7.2f %7.2f %7.2f\n';
    end

    % call function wbadcp_header

   % header_file=wbadcp_header(confile);
    

    fido = fopen(ofile_12h,'w');
    
   nn = 1;
   for nn = 1:12
    
      fprintf(fido,'%s',header_file(nn,1:38));
      fprintf(fido,'\n')
      nn = nn+1
 
   end
          
    binnum=['Bin                  = ',num2str(i)];
    mndp=num2str(mean(zfhh(2:end)));
    meandepth=['Depth                = ',mndp]; %mean depth
    fprintf(fido,binnum);
    fprintf(fido,'\n')
    fprintf(fido,meandepth);
    fprintf(fido,'\n')

    if i==1
        fprintf(fido,'Columns              = yy:mm:dd:hh:t:z:u:v\n');
    else
        fprintf(fido,'Columns              = yy:mm:dd:hh:z:u:v\n');
    end
    fprintf(fido,fmt,tab);
    fclose(fido);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
return

%%%%%%%%%%%%%%%%% now working with 12h resolution data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lv=length(vfhh);


% stack up bins vertically (will need to flip zd and uprof and vprof later
% for akima routine, which doesnt take descending values
clear zd % clear old zd
zd=dfhhdb;
uprof=ufhhdb;
vprof=vfhhdb;

rdi_zd=zd;
rdi_uprof=uprof;
rdi_vprof=vprof;
rdi_yy=yy;
rdi_mm=mm;
rdi_dd=dd;
rdi_hh=hh;

%save 366_cor_deps_vels.mat rdi_zd rdi_uprof rdi_vprof rdi_yy rdi_mm rdi_dd rdi_hh
return

