% get rodb_file with raw data format als variables
% added pressure into rodb header -- JM 2/05
% this script also uses jm_acm_edit_header function which reinserts header
% in the data file at the end of this script (no need to manually insert
% *.conf file which contains header file)
% removed w --JM 12/2005
cd 'C:\Documents and Settings\Jonathan\Desktop\UKMoorings_2007_2008\UK_2007_2008\wb2_5_200702\proc'
close all;
clear all;

ifile={'wb2_5_200702_426.use',... % 100m, C failure, use microcat 5239 for C data
    'wb2_5_200702_428.use',... % 175m, all fine, current meter at 400m (sn 438) no data--FLOODED)'wb2_5_200702_5239.use'} % 175m  microcat
    'wb2_5_200702_443.use',... % 800m, all fine
    'wb2_5_200702_444.use',... %1200m, all fine
    'wb2_5_200702_507.use'}; %2050m, all fine

ofile={'wb2_5_200702_426.edt',... % 100m, C failure, use microcat 5239 for C data
    'wb2_5_200702_428.edt',... % 175m, all fine, current meter at 400m (sn 438) no data--FLOODED)'wb2_5_200702_5239.edt'} % 175m  microcat
    'wb2_5_200702_443.edt',... % 800m, all fine
    'wb2_5_200702_444.edt',... %1200m, all fine
    'wb2_5_200702_507.edt'}; %2050m, all fine

confile={'wb2_5_200702_426.conf',... % 100m, C failure, use microcat 5239 for C data
    'wb2_5_200702_428.conf',... % 175m, all fine, current meter at 400m (sn 438) no data--FLOODED)'wb2_5_200702_5239.conf'} % 175m  microcat
    'wb2_5_200702_443.conf',... % 800m, all fine
    'wb2_5_200702_444.conf',... %1200m, all fine
    'wb2_5_200702_507.conf'}; %2050m, all fine

filfac=10;      % --- filter factor for despiking

for mr=1:length(ifile);
    fidi=[ifile{mr}];

    vars='Mag_Variation:YY:MM:DD:HH:P:T:U:V';
    [magdec,yy,mm,dd,h,p,t,uraw,vraw]=jm_rodbload(fidi,vars);

    % make time to julian;
    time=julian(yy,mm,dd,h);
    dt=(time(2)-time(1))*24    % sampling interval
    nt=round(40/dt);            % 40 hr lowpass

    % CORRECT FOR SOUND VELOCITY OF SEAWATER
    % FIRST CALC Sound velocity as fxn of T P S
 
    S(1:length(t),1)=35; 
    svel = sw_svel(S,t,p);
    svel_old=1500;
    corfac= (svel/svel_old);

    us = uraw .* corfac;
    vs = vraw .* corfac;

    % CORRECT FOR MAGNETIC DEVIATION GIVEN BY MAGDEC
    %use uvrot by Visbeck
    [u,v]=uvrot_visbeck(us,vs,magdec);

    % edit data for dummies before filtering
    % uncomment this for debug only
    % v(2)=NaN;% for debug
    % u(2)=NaN;% for debug
    % p(2)=NaN;% for debug
    % t(2)=NaN;% for debug

    iu=isnan(u);
    iv=isnan(v);
    ip=isnan(p);
    it=isnan(t);

    ius=sum(iu);
    ivs=sum(iv);
    ips=sum(ip);
    its=sum(it);

    if (ius >0 | ivs >0 |ips > 0 | its >0)
        display('Inspect NaNs in your records')
        return
        iu_idx=isnan(u)==0;
        iv_idx=isnan(v)==0;
        it_idx=isnan(t)==0;
        ip_idx=isnan(p)==0;
        u=u(iu_idx);
        v=v(iv_idx);
        t=t(it_idx);
        p=p(ip_idx);
        timef=time(iv_idx);
        ufix=interp1(timef,u,time); clear u
        u=ufix;
        vfix=interp1(timef,v,time); clear v
        v=vfix;
        tfix=interp1(timef,t,time); clear t
        t=tfix;
        pfix=interp1(timef,p,time); clear p
        p=pfix;
    end
    % u(iu)=meanmiss(u);
    % v(iv)=meanmiss(v);
    % p(ip)=meanmiss(p);
    % t(it)=meanmiss(t);

    % high pass data
    uhf=sqrt(mfilter(u,1,0,1/nt).^2);
    vhf=sqrt(mfilter(v,1,0,1/nt).^2);
    ii=find(uhf>filfac*mean(uhf) | vhf>filfac*mean(vhf));
    figure(1)
    plot(time,uhf,time(ii),uhf(ii),'*r'); gregtick;
    figure(2)
    plot(time,vhf,time(ii),vhf(ii),'or'); gregtick;

    % spikes identified

    % determine likely value at spike position first iteration
    uf=mfilter(u,1,1/nt,0);
    vf=mfilter(v,1,1/nt,0);
    pf=mfilter(p,1,1/nt,0);
    tf=mfilter(t,1,1/nt,0);

    % put value
    u(ii)=uf(ii); v(ii)=vf(ii); p(ii)=pf(ii); t(ii)=tf(ii);
    %depiked

    % low pass product again
    uf=mfilter(u,1,1/nt,0);
    vf=mfilter(v,1,1/nt,0);
    pf=mfilter(p,1,1/nt,0);
    tf=mfilter(t,1,1/nt,0);


    % 12h resolution
    %    tim=[fix(time(1)+0.5):0.5:fix(time(end))];
    tim=[fix(time(1)):0.5:fix(time(end))];
    vf1=interp1(time,vf,tim);
    uf1=interp1(time,uf,tim);
    pf1=interp1(time,pf,tim);
    tf1=interp1(time,tf,tim);

    figure(3)
    plot(tim,uf1); hold on;
    plot(tim,vf1,'-r'); gregtick;  hold off;


    % now put everything back in rodb

    [yy,mm,dd,hh]=gregorian(tim);

    tab = [yy;mm;dd;hh;tf1;pf1;uf1;vf1];

    % call function jm_acm_edit_header

    header_file=jm_acm_edit_header(confile{mr});

    fmt = '%4d %2d %2d %8.4f %7.3f %7.2f %7.2f %7.2f\n';
    fido = fopen(ofile{mr},'w');
    fprintf(fido,'%s',header_file{:});
    fprintf(fido,'Columns              = yy:mm:dd:hh:t:p:u:v\n');
    fprintf(fido,fmt,tab);
    fclose(fido);
    %print -dpsc test
    clear corfac d*  fidi fido fmt h* ii ip it iu iv magdec mm n*  p* q* r* s* S* t* u* v* w* x* y* z*
end
