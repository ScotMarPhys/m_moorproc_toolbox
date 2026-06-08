function jm_arg2rodb(filterSpec)
% ARG2RODB Convert Argonaut  raw data to RODB format.
% Input: CONF file (configuration file/header)
% Load: RAW file (stripped off data)
% Output: COR file (corrected for magnetic declination and sound speed) --JM 12/2005
% added pressure column --JM 2/2005
% allowed to load data that has no text in filename, eg. 123.dat --JM 3/2005
% changed column designations for temperature, pressure and w vel to corresspond
% to latest version of ViewArgonaut v3.20--JM 3/2005
% removed w vel from output--JM 12/2005
% added correction for sound speed--JM 12/2005

d = dir(filterSpec);
dirname = filterSpec(1:max(findstr(filterSpec,filesep)));

for k = 1:length(d)
    confile = fullfile(dirname,d(k).name);
    header = readconf(confile);
    var = rodbhead(str2mat(header{:}));
    i = find(var == rodbhead('Filename'));
    outfile = sscanf(header{i},'%*s = %s');

    i = find(var == rodbhead('Mag_Variation')); % Magnetic Declination to replace Magnetic Deviation (also changed in jm_rodbhead.m
    magdec = sscanf(header{i},'%*s = %f');

    i = max(findstr(d(k).name,'.'));
    rawfile = sprintf('%s.raw',d(k).name(1:i-1));
    rawdat = sprintf('%s',d(k).name(1:i-1));
    rawfile = fullfile(dirname,rawfile);

    % get the data by loading it
    eval(['load ' rawfile]);
    eval(['a=X' rawdat ';']); %use for filenames that dont have text in name eg, 123.dat
    %eval(['a=' rawdat ';']); %use for filenames that have text eg, M123.dat
    [m,n] = size(a);
    yy=a(:,1);mm=a(:,2);dd=a(:,3);
    hh=a(:,4)+a(:,5)/60+a(:,6)/3600;

    jdate=julian(yy,mm,dd,hh);

    % CORRECT FOR SOUND VELOCITY OF SEAWATER
    % FIRST CALC Sound Velocity as fxn of T P S
    T=a(:,29);
    P=a(:,30);
    S(1:length(T),1)=35; 
    svel = sw_svel(S,T,P);
    svel_old=1500;
    corfac= (svel/svel_old);
    uraw=a(:,7);
    vraw=a(:,8);
    us = uraw .* corfac;
    vs = vraw .* corfac;


    % CORRECT FOR MAGNETIC DEVIATION GIVEN BY MAGDEC
    %use uvrot by Visbeck
    [u,v]=uvrot_visbeck(us,vs,magdec);

    % be sure to have corresponding columns for variables (t,p,u,v) when using other versions of
    % ViewArgonaut, I used v3.20 --JM 3/2005
    tab = [yy,mm,dd,hh,a(:,29),a(:,30),u,v];

    fmt = '%4d %2d %2d %8.4f  %7.3f   %7.3f %7.2f %7.2f\n';
    fid = fopen(outfile,'w');
    fprintf(fid,'%s',header{:});
    fprintf(fid,'Columns              = yy:mm:dd:hh:t:p:u:v\n');
    fprintf(fid,fmt,tab');
    fclose(fid);

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Below for comparison of RZ's and Visbeck's uvrot script only and for data check
%     figure(1) % check T P S and svel
%     subplot(4,1,1);plot(jdate,T);jm_gregtick_noyr;hold on;
%     subplot(4,1,2);plot(jdate,P);jm_gregtick_noyr;
%     subplot(4,1,3);plot(jdate,S);jm_gregtick_noyr;
%     subplot(4,1,4);plot(jdate,svel);gregtick;
% 
%     figure(2) % look at velocity time series
%     [uo,vo]=uvrot(uraw,vraw,magdec); %old using RZ's uvrot
%     [ur,vr]=uvrot_visbeck(uraw,vraw,magdec); % rotated only using Visbeck's
%     subplot(2,1,1);plot(jdate,u,'r',jdate,uo,'b',jdate,us,'k',jdate,ur,'g');
%     ylabel('U (m/s)')
%     jm_gregtick_noyr
%     subplot(2,1,2);plot(jdate,v,'r',jdate,vo,'b',jdate,vs,'k',jdate,vr,'g');
%     ylabel('V (m/s)')
%     legend('Magnetic Declination and Sound Speed Corrected Velocity','RZ uvrot','Sound Corrected Only','Rotated Only')
%     jm_gregtick_noyr
%     
%     figure(3) % look at velocity sticks
%     vscale=0;
%     tax1 = julian(2003,9,1,0);
%     tax2 = julian(2005,4,1,0);
%     h=quiver(jdate(1:5),jdate(1:5).*0,ur(1:5),vr(1:5),vscale,'b'); hold on;
%     h2=quiver(jdate(1:5),jdate(1:5).*0,uo(1:5),vo(1:5),vscale,'k');
%     h3=quiver(jdate(1:5),jdate(1:5).*0,u(1:5),v(1:5),vscale,'r');
%     % set(h,'showarrowhead','off');
%     % set(h2,'showarrowhead','off');
%     % this plots a y=0 line on the current plot; same a yzero.m by D. Slater
%     % plot([tax1 tax2],[0 0],'k-')
%     box off
%     ylabel('Velocity (cm/s)','FontSize',7);
%     jm_gregtick_noyr
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%-------------------------------------------------------------------------------

function header = readconf(filename)
%READCONF
%  HEADER = READCONF(FILENAME)

fid = fopen(filename,'r');
i = 1;
while 1
    header{i} = fgets(fid);
    if header{i} == -1
        break
    end
    i = i + 1;
end
fclose(fid);
header = header(1:end-1);



