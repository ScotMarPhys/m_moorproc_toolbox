% function read_flag_raw_adcp(moor,'inpath','outpath','procpath')
% Function to read and flag RDI ADCP data .
% raw ADCP data in .mat format as exported by WINADCP
%
%
% required inputs: moor - mooring name e.g. 'wb2_3_200606'
%                  inpath - if not using standard paths - standarad paths
%                  not setup yet. so required input
% optional inputs: procpath - if not using standard paths for info.dat file
%                  outpath - if not using standard paths
%
% functions called: hms2h.m
%                   julian.m
%                   rodbload.m
%
%
% Houpert Loic, 7/10/16, adapt from Colin Griffiths script
% Houpert Loic 19/10/20, add panel with v component to be able to detect
% also no coherent bins
%

function read_flag_raw_adcp(moor, varargin)

if nargin==0
    help adcp2rodb_01
    return
end
display(' ')
display('------------------------------')
display('Starting stage 1 adcp2rodb_01')
display('------------------------------')
display(' ')

% check for optional arguments
a=strmatch('procpath',varargin,'exact');
if a>0
    procpath=char(varargin(a+1));
else
    procpath='/Users/hydrosea5/Desktop/RB1201/rapid/data/moor/proc';
end

a=strmatch('inpath',varargin,'exact');
if a>0
    inpath=char(varargin(a+1));
else
    inpath=['/Users/hydrosea5/Desktop/RB1201/rapid/data/moor/raw/' moor '/'];
end

a=strmatch('outpath',varargin,'exact');
if a>0
    outpath=char(varargin(a+1));
else
    outpath = './';
end


% --- get moring information from infofile 
%infofile =[procpath '/' moor '/' moor 'info.dat'];
infofile =[procpath moor '/' moor 'info.dat'];

if isunix
        [gash, operator]=system('whoami'); % This line will not work if run from a PC.
        % % gash is not used, but a second variable needs to be specified for the system command
else
    operator=input('Enter opeartor name: ','s');
end

% ----- read infofile / open logfile  ------------------------------------

infovar = 'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth'; 
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd]  =  rodbload(infofile,infovar);

fidlog   = fopen([outpath,moor,'_ADCP_stage1.log'],'a');
fprintf(fidlog,'Transformation of ADCP .mat data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));
fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n\n\n',lon);

e_d
s_d


bg = julian([s_d(:)' hms2h([s_t(:)' 0])]); %start
ed = julian([e_d(:)' hms2h([e_t(:)' 0])]); %end



vec=find((id>=319) & (id <=328)); % Possible ADCP codes - taken from IMP moorings package


% -------- load data --------------
for i = 1:length(vec)
    infile=[inpath num2str(sn(vec((i)))) '_data.mat']
    
    indep=z(vec(i));
    if length(indep)>1
        indep=indep(1);
    end
    fprintf(fidlog,'infile : %s\n',infile);
    fprintf(fidlog,'ADCP serial number  : %d\n',sn(vec(i)));

    plot_title = [moor ' ADCP s/n ' num2str(sn(vec(i)))];
    
    outfile=[moor '_' num2str(sn(vec((i)))) '_pg1pg4'];

    load(infile);
    deviation=0.0; % no correction applied
    maxs=80.0;   
    limit=64;top=100;ts=6;iplot=2;jplot=3; % Set plotting and O/P control variables (jplot >= iplot)
    averageamp=1;percentG=1;pg1plot=1;corplot=1;tsplot=1;amps=1;vertical=0;Emplot=0;Nmplot=0;
    averagecor=1;output=0;espec=1;nspec=1;scatter=1;histogram=0;pvd=1;ex=0;

    set(0,'defaultaxesfontname','TimesRoman');
    set(0,'defaultaxesfontsize',12);
    set(0,'defaulttextfontname','TimesRoman');
    set(0,'defaulttextfontsize',12);

    first=1;
    last = length(SerYear);
    good=last-first+1;
    mtime=datenum(SerYear(first:last),SerMon(first:last),SerDay(first:last),SerHour(first:last),SerMin(first:last),SerSec(first:last));

    bins  = length(SerBins);         % number of depth bins
    scans = length(AnBIT);  % number of ensembles in record
    scannber = 1:length(AnBIT(first:last));
    binnber  = SerBins;
    pressure=AnDepthmm/100.0;        % Pressure

    E=SerEmmpersec(first:last,:)/10.; % u relative to adcp
    N=SerNmmpersec(first:last,:)/10.; % v relative to adcp
    W=SerVmmpersec(first:last,:)/10.; % vertical velocity
    e=SerErmmpersec(first:last,:)/10.; % error velocity

    pg1=SerPG1(first:last,:)'; % percent good for various parameters - see figure titles
    pg2=SerPG2(first:last,:)';
    pg3=SerPG3(first:last,:)';
    pg4=SerPG4(first:last,:)';

    cor1=SerC1cnt(first:last,:)'; % beam correlations
    cor2=SerC2cnt(first:last,:)';
    cor3=SerC3cnt(first:last,:)';
    cor4=SerC4cnt(first:last,:)';
    cor=SerCAcnt(first:last,:);

    amp1=SerEA1cnt(first:last,:)'; % beam amplitudes
    amp2=SerEA2cnt(first:last,:)';
    amp3=SerEA3cnt(first:last,:)';
    amp4=SerEA4cnt(first:last,:)';
    amp=SerEAAcnt(first:last,:);

    % No qc control - relying on %good 4 to select number of good bins





    starttime=mtime(1);interval=(mtime(2)-mtime(1));

    mtime=mtime+interval/2.0; % Advance Times by half of recording interval

    % The RDI ADCP does not time centre each ensemble

    interval=interval*24.0;
    days=mtime-starttime;
% 
%     if pg1plot>0; 
% 
%     figure;
% 
%     maxpg1=max(max(pg1))
%     if maxpg1>0;
%      subplot(2,2,1);
%      pg1=pg1*100.0/maxpg1;
%      x=[1 scans];y=[1 bins];
%      color1=jet(100);color1(1,:)=0.0;
%      for colour=limit:top
%         color1(colour,:)=1.0;
%      end;
%      h=image(pg1,'xdata',x,'ydata',y);colormap(color1);colorbar;
%      text(100,1,['Max value  = ',num2str(maxpg1)]);
%      title(' % Good pg1 - 3 beam good');grid on;
%      xlabel('Scans');ylabel('bins');
%     end;
% 
%     maxpg2=max(max(pg2))
%     if maxpg2>0;
%      subplot(2,2,2);
%      pg2=pg2*100.0/maxpg2;
%      x=[1 scans];y=[1 bins];
%      color1=jet(100);color1(1,:)=0.0;
%      for colour=limit:top
%         color1(colour,:)=1.0;
%      end;
%      h=image(pg2,'xdata',x,'ydata',y);colormap(color1);colorbar;
%      text(100,1,['Max value  = ',num2str(maxpg2)]);
%      title(' % Good pg2 - error vel > threshold');grid on;
%      xlabel('Scans');ylabel('bins');
%     end;
% 
%     maxpg3=max(max(pg3))
%     if maxpg3>0
%      subplot(2,2,3);
%      pg3=pg3*100.0/maxpg3;
%      x=[1 scans];y=[1 bins];
%      color1=jet(100);color1(1,:)=0.0;
%      for colour=limit:top
%        color1(colour,:)=1.0;
%      end;
%      h=image(pg3,'xdata',x,'ydata',y);colormap(color1);colorbar;
%      text(100,1,['Max value  = ',num2str(maxpg3)]);
%      title(' % Good pg3 - 1+ beam bad - no vels');grid on;
%      xlabel('Scans');ylabel('bins');
%     end;
% 
%     maxpg4=max(max(pg4))
%     if maxpg4>0
%      subplot(2,2,4);
%      maxpg4=max(max(pg4));
%      pg4=pg4*100.0/maxpg4;
%      x=[1 scans];y=[1 bins];
%      color1=jet(100);color1(1,:)=0.0;
%      for colour=limit:top
%         color1(colour,:)=1.0;
%      end;
%      h=image(pg4,'xdata',x,'ydata',y);colormap(color1);colorbar;
%      text(100,1,['Max value  = ',num2str(maxpg4)]);
%      title(' % Good pg4 - 4 beam solution');grid on;
%      xlabel('Scans');ylabel('bins');
%     end
% 
% 
%     end;
% 
% 
% 
%     if percentG>0
%      figure;
%      maxpg4=max(max(pg4));
%      pg4=pg4*100.0/maxpg4;
%      x=[1 scans];y=[1 bins];
%      color1=jet(100);color1(1,:)=0.0;
%      for colour=limit:top
%         color1(colour,:)=1.0;
%      end;
%      h=image(pg4,'xdata',x,'ydata',y);colormap(color1);colorbar;
%      text(100,1,['Max value  = ',num2str(maxpg4)]);
%      title(' % Good pg4 - 4 beam solution');grid on;
%      xlabel('Scans');ylabel('bins');
%     end
% 
%     if corplot>0; 
% 
%        % correlation of 64 is used for for QA purposes
% 
%     figure;   
%     maxcor1=max(max(cor1))
%     if maxcor1>0;
%      subplot(2,2,1);
%      cor1=cor1*100.0/maxcor1;
%      x=[1 scans];y=[1 bins];
%      color1=jet(100);color1(1,:)=0.0;
%      for colour=limit:top
%         color1(colour,:)=1.0;
%      end;
%      h=image(cor1,'xdata',x,'ydata',y);colormap(color1);colorbar;
%      text(100,1,['Max value  = ',num2str(maxcor1)]);
%      title(' % cor1');grid on;
%      xlabel('Scans');ylabel('bins');
%     end;
% 
%     maxcor2=max(max(cor2));
%     if maxcor2>0;
%      subplot(2,2,2);
%      cor2=cor2*100.0/maxcor2;
%      x=[1 scans];y=[1 bins];
%      color1=jet(100);color1(1,:)=0.0;
%      for colour=limit:top
%         color1(colour,:)=1.0;
%      end;
%      h=image(cor2,'xdata',x,'ydata',y);colormap(color1);colorbar;
%      text(100,1,['Max value  = ',num2str(maxcor2)]);
%      title(' % cor2');grid on;
%      xlabel('Scans');ylabel('bins');
%     end;
% 
%     maxcor3=max(max(cor3));
%     if maxcor3>0
%      subplot(2,2,3);
%      cor3=cor3*100.0/maxcor3;
%      x=[1 scans];y=[1 bins];
%      color1=jet(100);color1(1,:)=0.0;
%      for colour=limit:top
%        color1(colour,:)=1.0;
%      end;
%      h=image(cor3,'xdata',x,'ydata',y);colormap(color1);colorbar;
%      text(100,1,['Max value  = ',num2str(maxcor3)]);
%      title(' % cor3');grid on;
%      xlabel('Scans');ylabel('bins');
%     end;
% 
%     maxcor4=max(max(cor4));
%     if maxcor4>0
%      subplot(2,2,4);
%      maxcor4=max(max(cor4));
%      cor4=cor4*100.0/maxcor4;
%      x=[1 scans];y=[1 bins];
%      color1=jet(100);color1(1,:)=0.0;
%      for colour=limit:top
%         color1(colour,:)=1.0;
%      end;
%      h=image(cor4,'xdata',x,'ydata',y);colormap(color1);colorbar;
%      text(100,1,['Max value  = ',num2str(maxcor4)]);
%      title(' % cor4');grid on;
%      xlabel('Scans');ylabel('bins');
%     end;
% 
%     end;
% 
%     if amps>0
%      figure;
%      subplot(2,2,1);
%       maxamp1=max(max(amp1));
%       minamp1=min(min(amp1));
%       amp1=(amp1-minamp1)*100.0/(maxamp1-minamp1);
%       x=[1 scans];y=[1 bins];
%       color1=jet(100);color1(1,:)=0.0;
%       h=image(amp1,'xdata',x,'ydata',y);colormap(color1);colorbar;
%       text(100,1,['Max value  = ',num2str(maxamp1)]);
%       text(100,3,['Min value  = ',num2str(minamp1)]);
%       title(' Amp1');grid on;
%       xlabel('Scans');ylabel('bins');
% 
%      subplot(2,2,2);
%       maxamp2=max(max(amp2));
%       minamp2=min(min(amp2));
%       amp2=(amp2-minamp2)*100.0/(maxamp2-minamp2);
%       x=[1 scans];y=[1 bins];
%       color1=jet(100);color1(1,:)=0.0;
%       h=image(amp2,'xdata',x,'ydata',y);colormap(color1);colorbar;
%       text(100,1,['Max value  = ',num2str(maxamp2)]);
%       text(100,3,['Min value  = ',num2str(minamp2)]);
%       title(' Amp2');grid on;
%       xlabel('Scans');ylabel('bins');
% 
%      subplot(2,2,3);
%       maxamp3=max(max(amp3));
%       minamp3=min(min(amp3));
%       amp3=(amp3-minamp3)*100.0/(maxamp3-minamp3);
%       x=[1 scans];y=[1 bins];
%       color1=jet(100);color1(1,:)=0.0;
%       h=image(amp3,'xdata',x,'ydata',y);colormap(color1);colorbar;
%       text(100,1,['Max value  = ',num2str(maxamp3)]);
%       text(100,3,['Min value  = ',num2str(minamp3)]);
%       title(' Amp3');grid on;
%       xlabel('Scans');ylabel('bins');
% 
%      subplot(2,2,4);
%       maxamp4=max(max(amp4));
%       minamp4=min(min(amp4));
%       amp4=(amp4-minamp4)*100.0/(maxamp4-minamp4);
%       x=[1 scans];y=[1 bins];
%       color1=jet(100);color1(1,:)=0.0;
%       h=image(amp4,'xdata',x,'ydata',y);colormap(color1);colorbar;
%       text(100,1,['Max value  = ',num2str(maxamp4)]);
%       text(100,3,['Min value  = ',num2str(minamp4)]);
%       title(' Amp4');grid on;
%       xlabel('Scans');ylabel('bins');
%     end  
% 
%     if averageamp>0
%     maxamp=max(max(amp));
%     if maxamp>0
%      figure;
%       amp=amp*100.0/maxamp;
%       x=[first last];y=[1 ts];
%       color1=jet(100);color1(1,:)=0.0;
%       h=image(amp','xdata',x,'ydata',y);colormap(color1);colorbar;
%       text(first+100,1,['Max value  = ',num2str(maxamp)]);
%       title(' Average Amp ');grid on;
%       xlabel('Scans');ylabel('bins');
%     end
% 
%     end
% 
% 
%     if averagecor>0
%     maxcor=max(max(cor));
%     if maxcor>0    
%       figure;
%       cor=cor*100.0/maxcor;
%       x=[first last];y=[1 ts];
%       color1=jet(100);color1(1,:)=0.0;
%       h=image(cor','xdata',x,'ydata',y);colormap(color1);colorbar;
%       text(first+100,1,['Max value  = ',num2str(maxcor)]);
%       title(' Average Cor ');grid on;
%       xlabel('Scans');ylabel('bins');
%     end
% 
%     end
% 
% 
%     if vertical>0
%       figure;
%       maxw=max(max(W));
%       minw=min(min(W));
%       Wm=(W-minw)*100.0/(maxw-minw);
%       maxwm=max(max(Wm));
%       minwm=min(min(Wm));
%       x=[1 scans];y=[1 bins];
%       color1=jet(100);color1(1,:)=0.0;
%       h=image(Wm','xdata',x,'ydata',y);colormap(color1);colorbar;
%       text(100,1,['Max value  = ',num2str(maxw,3)]);
%       text(100,3,['Min value  = ',num2str(minw,3)]);
%       title(' w Verical ');grid on;
%       xlabel('Scans');ylabel('bins');
%     end
% 
% 
%     if Emplot>0
%       figure;
%       maxe=max(max(E));
%       mine=min(min(E));
%       Em=(E-mine)*100.0/(maxe-mine);
%       x=[1 scans];y=[1 bins];
%       color1=jet(100);color1(1,:)=0.0;
%       h=image(Em','xdata',x,'ydata',y);colormap(color1);colorbar;
%       title(' East ');grid on;
%       xlabel('Scans');ylabel('bins');
%     end
% 
%     if Nmplot>0
%       figure;
%       maxn=max(max(N));
%       minn=min(min(N));
%       Nm=(N-minn)*100.0/(maxn-minn);
%       x=[1 scans];y=[1 bins];
%       color1=jet(100);color1(1,:)=0.0;
%       h=image(Nm','xdata',x,'ydata',y);colormap(color1);colorbar;
%       title(' North ');grid on;
%       xlabel('Scans');ylabel('bins');
%     end
% 
% 
% 
%     if tsplot>0;
%      figure;
%      for plots=tsplot:1:ts;
%       subplot(ts+1-tsplot,1,ts-plots+1);
%       [dirn,speed]=cart2pol(E(:,plots),N(:,plots));
%     %  dirn=90.0+deviation-dirn*180/pi;dirn=dirn-360.0e0*round(dirn/360.0e0-0.5e0);
%       plot(mtime,speed,'b-');
%       set(gca,'FontName','Times New Roman');
%       set(gca,'ylim',[0 maxs]);
%       set(gca,'YTick',[0 maxs*0.5 maxs]);
%       i=find(~isnan(speed));
%       chan=speed(i);
%       plots
%       max(chan)
%       text(mtime(1),maxs*0.75,['Mean Speed  bin ',num2str(plots,2),' = ',num2str(mean(chan),2)]);
%       datetick('x');
%       grid on;hold on;
%      end;
%      title(plot_title);
%     end
% 
%     if output>0;
% 
%        for bin=1:ts;
% 
%         outfile=[outname,num2str(bin),'.dat'];
%         fidstring=['fid',num2str(bin+10)];
%         openstring=[fidstring,'=fopen(outfile,''w'');'];
%         closestring=['status=fclose(',fidstring,');'];
%         formatstring=['%5d %3d %3d %3d %3d %4.1f %5d %8.2f %8.2f %8.2f %8.2f %8.2f\n'];
%         varstring1=['ct(1),ct(2),ct(3),ct(4),ct(5),ct(6),scan,'];
%         varstring2=['speed(i),dirn(i),east(i),north(i),deviation'];
%         fprintstring=['fprintf(',fidstring,',''',formatstring,''',',varstring1,varstring2,');'];   
% 
%         eval(openstring);
% 
%      scan=0;
%      east=E(:,bin);north=N(:,bin);
%      [dirn,speed]=cart2pol(east,north);
%      dirn=dirn*180/pi;dirn=90.0-dirn;dirn=dirn+deviation; % rotate sense of 'dirn' from E C/C to N C
%      dirn=dirn-360.0e0*round(dirn/360.0e0-0.5e0); % range dirn within 0>360
%      dirn1=90.0-dirn;dirn1=dirn1*pi/180.0; % reverse sense
%      [east,north]=pol2cart(dirn1,speed); % recompute E & N now true not magnetic
% 
%       for i=1:length(speed);
%        scan=scan+1;
%        ct=datevec(mtime(i));
%        ct(6)=fix(ct(6));
%         eval(fprintstring);
%       end;
% 
%       eval(closestring);
%      end;
% 
%     end
% 
%     if espec>0
%      figure;
%      kplot=0;
%       for i=1:iplot;
%        for j=1:jplot;
%           kplot=kplot+1;subplot(iplot,jplot,kplot);
%           Y=fft(E(:,kplot));M=length(Y);Y(1)=[];power=abs(Y(1:fix(M/2))).^2;nyquist=1/2;
%           freq=(1:fix(M/2))/(M/2)*nyquist;plot(freq,power);axis([0.0 0.5 0 max(power)]);grid on;
%           ylabel(' East Energy');xlabel('frequency');title(['bin ',num2str(kplot)]);
%        end;
%      end;
%     end;
% 
%     if nspec>0
%      figure;
%      kplot=0;
%       for i=1:iplot;
%        for j=1:jplot;
%           kplot=kplot+1;subplot(iplot,jplot,kplot);
%           Y=fft(N(:,kplot));M=length(Y);Y(1)=[];power=abs(Y(1:fix(M/2))).^2;nyquist=1/2;
%           freq=(1:fix(M/2))/(M/2)*nyquist;plot(freq,power);axis([0.0 0.5 0 max(power)]);grid on;
%           ylabel(' North Energy');xlabel('frequency');title(['bin ',num2str(kplot)]);
%        end;
%       end;
%     end;
% 
%     if scatter>0
% 
%      figure;
%      kplot=0;
%      set(gca,'FontName','Times New Roman'); 
%       for i=1:iplot;
%         for j=1:jplot;
%          kplot=kplot+1;subplot(iplot,jplot,kplot);
% 
%           east=E(:,kplot);north=N(:,kplot);
%           [dirn,speed]=cart2pol(east,north);
%           dirn=dirn*180/pi;dirn=90.0-dirn;dirn=dirn+deviation; % rotate sense of 'dirn' from E C/C to N C
%           dirn=dirn-360.0e0*round(dirn/360.0e0-0.5e0); % range dirn within 0>360
%           dirn1=90.0-dirn;dirn1=dirn1*pi/180.0; % reverse sense
%           [east,north]=pol2cart(dirn1,speed); % recompute E & N now true not magnetic
% 
%          plot(east,north,'b+','Markersize',1);grid on;axis equal;
%          set(gca,'ylim',[-maxs maxs]);
%          set(gca,'xlim',[-maxs maxs]);
%          set(gca,'ytick',[-maxs -maxs/2 0 maxs/2 maxs]);
%          set(gca,'xtick',[-maxs -maxs/2 0 maxs/2 maxs]);
%          xlabel('east cm/s');
%          ylabel('north cm/s');
%          title([' bin ',num2str(kplot)]);
%         end;
%       end;
%     end;
% 
% 
%     if histogram>0
%      figure;
%      kplot=0;
%       for i=1:iplot;
%         for j=1:jplot;
%          kplot=kplot+1;subplot(iplot,jplot,kplot);
%          east=E(:,kplot);north=N(:,kplot);
%          [dirn,speed]=cart2pol(east,north);
%       %   dirn=90.0+deviation-dirn*180/pi;dirn=dirn-360.0e0*round(dirn/360.0e0-0.5e0);
%          hist(speed,150);xlabel('speed cm/s');
%          set(gca,'FontName','Times New Roman'); 
%          set(gca,'ylim',[0 2500]);
%          set(gca,'xlim',[0 maxs]);
%          set(gca,'ytick',[0 500 1000 1500 2000 2500]);
%          set(gca,'xtick',[0 maxs/2 maxs]);
%          title([' bin ',num2str(kplot)]);
%         end;
%      end;
% 
%      figure;
%      maxd=360;
%      kplot=0;
%       for i=1:iplot;
%         for j=1:jplot;
%          kplot=kplot+1;subplot(iplot,jplot,kplot);
%          east=E(:,kplot);north=N(:,kplot);
%          [dirn,speed]=cart2pol(east,north);
%          dirn=90.0+deviation-dirn*180/pi;dirn=dirn-360.0e0*round(dirn/360.0e0-0.5e0);
%          hist(dirn,360);xlabel('Direction (True)');
%          set(gca,'FontName','Times New Roman'); 
%          set(gca,'ylim',[0 300]);
%          set(gca,'xlim',[0 maxd]);
%          set(gca,'ytick',[0 100 200 300]);
%          set(gca,'xtick',[0 90 180 270 360]);
%          title([' bin ',num2str(kplot)]);
%         end;
%      end;
%     end;
% 
% 
% 
%     if pvd>0
%      figure;
%      kplot=0;
%        scale=interval*60*60/100000; % convert scale from cm/s to km/day - interval is in hours
%       for i=1:iplot;
%         for j=1:jplot;
%          kplot=kplot+1;subplot(iplot,jplot,kplot);
%          tdistx=[];tdisty=[];distx=0.0;disty=0.0;
%          text(0.0,0.0,'Start');grid on;hold on;
%          xlabel('West to East (km)');
%          ylabel('South to North (km)');
%          title([' PVD bin ',num2str(kplot)]);
% 
%          east=E(:,kplot);north=N(:,kplot);
%          [dirn,speed]=cart2pol(east,north);
%          dirn=dirn*180/pi;dirn=90.0-dirn;dirn=dirn+deviation; % rotate sense of 'dirn' from E C/C to N C
%          dirn=dirn-360.0e0*round(dirn/360.0e0-0.5e0); % range dirn within 0>360
%          dirn1=90.0-dirn;dirn1=dirn1*pi/180.0; % reverse sense (Mathematical direction)
%          [east,north]=pol2cart(dirn1,speed); % recompute E & N, now true, not magnetic
% 
%          for n=1:good;
%            distx=distx+east(n)*scale;
%            disty=disty+north(n)*scale;
%            tdistx=[tdistx distx];
%            tdisty=[tdisty disty];
%           end;
%          plot(tdistx,tdisty);
%          clear tdistx tdisty distx disty east north;
%         end 
%        end;
% 
%     end;
% 
% 
%     if ex>0 % just o/p exceedances for bin 2
%     one=1;
%     zero=0;
%     fid21=fopen([mooring,'bin2.hst'],'w');
%     for bin=2:2;
%      east=E(:,bin);north=N(:,bin);
%      [dirn,speed]=cart2pol(east,north);
%      maxspeed=fix(max(speed))+1;
%      l=length(speed);total=0;
% 
%       for index=1:maxspeed
%         sh(index)=0.0;cum(index)=0.0;
%       end;
% 
%       for i=1:l
%          if speed(i) > 0.00;
%          ind=fix(speed(i)+1);
%          sh(ind)=sh(ind)+1;total=total+1;
%          end;
%       end;
%      per(1)=100.0*sh(1)/l;ex(1)=per(1);cum(1)=1.0-sh(1)/l;
%      fprintf(fid21,'%4d %4d %5d %7.3f %8.2f %10.8f\n',...
%           bin,one,sh(1),per(1),ex(1),cum(1));
%       for i=2:maxspeed;
%         per(i)=100.0*sh(i)/l;
%         cum(i)=cum(i-1)-sh(i)/l;
%         ex(i)=per(i)+ex(i-1);
%         fprintf(fid21,'%4d %4d %5d %7.3f %8.2f %10.8f\n',...
%           bin,i,sh(i),per(i),ex(i),cum(i));
%       end;
%     end;
% 
% 
%     status=fclose(fid21);
% 
%     end

    if percentG>0
     figure('position',get(0,'ScreenSize')); 
     sp1=subplot(2,1,1)
     hold on; 
         %maxpg4=max(max(pg4));
         %pg4=pg4*100.0/maxpg4;
         x=[1 scans];y=[1 bins];
         Vvel = N';
         Vvel(pg1+pg4<70)=nan;
         h=pcolor(scannber,SerBins,Vvel);colormap(sp1,jet(100));colorbar;
         set(gca,'clim',[-50 50])
         %text(100,1,['Max value  = ',num2str(maxpg4)]);
         title(' % V component');grid on;
         xlabel('Scans');ylabel('bins');
         xlimpg4 = get(gca,'xlim');
         ylimpg4 = get(gca,'ylim');
         set(gca,'ydir','reverse')
         set(gca,'xlim',[xlimpg4(1)-100 xlimpg4(2)+100]);
         set(gca,'ylim',[ylimpg4(1)-5 ylimpg4(2)+5]);     
        shading interp
  
     sp2=subplot(2,1,2)
         hold on; 
         %maxpg4=max(max(pg4));
         %pg4=pg4*100.0/maxpg4;
         x=[1 scans];y=[1 bins];
         color1=jet(100);color1(1,:)=0.0;
         for colour=75:top
            color1(colour,:)=1.0;
         end;
         h=image(pg4 + pg1,'xdata',x,'ydata',y);colormap(sp2,color1);colorbar;
         %text(100,1,['Max value  = ',num2str(maxpg4)]);
         title(' % Good pg1 + pg4');grid on;
         xlabel('Scans');ylabel('bins');
         xlimpg4 = get(gca,'xlim');
         ylimpg4 = get(gca,'ylim');
         set(gca,'ydir','reverse')
         set(gca,'xlim',[xlimpg4(1)-100 xlimpg4(2)+100]);
         set(gca,'ylim',[ylimpg4(1)-5 ylimpg4(2)+5]);
         
         
      set(gcf,'CurrentAxes',sp1)
      
         disp(' '); disp('Select area that you want to be flagged as bad'); disp(' ');    
         matpgfourfile = [outpath outfile '_manualflag.mat'];
         
         vardescription.mtime = 'time, matlabformat';
         vardescription.pg4   = 'percent good 4: The percentage of measurements with four-beam solutions';
         vardescription.pg1   = 'percent good 1: The percentage of measurements with three-beam solutions';    
         vardescription.badarea = 'Outside this area the velocity value are set to nan, in order to not have gaps in the different bin time-series';       

         
         if exist(matpgfourfile)==2
            ierase=input('\n A mat file already exists, do you want to erase it? y/n:  ','s');
            if (strcmp(ierase,'y')||strcmp(ierase,'Y')||strcmp(ierase,'yes')||strcmp(ierase,'Yes')||strcmp(ierase,'YES'))
                badarea = ginput;
                ibaddata= ~inpolygon(repmat(scannber,length(binnber),1),repmat(binnber,1,length(scannber)),badarea(:,1),badarea(:,2));
                save(matpgfourfile, 'ibaddata','scannber','binnber','mtime','pg4','pg1','N');

                plot(badarea(:,1),badarea(:,2),'-r');          
            else
                continue
            end
         else
             badarea = ginput;
             ibaddata= ~inpolygon(repmat(scannber,length(binnber),1),repmat(binnber,1,length(scannber)),badarea(:,1),badarea(:,2));
             save(matpgfourfile, 'ibaddata','scannber','binnber','mtime','pg4','pg1','N');

             plot(badarea(:,1),badarea(:,2),'-r');    
   
         end
         
       set(gcf,'CurrentAxes',sp2)  
       plot(badarea(:,1),badarea(:,2),'-r');    

             print(gcf,'-dpng',[outpath outfile])          
    end
end 


end

