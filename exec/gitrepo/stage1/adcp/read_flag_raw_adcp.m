% function read_flag_raw_adcp(moor,varargin)
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

global MOORPROC_G

basedir = MOORPROC_G.moordatadir;
cruise= MOORPROC_G.cruise;

if nargin==0
    help read_flag_raw_adcp
    return
end

if nargin==1
    pd = moor_inoutpaths('adcp',moor);
else
    pd = varargin{1};
end

operator = MOORPROC_G.operator;

% ----- read infofile / open logfile  ------------------------------------

infovar = 'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth'; 
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd]  =  rodbload(pd.infofile,infovar);

fidlog   = fopen(pd.stage1log,'a');
fprintf(fidlog,'Transformation of ADCP .mat data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));
fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n\n\n',lon);

bg = julian([s_d(:)' hms2h([s_t(:)' 0])]); %start
ed = julian([e_d(:)' hms2h([e_t(:)' 0])]); %end



vec=find((id>=319) & (id <=328)); % Possible ADCP codes - taken from IMP moorings package

serial_nums=sn(vec)

% -------- load data --------------
for i = 1:length(vec)

    fprintf('Processing sn %d',serial_nums(i))
    infile=fullfile(pd.rawpath,sprintf('%d_data.mat',serial_nums(i)));
    
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


    if percentG>0
     figure('units','normalized','outerposition',[0.01 0.01 .99 .99]); 
     sp1=subplot(2,1,1);
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
  
     sp2=subplot(2,1,2);
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
         matpgfourfile = fullfile(pd.stage1path,[outfile '_manualflag.mat']);
         
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

             print(gcf,'-dpng',fullfile(pd.stage1path,outfile))         
    end
end 
close all
end

