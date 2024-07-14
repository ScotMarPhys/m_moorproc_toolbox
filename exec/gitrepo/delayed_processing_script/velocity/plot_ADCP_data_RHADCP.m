% Function to load and plot ADCP data at the rodb file format 
% 
%       by loh, SAMS, 15/10/2015
%       by loh, 20/10/20 add upward and downward looking on same plot

close all
global MOORPROC_G
clearvars -except MEXEC_G MOORPROC_G
cruise   = MOORPROC_G.cruise;
operator = MOORPROC_G.operator;
moor = input('mooring deployment (e.g. ebh2_15_2022) to process:   ','s');

pd = moor_inoutpaths('adcp',moor);

datatype = 'lowedt' %'use';%'lowedt'; %"edt" for no filtered data and "lowedt" for lowpass data 

cd(pd.stage3path)

% --- get moring information from infofile 
infofile =pd.infofile;

%----------------------------------------------------
% Graphical parameter
vsblfig = 'on';
zsc = get(0,'MonitorPositions');
scrsz = [1 1 1600 1200]/2; % zsc(1,:);
figpos = scrsz/30; %[0 0 27 21];
fs1 = 12;
fs2=10;
set(0,'DefaultAxesFontName', 'Helvetica')
set(0,'DefaultAxesFontSize', fs1)
set(0,'DefaultTextFontname', 'Helvetica')
set(0,'DefaultTextFontSize', fs2)
%-----------------------------------------------------
%-----------------------------------------------------

%-----------------------------------------------------
%-----------------------------------------------------
% Load data

[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

vec=find((id>=319) & (id <=328)) % Possible ADCP codes - taken from IMP moorings package

sn=sn(vec);
z = z(vec);
zadcp = z;

cmap =[ ...
          [   0         0         1.0000 ];
          [   0.0313    0.0323    0.9995 ];
          [   0.0626    0.0645    0.9990 ];
          [   0.0938    0.0968    0.9985 ];
          [   0.1251    0.1290    0.9980 ];
          [   0.1564    0.1613    0.9976 ];
          [   0.1877    0.1935    0.9971 ];
          [   0.2190    0.2258    0.9966 ];
          [   0.2502    0.2581    0.9961 ];
          [   0.2815    0.2903    0.9956 ];
          [   0.3128    0.3226    0.9951 ];
          [   0.3441    0.3548    0.9946 ];
          [   0.3754    0.3871    0.9941 ];
          [   0.4066    0.4194    0.9936 ];
          [   0.4379    0.4516    0.9932 ];
          [   0.4692    0.4839    0.9927 ];
          [   0.5005    0.5161    0.9922 ];
          [   0.5318    0.5484    0.9917 ];
          [   0.5630    0.5806    0.9912 ];
          [   0.5943    0.6129    0.9907 ];
          [   0.6256    0.6452    0.9902 ];
          [   0.6569    0.6774    0.9897 ];
          [   0.6882    0.7097    0.9892 ];
          [   0.7195    0.7419    0.9888 ];
          [   0.7507    0.7742    0.9883 ];
          [   0.7820    0.8065    0.9878 ];
          [   0.8133    0.8387    0.9873 ];
          [   0.8446    0.8710    0.9868 ];
          [   0.8759    0.9032    0.9863 ];
          [   0.9071    0.9355    0.9858 ];
          [   0.9384    0.9677    0.9853 ];
          [   0.9697    1.0000    0.9848 ];
          [   0.9707    0.9677    0.9531 ];
          [   0.9717    0.9355    0.9213 ];
          [   0.9726    0.9032    0.8895 ];
          [   0.9736    0.8710    0.8578 ];
          [   0.9746    0.8387    0.8260 ];
          [   0.9756    0.8065    0.7942 ];
          [   0.9765    0.7742    0.7625 ];
          [   0.9775    0.7419    0.7307 ];
          [   0.9785    0.7097    0.6989 ];
          [   0.9795    0.6774    0.6672 ];
          [   0.9804    0.6452    0.6354 ];
          [   0.9814    0.6129    0.6036 ];
          [   0.9824    0.5806    0.5718 ];
          [   0.9834    0.5484    0.5401 ];
          [   0.9844    0.5161    0.5083 ];
          [   0.9853    0.4839    0.4765 ];
          [   0.9863    0.4516    0.4448 ];
          [   0.9873    0.4194    0.4130 ];
          [   0.9883    0.3871    0.3812 ];
          [   0.9892    0.3548    0.3495 ];
          [   0.9902    0.3226    0.3177 ];
          [   0.9912    0.2903    0.2859 ];
          [   0.9922    0.2581    0.2542 ];
          [   0.9932    0.2258    0.2224 ];
          [   0.9941    0.1935    0.1906 ];
          [   0.9951    0.1613    0.1588 ];
          [   0.9961    0.1290    0.1271 ];
          [   0.9971    0.0968    0.0953 ];
          [   0.9980    0.0645    0.0635 ];
          [   0.9990    0.0323    0.0318 ];
          [   1.0000         0         0 ];
    ];
cmapvel = cmap; %jet;
    
for proc = 1 : length(vec) % loop for multiple instruments on moooring
    
    % first determine how many bin files are to be processed
    % trying to do automatically
    num_bins=dir(fullfile(pd.stage2path,...
        [sprintf(pd.stage2inform,sn(proc)),'*.',datatype]));
    num_bins=length(num_bins);
  
    disp(['ADCP Serial number: ' num2str(sn(proc))])

    %---------------------------------
    % definition length of the files
    j=2;
    columns = ['YY:MM:DD:HH:Z:T:U:V:W:HDG:PIT:ROL:CS:CD:BEAM1SS:BEAM2SS:BEAM3SS'...
            ':BEAM4SS:BEAM1COR:BEAM2COR:BEAM3COR:BEAM4COR:EV:BEAM1PGP:BEAM2PGP:BEAM3PGP:BEAM4PGP'];
        indep  = zadcp(proc);    
        infile  = fullfile(pd.stage3path,...
                    sprintf([moor '_%d_bin%02.f.%s'],sn(proc),j,datatype)); 
        [YY,MM,DD,HH,z,t,u,v,w,heading,pitch,roll,spd,direction,Amp1,Amp2,Amp3,Amp4,...
                Beam1Cor,Beam2Cor,Beam3Cor,Beam4Cor,err,PG1,PG2,PG3,PG4] = ...
                rodbload(infile,[columns]);    
    
    ADCPdata(proc).time=nan(num_bins,length(z));
    ADCPdata(proc).z=nan(num_bins,length(z));   
    ADCPdata(proc).t=nan(num_bins,length(z)); 
    ADCPdata(proc).u=nan(num_bins,length(z));
    ADCPdata(proc).v=nan(num_bins,length(z));   
    ADCPdata(proc).w=nan(num_bins,length(z));        
    ADCPdata(proc).heading=nan(num_bins,length(z));
    ADCPdata(proc).pitch=nan(num_bins,length(z));   
    ADCPdata(proc).roll=nan(num_bins,length(z));    
    %---------------------------------   
    
    namefigadcp1 = [moor '_adcp_' num2str(proc,'%02d') '_uv_surf'];
    namefigadcp2 = [moor '_adcp_' num2str(proc,'%02d') '_pgood_surf'];    
    
    
    for j=1:num_bins % loop for total number of bins
        
        columns = ['YY:MM:DD:HH:Z:T:U:V:W:HDG:PIT:ROL:CS:CD:BEAM1SS:BEAM2SS:BEAM3SS'...
            ':BEAM4SS:BEAM1COR:BEAM2COR:BEAM3COR:BEAM4COR:EV:BEAM1PGP:BEAM2PGP:BEAM3PGP:BEAM4PGP'];
        
        infile  = fullfile(pd.stage3path,...
                    sprintf([moor '_%d_bin%02.f.%s'],sn(proc),j,datatype));
        
        if exist(infile,'file')==0
            disp(['infile: ' infile ' does not exist.'])
        
        else
            [YY,MM,DD,HH,z,t,u,v,w,heading,pitch,roll,spd,direction,Amp1,Amp2,Amp3,Amp4,...
                Beam1Cor,Beam2Cor,Beam3Cor,Beam4Cor,err,PG1,PG2,PG3,PG4] = ...
                rodbload(infile,[columns]);
            
            ADCPdata(proc).time(j,:) = datenum(YY,MM,DD,HH,0*HH,0*HH);
            ADCPdata(proc).z(j,:) = z;
            ADCPdata(proc).t(j,:) = t;          
            ADCPdata(proc).u(j,:) = u;
            ADCPdata(proc).v(j,:) = v;       
            ADCPdata(proc).w(j,:) = w;     
          
         switch(datatype)
         	case('use')
            		ADCPdata(proc).PG1(j,:) = PG1;       
            		ADCPdata(proc).PG2(j,:) = PG2;    
            		ADCPdata(proc).PG3(j,:) = PG3;       
            		ADCPdata(proc).PG4(j,:) = PG4;
         end 
         
                        
        end 
 end
    
    
switch(datatype)
    case('use')   
ibad= ADCPdata(proc).PG4<40 ;
ADCPdata(proc).u(ibad)=nan;
ADCPdata(proc).v(ibad)=nan;
ADCPdata(proc).w(ibad)=nan;
end 

ADCPdata(proc).u(ADCPdata(proc).u==-9999)=nan;
ADCPdata(proc).v(ADCPdata(proc).v==-9999)=nan;
ADCPdata(proc).w(ADCPdata(proc).w==-9999)=nan;
    
    
% save(['adcp_' datatype '_data_' moor],'ADCPdata')

end

fig=figure('visible',vsblfig,'position',scrsz);
set(fig,'PaperUnits','centimeters','PaperOrientation','portrait',... 
            'Paperposition',figpos)
subplot(2,1,1)
hold on
surf(ADCPdata(1).time, ADCPdata(1).z,ADCPdata(1).u);
shading flat
view([0 -90])
caxis([-100 100])
title(['U -- ' moor ' ADCPs '  num2str(sn) ' (up) '],'interpreter','none'  )
datetick('x','mmm yy')
colormap(cmapvel)
cb = colorbar;
cb.Label.String='cm/s';
ylabel('depth (m)')

subplot(2,1,2)
hold on
surf(ADCPdata(1).time, ADCPdata(1).z,ADCPdata(1).v);
shading flat
view([0 -90])
caxis([-100 100])
title(['V -- ' moor ' ADCPs '  num2str(sn(1)) ' (up) ' ],'interpreter','none'  )
datetick('x','mmm yy')
colormap(cmapvel)
cb = colorbar;
cb.Label.String='cm/s';
ylabel('depth (m)')

print(gcf,'-dpng',['adcp_' datatype '_data_' moor '_up']);	

