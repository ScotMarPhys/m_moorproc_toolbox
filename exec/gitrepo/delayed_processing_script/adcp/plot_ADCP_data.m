% Function to load and plot ADCP data at the rodb file format 
% 
%       by loh, SAMS, 15/10/2015

close all, clear all

%----------------------------------------------------
% Graphical parameter
vsblfig = 'off';
zsc = get(0,'MonitorPositions');
scrsz = [1 1 1900 1400]; % zsc(1,:);
figpos = scrsz/60; %[0 0 27 21];
fs1 = 12;
fs2=10;
set(0,'DefaultAxesFontName', 'Helvetica')
set(0,'DefaultAxesFontSize', fs1)
set(0,'DefaultTextFontname', 'Helvetica')
set(0,'DefaultTextFontSize', fs2)
%-----------------------------------------------------
%-----------------------------------------------------
% Mooring parameters
cruise   = 'pe399'; %'pe399';
operator = 'LoicH';
%moor  = 'rtadcp1_01_2014';
moor ='ib5_01_2018'
datatype = 'lowedt' %'use';%'lowedt'; %"edt" for no filtered data and "lowedt" for lowpass data 


basedir  = '~/osnap/data/moor/';
procpath = [basedir 'proc/'];
% outpath  = [procpath mooring '/adcp/'];
outpath  = [procpath moor '/adp/'];

cd(outpath)
% -- set path for data input
inpath  = [procpath moor '/'];

% --- get moring information from infofile 
infofile =[procpath moor '/' moor 'info.dat'];
%-----------------------------------------------------
%-----------------------------------------------------
% Load data

[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

vec=find((id>=319) & (id <=328)) % Possible ADCP codes - taken from IMP moorings package

sn=sn(vec);
z = z(vec);
zadcp = z;
cmapvel    = jet;%cbrewer('div','RdYlBu',20);

for proc = 1 : length(vec) % loop for multiple instruments on moooring
    
    % first determine how many bin files are to be processed
    % trying to do automatically
    num_bins=dir([inpath,'adp/*',num2str(sn(proc)) '*.use'])
    num_bins=length(num_bins)
  
    disp(['ADCP Serial number: ' num2str(sn(proc))])

    
    num_bins
    %---------------------------------
    % definition length of the files
    j=2;
    columns = ['YY:MM:DD:HH:Z:T:U:V:W:HDG:PIT:ROL:CS:CD:BEAM1SS:BEAM2SS:BEAM3SS'...
            ':BEAM4SS:BEAM1COR:BEAM2COR:BEAM3COR:BEAM4COR:EV:BEAM1PGP:BEAM2PGP:BEAM3PGP:BEAM4PGP'];
        indep  = zadcp(proc);    
        infile  = [inpath,'adp/',moor,'_',num2str(sn(proc)),'_bin0',num2str(j),'.',datatype];
        rodbfile= [moor,'_',num2str(sn(proc)),'_bin0',num2str(j),'.',datatype]; 
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
        
        if j<=9
            infile  = [inpath,'adp/',moor,'_',num2str(sn(proc)),'_bin0',num2str(j),'.',datatype];
        else
            infile  = [inpath,'adp/',moor,'_',num2str(sn(proc)),'_bin',num2str(j),'.',datatype];
        end
        if exist(infile,'file')==0
            disp(['infile: ' infile ' does not exist.'])
 
        elseif exist(infile,'file')   > 0 
            if j<=9
                rodbfile= [moor,'_',num2str(sn(proc)),'_bin0',num2str(j),'.',datatype]; 
            else
                rodbfile= [moor,'_',num2str(sn(proc)),'_bin',num2str(j),'.',datatype]; 
            end
            outfile = [outpath,rodbfile];


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


		fig=figure('visible',vsblfig,'position',scrsz);
		set(fig,'PaperUnits','centimeters','PaperOrientation','portrait',... 
				    'Paperposition',figpos)
subplot(3,1,1)
surf(ADCPdata(proc).time, ADCPdata(proc).z,ADCPdata(proc).u);
shading flat
view([0 -90])
caxis([-45 45])
title('U (cm.s-1)')
datetick
colormap(cmapvel)
colorbar


subplot(3,1,2)
surf(ADCPdata(proc).time, ADCPdata(proc).z,ADCPdata(proc).v);
shading flat
view([0 -90])
caxis([-45 45])
title('V (cm.s-1)')
datetick
colormap(cmapvel)
colorbar

if isfield(ADCPdata,'PG4')
subplot(3,1,3)
surf(ADCPdata(proc).time, ADCPdata(proc).z,ADCPdata(proc).PG4);
shading flat
view([0 -90])
caxis([0 100])
title('BEAM4PGP')
datetick
colormap(gca,jet)
colorbar
end

		print(gcf,'-dpng',['adcp_' datatype '_data_' moor '_' num2str(sn(proc))]);	
    
save(['adcp_' datatype '_data_' moor '_' num2str(sn(proc)) '.mat'], 'ADCPdata') 
end


