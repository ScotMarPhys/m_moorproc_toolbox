% Find the anchor position from more than 2 positions and ranges
% run from the same directory as the position info file
% position info file contains 3 columns. One for range (m), one
% for lat (dd.dddd) and one for long (dd.dddd).

% modified by Loic Houpert 26/06/2015, add print epsc and png, remove the
% calculation of the corrected water depth (at the anchor drop) in the script. The user will
% have to give the corrected water depth in the input data file.
% by using for example the function get_position_depth_pe399
% Modified for DY053, June 2016

%modified by Loic Houpert and James Coogan 08/07/2018, added distance to
%target depth.

clear all ; close all ;

% addpath('../../carter_corr') 
addpath('/data/pstar/programs/gitvcd/m_moorproc_toolbox/other_software/carter_corr/')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA:
% first row is target position
% second row is uncorrected water depth at anchor drop : Release height above bottom : transducer depth
% third row is position at anchor drop
[posfile1,posfile1_pathname] = uigetfile('*.*','Load trilateration file : ');
cd(posfile1_pathname);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% posfile1=input('Enter name of position file\n','s');
% wd=input('Enter uncorrected water depth at anchor drop in metres\n');
% rht=input('Enter approximate release height above seabed in metres\n');
% td=input('Enter approximate transducer depth in metres\n');

posfile=load(posfile1);
no_fixes=size(posfile,1);
lat_tgt = posfile(1,2);
lon_tgt = posfile(1,3);
wd = posfile(2,1); rht = posfile(2,2); td = posfile(2,3) ;

lat_drop = posfile(3,2) ; lon_drop = posfile(3,3) ;
wd_corr = carter(lat_drop,lon_drop,wd);

for i=4:no_fixes
    range(i) = posfile(i,1);
    lat(i) = posfile(i,2);
    lon(i) = posfile(i,3);
    range(i) = carter(lat(i),lon(i),range(i)); % This is correct but not that adams slant ranges are already corrected so this should not be applied to his data
    rangeh(i) = sqrt((range(i))^2 - (wd_corr-rht-td)^2); % horizontal range from ship to position vertically above release
    
end


figure(1);clf;hold on ;grid on ;
plot(lon_drop,lat_drop,'k^','markerfacecolor','b','markersize',5)
plot(lon_tgt,lat_tgt,'ks','markersize',5)
xdegr = -180:2:180;
for i=4:no_fixes
    [clon,clat,az3] = m_fdist(lon(i),lat(i),xdegr,rangeh(i));
    clon = clon - 360;
    plot(clon,clat,'b');   hold on ; grid on ;
    plot(lon(i),lat(i),'b.','markersize',10)
end
c=0.05;

axis([min(lon(4:end))-c,max(lon(4:end))+c,min(lat(4:end))-c,max(lat(4:end))+c]);
scale=1/cos((pi()/180)*mean(lat(4:end)));
daspect([scale 1 1]);
titletext1=['Trilateration Survey using: ',posfile1];
titletext2=['Corrected water depth: ',num2str(wd_corr),'m. Release Height: ',num2str(rht),'m. Transducer depth: ',num2str(td),'m.'];
title({titletext1;titletext2},'interpreter','none');
xlabel('longitude'); ylabel('latitude')


xlims=get(gca,'xlim');
ylims=get(gca,'ylim');


% determine anchor position from figure(1)
disp('Use figure 1 to determine anchor seabed position, step 1: zoom in on fig to find position, then hit return')
pause;
disp('Use figure 1 to determine anchor seabed position, step2: click on position in fig')
[APlon,APlat] = ginput(1);
%APlat = input('Latitude = ');
%APlon = input('Longitude = ');

dist_tgt = sw_dist([lat_tgt,APlat], [lon_tgt, APlon], 'km')*1000;

plot(APlon,APlat,'r+');
plot(APlon,APlat,'ro');
titletext3=['Red = anchor seabed position. ',num2str(APlat,'%2.6f'),'N ',num2str(APlon,'%2.6f'),'W.'];
titletext4= ['Distance to target: ',num2str(dist_tgt,'%5.0f'),' m'];
title({titletext1;titletext2;titletext3;titletext4});

set(gca,'xlim',xlims)
set(gca,'ylim',ylims)


fallback=sw_dist([lat_drop APlat],[lon_drop APlon],'km')*1000;
fallback=sprintf('%5.0f',fallback);
disp(['Fallback :',fallback,' m']);
disp(['Distance to target:', num2str(dist_tgt,'%5.0f'),' m']);
text((xlims(2)-xlims(1))*0.75+xlims(1),(ylims(2)-ylims(1))*0.95+ylims(1),['Fallback = ' fallback 'm'],'color','r')

print('-depsc2', [posfile1(1:end-4) '_results'])
print('-dpng',[posfile1(1:end-4) '_results'])




