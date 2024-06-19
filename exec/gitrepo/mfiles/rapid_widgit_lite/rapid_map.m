function rapid_map(this_moor,mlat,mlon,mdepth,north,south,east,west,cont_int)
% This routine plots matp
figure

% Set up m_map plot
m_proj('lambert','lon',[west, east],'lat',[south, north]);
m_grid('box','on','color','k','linewidth',[1],'fontsize',[12]);
hold on

% Using user supplied mooring depth unless otherwise chagned
dpfl = 'user';

% load data covering region
dd = 0.1*abs(north-south);

% See if we have swath data - defaul t is yes
iswath = 1;
if strmatch(this_moor,char([{'mar3','mar3l','mar3l1','mar3l2','nog'}])) > 0
	topo_map = 'mar34.mat';
elseif strmatch(this_moor,char([{'mar2'},{'mar1'},{'mar1l'},{'mar1l1'},{'mar1l2'}])) > 0
	topo_map = 'mar12.mat';
elseif strmatch(this_moor,char([{'mar0'}])) > 0
% topo_map = 'mar0_JC064_swath.mat';
    topo_map = 'mar0_jc103.mat';
elseif strmatch(this_moor,char([{'wbadcp'},{'wbal'},{'wbal1'},{'wbal2'},{'wb1'}, ...
                 {'wb2'},{'wb2l'},{'wb2l1'},{'wb2l2'},{'wbh2'}])) > 0
    topo_map = 'gr_kn182_plot.mat';
else
	iswath = 0;
    fprintf(1,'Using Smith and Sandwell \n')
    [elevations,elat,elon]=mygrid_sand([south-dd,north+dd,west-dd,east+dd],1);
    elon = elon - 360;
end

% If swath data available lets use it
if iswath == 1
	fprintf(1,'Using swath data %s \n',topo_map)
    load(['/noc/users/pstar/rpdmoc/rapid/data/exec/jc103/mfiles/rapid_widgit_v2/data/' topo_map]);
	ixlo = lon > west-dd & lon < east+dd;
	iyla = lat > south-dd & lat < north+dd;
    elon = lon(ixlo);
    elat = lat(iyla); 
    elevations = -depth(iyla,ixlo);
    mdepthi = interp2(lon,lat,depth,mlon,mlat);
	if ~isnan(mdepthi)
		dpfl = 'swath';
		mdepth = mdepthi;
	end
    clear lat lon depth;
end

mxdp = max(max(-elevations));
mndp = min(min(-elevations));

if strmatch(this_moor,'ship') == 1
 mdepth = 0.5*(mxdp+mndp);
end

if mdepth > mxdp | mdepth < mndp
 disp('Somthing not right about depths')
end

% Plot some contours but first work out best interval
mdepth0 = 0.5*(mxdp+mndp);
mdepth0 = -elevations(floor(size(elevations,1)/2),floor(size(elevations,2)/2));
mdepth0 = 10*floor(mdepth/10);
if cont_int < 0 
	cnti = abs(mxdp-mndp)/8;
	cnti = 10*floor(cnti/10);
	cnti = max(cnti,10);
	cnti = min(cnti,200);
else
	cnti = cont_int;
end

ncn1 = ceil((mxdp-mdepth0)/cnti);
ncn2 = ceil((mdepth0-mndp)/cnti);
conts = -mdepth0 + cnti*[-ncn1:0 0:ncn2];
m_contour(elon,elat,elevations,conts(1:ncn1),'b--','LineWidth',0.75);
m_contour(elon,elat,elevations,conts(ncn1+1:ncn1+2),'c--','LineWidth',1.25);
m_contour(elon,elat,elevations,conts(ncn1+3:ncn1+ncn2+2),'g--','LineWidth',0.75);

fprintf(1,'Water depth at centre contour is %6.1f m \n',mdepth0)
fprintf(1,' ** Water depth at mooring is %6.1f m from %s ** \n',mdepth,dpfl)
fprintf(1,'Depth range %6.1f to %6.1f \n',mndp,mxdp)

% Plot text below plot
cnttxt = sprintf('Contour interval = %5.1f m \nCyan contour depth = %6.1f m',cnti,mdepth0);
xt = west  + 0.15*(east-west);
yt = south - 0.2*(north-south);
hb = m_text(xt,yt,cnttxt,'FontSize',12);

% Plot location of mooring
ha = m_plot(mlon,mlat,'k+');
set(ha,'MarkerSize',12,'LineWidth',1.5)

%CContour info
cnttxt = sprintf('Contour interval = %5.1f m \nCyan contour depth = %6.1f m',cnti,mdepth0);
xt = west  + 0.15*(east-west);
yt = north - 0.15*(north-south);
hb = m_text(xt,yt,cnttxt,'FontSize',12);
fprintf(1,'Hello world')
keyboard