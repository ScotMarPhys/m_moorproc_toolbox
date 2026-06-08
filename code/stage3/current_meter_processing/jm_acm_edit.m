% get rodb_file with raw data format als variables
% added pressure into rodb header -- JM 2/05
% this script also uses jm_acm_edit_header function which reinserts header
% in the data file at the end of this script (no need to manually insert
% *.conf file which contains header file)
% removed w --JM 12/2005

close all;
clear all;

ifile={'37102001.COR',...
        '37104001.COR'};

ofile={'37102001.edt',...
        '37104001.edt'};
    
confile={'37102001.conf',...
        '37104001.conf'};
        
filfac=10;      % --- filter factor for despiking


for mr=1:length(ifile);
fidi=[ifile{mr}];


%[fn,moor,inst,sn,insdp,wd,lat,lon,magdev,yy,mm,dd,hh,t,u,v,w]=rodbload(ifile,...
%       'Filename:Mooring:Instrument:Serial_Number:InstrDepth:WaterDepth:Latitude:Longitude:Mag_Deviation:yy:mm:dd:hh:t:u:v:w');

[fn,moor,yy,mm,dd,h,t,p,u,v]=jm_rodbload(fidi,...
       'Filename:Mooring:yy:mm:dd:hh:t:p:u:v');


% make time to julian;

time=julian(yy,mm,dd,h);
dt=(time(2)-time(1))*24    % sampling interval
nt=round(40/dt);            % 40 hr lowpass

% edit data for dummies before filtering
iu=isnan(u);
u(iu)=meanmiss(u);
iv=isnan(v);
v(iv)=meanmiss(v);
ip=isnan(p);
p(ip)=meanmiss(p);
it=isnan(t);
t(it)=meanmiss(t);


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
tim=[fix(time(1)+0.5):0.5:fix(time(end))];

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
end
