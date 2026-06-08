% function merge_TS_current_eddy(moor,'procpath','proclvl','layout','plot_interval','unfiltered')
%

%
function merge_TS_current_eddy(moor,varargin)
if nargin <1
    help merge_TS_current_eddy
    return
end

% check for optional arguments
a=strmatch('layout',varargin,'exact');
if a>0
    layout=char(varargin(a+1));
else
    layout='portrait';
end

a=strmatch('procpath',varargin,'exact');
if a>0
    procpath=char(varargin(a+1));
else
    data_report_tools_dir=which('data_report_tools');
    b=strfind(data_report_tools_dir,'/');
    data_report_tools_dir=data_report_tools_dir(1:b(end));
    procpath=[data_report_tools_dir '../../../moor/proc/']; %DR changed to be relative paths now that data_report_tools are on the network 19/2/12
end



a=strmatch('proclvl',varargin,'exact');
if a>0
    proclvlstr0=char(varargin(a+1));
    proclvl   = str2num(proclvlstr0);
else
    proclvl=2;
    proclvlstr0 = num2str(proclvl);
end

a=strmatch('plot_interval',varargin,'exact');
if a>0
    plot_interval=eval(varargin{a+1});
else
    plot_interval=0;
end

a=strmatch('unfiltered',varargin,'exact');
if a>0
    unfilt=1;
    proclvlstr = [proclvlstr0 '_unfilt'];    
else
    unfilt=0;
     proclvlstr = [proclvlstr0 '_lpfilt'];       
end


if isunix
    infofile=[procpath,moor,'/',moor,'info.dat'];
elseif ispc
    infofile=[procpath,moor,'\',moor,'info.dat'];
end

% Load vectors of mooring information
% id instrument id, sn serial number, z nominal depth of each instrument
% s_t, e_t, s_d, e_d start and end times and dates
% lat lon mooring position, wd corrected water depth (m)
% mr mooring name
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,...
    'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');


% JULIAN Convert Gregorian date to Julian day.
% JD = JULIAN(YY,MM,DD,HH) or JD = JULIAN([YY,MM,DD,HH]) returns the Julian
% day number of the calendar date specified by year, month, day, and decimal
% hour.
% JD = JULIAN(YY,MM,DD) or JD = JULIAN([YY,MM,DD]) if decimal hour is absent,
% it is assumed to be zero.
% Although the formal definition holds that Julian days start and end at
% noon, here Julian days start and end at midnight. In this convention,
% Julian day 2440000 began at 00:00 hours, May 23, 1968.
jd_start = julian([s_d' hms2h([s_t;0]')]);
jd_end   = julian([e_d' hms2h([e_t;0]')]);

disp('z : instrument id : serial number')
for i = 1:length(id)
    disp([z(i),id(i),sn(i)])
end


% ------------------------------------------------
% Determine plot_interval if not input to function
if plot_interval==0
    plot_interval = zeros(2,4);
    plot_interval(1,1) = s_d(1); plot_interval(1,2) = s_d(2)-1; plot_interval(1,3) = 1; plot_interval(1,4) = 0;
    plot_interval(2,1) = e_d(1); plot_interval(2,2) = e_d(2)+1; plot_interval(2,3) = 1; plot_interval(2,4) = 0;
    if plot_interval(1,2)==0
        plot_interval(1,2)=12; plot_interval(1,1)=plot_interval(1,1)-1;
    end
    if plot_interval(2,2)==13
        plot_interval(2,2)=1; plot_interval(2,1)=plot_interval(2,1)+1;
    end
end


% create xtick spacings based on start of months     
check=0;
i=2;
xticks(1,:)=plot_interval(1,:);
while check~=1
    xticks(i,:)=xticks(i-1,:);
    if xticks(i,2)<12
        xticks(i,2)=xticks(i-1,2)+1;
    else
        xticks(i,2)=1;
        xticks(i,1)=xticks(i-1,1)+1;
    end
    if xticks(i,:)>=plot_interval(2,:)
        check = 1;
    end
    i=i+1;
end

if i<4
   jdxticks =julian(plot_interval(1,:)):(julian(plot_interval(2,:))-julian(plot_interval(1,:)))/5:julian(plot_interval(2,:));
   gxticks = gregorian(jdxticks);
   xticks = gxticks(:,1:4);
   xticklabels = datestr(gxticks,'dd mmm');
else
	jdxticks=julian(xticks);
	% create xticklabels from xticks
	months=['Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'];
	xticklabels=months(xticks(:,2),1:3);
end


% cannot have multi-line xticklabels so have to use manual label command
% this is not really a problem as only want to display years on bottom plot
year_indexes =[];
for i=1:length(xticklabels)
    if find(strfind(xticklabels(i,1:3),'Jan'))
        year_indexes=[year_indexes; i];
    end
end
% use year_indexes later for plotting on bottom graph

jd1 = julian(plot_interval(1,:));
jd2 = julian(plot_interval(2,:)); 

num_to_plot=2; % num_to_plot is the number of samples per day to plot and can be adjusted accordingly

% find the index number of Microcats
iiMC = find(id == 337);
vecMC = sn(iiMC);
% find the index number of RBRs
iiRBR = find(id == 330);
vecRBR = sn(iiRBR);
% find the index number of Idronauts
iiIDR = find(id == 339);
vecIDR = sn(iiIDR);
% find the index number of S4s
iiS4 = find(id == 302);
vecS4 = sn(iiS4);
% and find index number of RCM11s
iiRCM11 = find(id == 310);
vecRCM11 = sn(iiRCM11);
% and find index number of Sontek Argonauts
iiARG = find(id == 366);
vecARG = sn(iiARG);
% and find index number of Nortek Aquadopps
iiNOR = find((id == 368|id==370));
vecNOR = sn(iiNOR);
% and find index number of AADI Seaguards
iiSG = find(id == 301);
vecSG = sn(iiSG);

depths(:,1) = id([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG]);
depths(:,2) = z([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG]);
depths=sortrows(depths,2);
iiiMC=find(depths(:,1)==337);
iiiRBR=find(depths(:,1)==330);
iiiIDR=find(depths(:,1)==339);
iiiS4=find(depths(:,1)==302);  
iiiRCM11=find(depths(:,1)==310);  
iiiARG=find(depths(:,1)==366 | depths(:,1)==366337); 
iiiNOR=find(depths(:,1)==368|depths(:,1)==370);
iiiSG=find(depths(:,1)==301);
iii=[iiiS4;iiiRCM11;iiiARG;iiiNOR;iiiMC;iiiRBR;iiiIDR;iiiSG];

%set figure size on screen for better viewing
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

%set print area of figure
pos1  = [1/8*scnsize(3),8*bdwidth,1/2*scnsize(3),(scnsize(4) - 30*bdwidth)];
pressure_plot=figure('Position',pos1);
set(pressure_plot,'PaperUnits','centimeters');
set(pressure_plot, 'PaperType', 'A4');
set(pressure_plot, 'PaperOrientation',layout);
papersize = get(pressure_plot,'PaperSize');
width=17; height=26; left = (papersize(1)- width)/2; bottom = (papersize(2)- height)/2;
figuresize = [left, bottom, width, height];
set(pressure_plot, 'PaperPosition', figuresize);

plot_string={};

% -----------------------------------
% START OF READING IN INSTRUMENT DATA
% -----------------------------------
%--------------------------------------
% Now read in Microcat data if required
%--------------------------------------
if iiMC>0
    j=1;
    
    % loop to read one file at a time

    for i=1:length(vecMC);
       serialno = vecMC(i);
       disp('*************************************************************')
       disp(['Reading MICROCAT - ',num2str(serialno)])
       disp('*************************************************************')
	if proclvl==2
	       if isunix
        	   infile = [procpath,moor,'/microcat/',moor,'_',sprintf('%0.4d',vecMC(i)),'.use'];
       		elseif ispc
        	   infile = [procpath,moor,'\microcat\',moor,'_',sprintf('%0.4d',vecMC(i)),'.use'];
       		end
	elseif proclvl==3
	       if isunix
        	   infile = [procpath,moor,'/microcat/',moor,'_',sprintf('%0.3d',i),'.microcat'];
       		elseif ispc
        	   infile = [procpath,moor,'\microcat\',moor,'_',sprintf('%0.3d',i),'.microcat'];
       		end	
	end
	
       % check if file exists
       fileopen=fopen(infile,'r');
       
       if fileopen>0
           % read data into vectors and then into structure array

           [yy,mm,dd,hh,t,c,p] = ...
               rodbload(infile,'yy:mm:dd:hh:t:c:p');
           jd=julian(yy,mm,dd,hh);

           bad_data=find(p==-9999); p(bad_data)=NaN;

           eval_string(iiiMC(j))={['MC_' num2str(serialno)]};

           eval([char(eval_string(iiiMC(j))) '.jd=jd;']);
           eval([char(eval_string(iiiMC(j))) '.p=p;']);

           if unfilt==0
           % Apply a butterworth filter to the data using auto_filt and use for
           % plots
           sampling_rate = 1/median(diff(jd));

            ii = eval(['find(~isnan(' char(eval_string(iiiMC(j))) '.p));']); 
            eval([char(eval_string(iiiMC(j))) '.p(ii)=auto_filt(' char(eval_string(iiiMC(j)))...
                  '.p(ii), sampling_rate, 1/2,''low'',4);']);
           end
        end
       j=j+1;
    end
end
%--------------------------------------
% Now read in Aquadopp data if required
%--------------------------------------
if iiNOR>0
    j=1;
    
    % loop to read one file at a time

    for i=1:length(vecNOR);
       serialno = vecNOR(i);
       disp('*************************************************************')
       disp(['Reading AQUADOPP - ',num2str(serialno)])
       disp('*************************************************************')

	if proclvl==2
	       if isunix
        	   infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.use'];
       		elseif ispc
        	   infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.use'];
       		end
	elseif proclvl==3
	       if isunix
        	   infile = [procpath,moor,'/nor/',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
       		elseif ispc
        	   infile = [procpath,moor,'\nor\',moor,'_',sprintf('%3.3d',vecNOR(i)),'.edt'];
       		end	
	end       
       
       % check if file exists
       fileopen=fopen(infile,'r');
      
       if fileopen>0
           % read data into vectors and then into structure array

           [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = ...
               rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd');
           jd=julian(yy,mm,dd,hh);

           bad_data=find(p==-9999); p(bad_data)=NaN;

           eval_string(iiiNOR(j))={['NOR_' num2str(serialno)]};

           eval([char(eval_string(iiiNOR(j))) '.jd=jd;']);
           eval([char(eval_string(iiiNOR(j))) '.p=p;']);
           if unfilt==0
               % Apply a butterworth filter to the data using auto_filt and use for
               % plots
               sampling_rate = 1/median(diff(jd));

                ii = eval(['find(~isnan(' char(eval_string(iiiNOR(j))) '.p));']); 
                eval([char(eval_string(iiiNOR(j))) '.p(ii)=auto_filt(' char(eval_string(iiiNOR(j)))...
                      '.p(ii), sampling_rate, 1/2,''low'',4);']);
           end
       end
       
       j=j+1;
    end
end



%keyboard
print('-dpng',[moor '_merge_TS_current_proclvl_' proclvlstr])
savefig([moor '_merge_TS_current_proclvl_' proclvlstr])
