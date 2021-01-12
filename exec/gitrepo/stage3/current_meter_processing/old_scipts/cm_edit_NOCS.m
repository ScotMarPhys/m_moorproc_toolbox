
function cm_edit_NOCS(moor,varargin)
%
% Function for plotting directions of a mooring overlayed on the same axes
% and speeds of a mooring overlayed on a second pair of axes.
%
% function cm_edit_NOCS('moor','filfac','procpath')
%
% required inputs:-
%   moor: complete mooring name as string. e.g. 'wb1_1_200420'
%
% optional inputs:-
%   filfac: filter factor used for spike removal (default = 10)
%   procpath: can specify exact path to proc directory if not using 
%             standard data paths. 
%             e.g. '/Volumes/noc/mpoc/rpdmoc/rapid/data/moor/proc/'
% 
% Optional process to manually remove spikes with Data Brushing graphic tool
% that is part of Matlab.
%
% functions called:-
%   rodbload, julian
%   from .../exec/moor/tools and .../exec/moor/rodb paths
%   mfilter, uvrot, replace_spike and manual_spike_replace - entered as subfunctions of this script
% 

% Routine written by Darren Rayner August 2009.
% created from Jon Molina's jm_acm_edit_UK_v2008.m script for processing
% RCM11s.
% Aim to unify for all current meters into one function and run from whole
% mooring instead of for each instrument individually.
% changed process slightly so that replaces spikes with NaNs then low-pass
% filters to obtain likely value and then replaces the NaNs with this new
% value.
% Includes request for speed of sound correctinos either using fixed value
% or for RCM11s possibly using as measured conductivity if record looks
% reasonable.


if nargin <1
    help cm_edit_NOCS
    return
end

varargin_string=varargin;
for i=1:length(varargin) % need to change numeric values in varargin to string so can search with strmatch below
    if isnumeric(varargin{i})
        varargin_string{i}=num2str(varargin{i});
    end
end

% check for optional arguments
a=strmatch('filfac',varargin_string,'exact'); % check for filfac first
if a>0
    filfac=str2double(varargin{a+1});
else
    filfac=10; % --- filter factor for despiking
end
a=strmatch('procpath',varargin_string,'exact');
if a>0
    procpath=char(varargin{a+1});
else
    procpath='/noc/mpoc/rpdmoc/rapid/data/moor/proc/'
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
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr,magdec]  =  rodbload(infofile,...
    'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring:MagDeviation');

% JULIAN Convert Gregorian date to Julian day.
% JD = JULIAN(YY,MM,DD,HH) or JD = JULIAN([YY,MM,DD,HH]) returns the Julian
% day number of the calendar date specified by year, month, day, and decimal
% hour.
% JD = JULIAN(YY,MM,DD) or JD = JULIAN([YY,MM,DD]) if decimal hour is absent,
% it is assumed to be zero.
% Although the formal definition holds that Julian days start and end at
% noon, here Julian days start and end at midnight. In this convention,
% Julian day 2440000 began at 00:00 hours, May 23, 1968.

%jd_start = julian([s_d' hms2h([s_t;0]')']);
%jd_end   = julian([e_d' hms2h([e_t;0]')']);

disp('z : instrument id : serial number')
for i = 1:length(id)
    disp([z(i),id(i),sn(i)])
end


% find the index number of S4s
iiS4 = find(id == 302);
vecS4 = sn(iiS4);
% and find index number of RCM11s
iiRCM11 = find(id == 310);
vecRCM11 = sn(iiRCM11);
% and find index number of Sontek Argonauts
iiARG = find(id == 366);
vecARG = sn(iiARG);
% and find index number of Norteks
iiNOR = find(id == 368);
vecNOR = sn(iiNOR);

depths(:,1) = id([iiS4;iiRCM11;iiARG;iiNOR]);
depths(:,2) = z([iiS4;iiRCM11;iiARG;iiNOR]);
depths=sortrows(depths,2);
iiiS4=find(depths(:,1)==302);  
iiiRCM11=find(depths(:,1)==310);  
iiiARG=find(depths(:,1)==366 | depths(:,1)==366337);  
iiiNOR=find(depths(:,1)==368);


%%
%%%%% Determine xticks
% ------------------------------------------------
plot_interval = zeros(2,4);
plot_interval(1,1) = s_d(1); plot_interval(1,2) = s_d(2)-1; plot_interval(1,3) = 1; plot_interval(1,4) = 0;
plot_interval(2,1) = e_d(1); plot_interval(2,2) = e_d(2)+1; plot_interval(2,3) = 1; plot_interval(2,4) = 0;
if plot_interval(1,2)==0
    plot_interval(1,2)=12; plot_interval(1,1)=plot_interval(1,1)-1;
end
if plot_interval(2,2)==13
    plot_interval(2,2)=1; plot_interval(2,1)=plot_interval(2,1)+1;
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
    if xticks(i,:)==plot_interval(2,:)
        check = 1;
        plot_x_labels=1; % toggle value to be used later
    elseif i<3
        if (plot_interval(i,1)==plot_interval(i-1,1)) && (plot_interval(i,2)==plot_interval(i-1,2))
            check = 1;
            disp('short plot_interval so xticks may be limited')
            plot_x_labels=0; % toggle value to be used later
        end
    end
    i=i+1;
end

jdxticks=julian(xticks);

% create xticklabels from xticks
months=['Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'];
xticklabels=months(xticks(:,2),1:3);

if plot_x_labels>0

    % cannot have multi-line xticklabels so have to use manual label command
    % this is not really a problem as only want to display years on bottom plot
    year_indexes =[];
    for i=1:length(xticklabels)
        if find(strfind(xticklabels(i,1:3),'Jan'))
            year_indexes=[year_indexes; i];
        end
    end
    % use year_indexes later for plotting on bottom graph
end

jd1 = julian(plot_interval(1,:));
jd2 = julian(plot_interval(2,:)); 
% end of determining xticks
%%
%setup figures
%set figure size on screen for better viewing
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

%set print area of figure
pos1  = [0.1*scnsize(3),0.5*scnsize(4),0.4*scnsize(3),0.4*(scnsize(4))];
pos2  = [0.5*scnsize(3),0.5*scnsize(4),0.4*scnsize(3),0.4*(scnsize(4))];
pos3  = [0.1*scnsize(3),0.05*scnsize(4),0.4*scnsize(3),0.4*(scnsize(4))];

spikes_u=figure('Position',pos1);
spikes_v=figure('Position',pos2);
low_pass_fig=figure('Position',pos3);




%%
% --------------------
% Read in S4 data if required.
% --------------------
if iiS4>0 
    
    % loop to read one file at a time
    j=1;
    for i=1:length(vecS4);
       serialno = vecS4(i);
       disp('*************************************************************')
       disp(['Reading S4 - ',num2str(serialno)])
       disp('*************************************************************')

       if isunix
           infile = [procpath,moor,'/s4/',moor,'_',sprintf('%4.4d',vecS4(i)),'.use'];
           outfile = ['./',moor,'_',sprintf('%4.4d',vecS4(i)),'.edt'];
       elseif ispc
           infile = [procpath,moor,'\s4\',moor,'_',sprintf('%4.4d',vecS4(i)),'.use'];
           outfile = ['.\',moor,'_',sprintf('%4.4d',vecS4(i)),'.edt'];
       end

       % read data into vectors and then into structure array
       [yy,mm,dd,hh,u,v,t,c,p,hdg] = rodbload(infile,'yy:mm:dd:hh:u:v:t:c:p:hdg');
       jd=julian(yy,mm,dd,hh);

       u(u==-9999)=NaN;
       v(v==-9999)=NaN;
       t(t==-9999)=NaN;
       c(c==-9999)=NaN;
       p(p==-9999)=NaN;
       hdg(hdg==-9999)=NaN;
       
       % CORRECT FOR MAGNETIC DEVIATION GIVEN BY MAGDEC
       %use uvrot by Visbeck
       [u,v]=uvrot(u,v,magdec);

       sampling_interval = 24*median(diff(jd)); % in hours
       nt=round(40/sampling_interval); % number of samples in 40 hours - used for 40-hour lowpass filtering
       
       % high pass data for spike identification
       uhf=sqrt(mfilter(u(~isnan(u)),1,0,1/nt).^2);
       vhf=sqrt(mfilter(v(~isnan(v)),1,0,1/nt).^2);
       ii=find(uhf>filfac*mean(uhf) | vhf>filfac*mean(vhf)); % line that identifies spikes
       
       figure(spikes_u)
       plot(jd-jd(1),uhf,jd(ii)-jd(1),uhf(ii),'*r'); 
       
       title(['S4 SN: ',sprintf('%4.4d',vecS4(i)), ' high pass filtered u-component'])
       % end of automatic spike identification for u
       
       % replace spikes with NaNs
       u(ii)=NaN;
       uhf(ii)=NaN;
       
       figure(spikes_v)
       plot(jd-jd(1),vhf,jd(ii)-jd(1),vhf(ii),'or');
       
       title(['S4 SN: ',sprintf('%4.4d',vecS4(i)), ' high pass filtered v-component'])
       % end of automatic spike identification for v
       
       % replace spikes with NaNs
       v(ii)=NaN;
       vhf(ii)=NaN;
       
       figure(spikes_u)
       % manually replace spikes and bad data in u if required
       disp('u-component despiking:')
       u=manual_spike_replace(u,jd,uhf,spikes_u);
       % replace corresponding values in other timeseries
       a=find(isnan(u));
       v(a)=NaN; p(a)=NaN; t(a)=NaN; c(a)=NaN; uhf(a)=NaN; vhf(a)=NaN;
       % replot to show edited high frequency series
       plot(jd-jd(1),uhf);
       figure(spikes_v)
       plot(jd-jd(1),vhf);
       
       figure(spikes_v)
       disp('v-component despiking:')
       % manually replace spikes and bad data in v if required
       v=manual_spike_replace(v,jd,vhf,spikes_v);
       % replace corresponding values in other timeseries
       a=find(isnan(v));
       u(a)=NaN; p(a)=NaN; t(a)=NaN; c(a)=NaN; uhf(a)=NaN; vhf(a)=NaN;
       % replot to show edited high frequency series
       plot(jd-jd(1),vhf);
       figure(spikes_u)
       plot(jd-jd(1),uhf);
       
       % interpolate over NaNs to fill in gaps in time series if less than
       % 24 hour gaps
       % this replace_spike function replaces Jon Molina's code for spike 
       % replacement which just used low pass filtering which can be skewed
       % if have very large spikes or spikes near either end of the series
       u=replace_spike(u,24/sampling_interval);
       v=replace_spike(v,24/sampling_interval);
       p=replace_spike(p,24/sampling_interval);
       t=replace_spike(t,24/sampling_interval);
       c=replace_spike(c,24/sampling_interval);
       
%        % Jon Molina's original method for spike replacement
%        % determine likely value at spike position first iteration
%        uf=mfilter(u,1,1/nt,0);
%        vf=mfilter(v,1,1/nt,0);
%        pf=mfilter(p,1,1/nt,0);
%        tf=mfilter(t,1,1/nt,0);
%        cf=mfilter(c,1,1/nt,0);
% 
%        % replace spike value
%        u(ii)=uf(ii); v(ii)=vf(ii); p(ii)=pf(ii); t(ii)=tf(ii); c(ii)=cf(ii);
%        %despiked

       % NEED TO REMOVE NANS BEFORE FILTERING. OTHERWISE REPLACES WHOLE SERIES WITH NaNs.
       abc=find(~isnan(u));

       % low pass product again
       uf=mfilter(u(abc),1,1/nt,0);
       vf=mfilter(v(abc),1,1/nt,0);
       pf=mfilter(p(abc),1,1/nt,0);
       tf=mfilter(t(abc),1,1/nt,0);
       cf=mfilter(c(abc),1,1/nt,0);

       % 12h resolution
       tim=fix(jd(1)):0.5:fix(jd(end));
       vf1=interp1(jd(abc),vf,tim);
       uf1=interp1(jd(abc),uf,tim);
       pf1=interp1(jd(abc),pf,tim);
       tf1=interp1(jd(abc),tf,tim);
       cf1=interp1(jd(abc),cf,tim);

       figure(low_pass_fig)
       plot(tim-jd1,uf1); hold on;
       plot(tim-jd1,vf1,'-r'); 
       xlim([0 jd2-jd1]);
       set(gca,'xTickLabel',xticklabels);
       set(gca,'XTick',jdxticks-jd1);
       if plot_x_labels>0
            ylimits=get(gca,'ylim');
            y_label_pos=ylimits(1)-abs(ylimits(1)-ylimits(2))/10;
            for k=1:length(year_indexes)
                text((jd2-jd1)*(year_indexes(k)-1)/(length(xticklabels)-1),y_label_pos,num2str(xticks(year_indexes(k),1)),'FontSize',10);
            end
       end
       hold off;
       title(['S4 SN: ',sprintf('%4.4d',vecS4(i)), ' 12-hour resolution low-pass filtered currents'])
       legend({'u-component','v-component'})
       hold off
       

       %replace NANs with -9999
       tf1(isnan(tf1))=-9999;
       pf1(isnan(pf1))=-9999;
       uf1(isnan(uf1))=-9999;
       vf1(isnan(vf1))=-9999;
       cf1(isnan(cf1))=-9999;
       
       % Save to rodb format
       gt=gregorian(tim);
       YY=gt(:,1); MM=gt(:,2); DD=gt(:,3); 
       if size(gt,2) == 6
           HH=hms2h(gt(:,4),gt(:,5),gt(:,6)); 
       else 
           HH= gt(:,4);
       end  
       fort = '%4d  %2d  %2d  %8.4f  %6.3f  %6.1f  %6.4f  %7.2f  %7.2f';
       cols = 'YY:MM:DD:HH:T:P:C:U:V';  % NB: not including heading as not straight forward to filter.
       data=[YY MM DD HH tf1' pf1' cf1' uf1' vf1'];
       infovar=['Mooring:Latitude:Longitude:Columns:Start_Date:Start_Time:End_Date:End_Time:SerialNumber:' ...
           'WaterDepth:InstrDepth:MagDeviation'];
       rodbsave(outfile,infovar,fort,moor,lat,lon,cols,s_d,s_t,e_d,e_t,serialno,wd,depths(iiiS4(i),2),magdec,...
                data);
     
       j=j+1;
       disp('PAUSED - press any key to continue to next instrument')
       pause
       
    end
end
%%
% --------------------
% Read in RCM11 data if required.
% --------------------
if iiRCM11>0 
    
    % loop to read one file at a time
    j=1;
    for i=1:length(vecRCM11);
       serialno = vecRCM11(i);
       disp('*************************************************************')
       disp(['Reading RCM11 - ',num2str(serialno)])
       disp('*************************************************************')

       if isunix
           infile = [procpath,moor,'/rcm/',moor,'_',sprintf('%3.3d',vecRCM11(i)),'.use'];
           outfile = ['./',moor,'_',sprintf('%3.3d',vecRCM11(i)),'.edt'];
       elseif ispc
           infile = [procpath,moor,'\rcm\',moor,'_',sprintf('%3.3d',vecRCM11(i)),'.use'];
           outfile = ['.\',moor,'_',sprintf('%3.3d',vecRCM11(i)),'.edt'];
       end

       % read data into vectors and then into structure array
       [yy,mm,dd,hh,ref,u,v,t,c,p,tlt,mss] = rodbload(infile,'yy:mm:dd:hh:ref:u:v:t:c:p:tlt:mss');
       jd=julian(yy,mm,dd,hh);

       t(t==-9999)=NaN;
       p(p==-9999)=NaN;
       u(u==-9999)=NaN;
       v(v==-9999)=NaN;
       c(c==-9999)=NaN;
       tlt(tlt==-9999)=NaN;
       mss(mss==-9999)=NaN;
       
       %%
       % CORRECT FOR SOUND VELOCITY OF SEAWATER
       % first display plot of conductivity and ask operator if wants to
       % use conductivity as measured or a fixed value.
       % currently automatically uses measured temperature
       figure(999)
       plot(jd-jd1,c)
       xlim([0 jd2-jd1]);
       set(gca,'xTickLabel',xticklabels);
       set(gca,'XTick',jdxticks-jd1);
       if plot_x_labels>0
            ylimits=get(gca,'ylim');
            y_label_pos=ylimits(1)-abs(ylimits(1)-ylimits(2))/10;
            for k=1:length(year_indexes)
                text((jd2-jd1)*(year_indexes(k)-1)/(length(xticklabels)-1),y_label_pos,num2str(xticks(year_indexes(k),1)),'FontSize',10);
            end
       end
       title(['Conductivity as measured by RCM11 serial number ', sprintf('%3.3d',vecRCM11(i))])
       use_c=input(['Do you want to use the as measured conductivity for speed\n' ...
                   'of sound correction (see plot) or a fixed salinity value?\n' ...
                   'for measured enter 0, for fixed, enter value: ']);
       close 999
       t68=t*1.00024; % convert T90 temp values to T68 for use with CSIRO routines
       if use_c==0
            c3515 = sw_c3515;
            c_ratio=c/c3515;
            s=sw_salt(c_ratio,t68,p);
       else
           s(1:length(t),1)=use_c; %creates vector of same length as t filled with use_c value for fixed salinity.
       end
       svel = sw_svel(s,t68,p);
       svel_old=1500; % 1500m/s specific for RCM11s
       corfac= (svel./svel_old);
       us = u .* corfac;
       vs = v .* corfac;
       % end of correcting for speed of sound
       %%
       
       % CORRECT FOR MAGNETIC DEVIATION GIVEN BY MAGDEC
       %use uvrot by Visbeck
       [u,v]=uvrot(us,vs,magdec);

        
       sampling_interval = 24*median(diff(jd)); % in hours
       nt=round(40/sampling_interval); % number of samples in 40 hours - used for 40-hour lowpass filtering
       
       
       % high pass data for spike identification
       % remove NaNs during filter step
       uhf=sqrt(mfilter(u(~isnan(u)),1,0,1/nt).^2);
       vhf=sqrt(mfilter(v(~isnan(v)),1,0,1/nt).^2);
       ii=find(uhf>filfac*mean(uhf) | vhf>filfac*mean(vhf)); % line that identifies spikes - doesn't identify spikes in T and P but uses
                                                             % spikes found in U and V to position spikes for T and P.
       
       figure(spikes_u)
       plot(jd-jd(1),uhf,jd(ii)-jd(1),uhf(ii),'*r'); 
       
       title(['RCM11 SN: ',sprintf('%3.3d',vecRCM11(i)), ' high pass filtered u-component'])
       % replace spikes with NaNs
       u(ii)=NaN;
       uhf(ii)=NaN;
       
       figure(spikes_v)
       plot(jd-jd(1),vhf,jd(ii)-jd(1),vhf(ii),'or');
      
       title(['RCM11 SN: ',sprintf('%3.3d',vecRCM11(i)), ' high pass filtered v-component'])
       % end of spike identification
       % replace spikes with NaNs
       v(ii)=NaN;
       vhf(ii)=NaN;
       
       figure(spikes_u)
       disp('u-component despiking:')
       % manually replace spikes and bad data in u if required
       u=manual_spike_replace(u,jd,uhf,spikes_u);
       % replace corresponding values in other timeseries
       a=find(isnan(u));
       v(a)=NaN; p(a)=NaN; t(a)=NaN; c(a)=NaN; uhf(a)=NaN; vhf(a)=NaN; mss(a)=NaN; tlt(a)=NaN;
       % replot to show edited high frequency series
       plot(jd-jd(1),uhf);
       figure(spikes_v)
       plot(jd-jd(1),vhf);
       
       figure(spikes_v)
       disp('v-component despiking:')
       % manually replace spikes and bad data in v if required
       v=manual_spike_replace(v,jd,vhf,spikes_v);
       % replace corresponding values in other timeseries
       a=find(isnan(v));
       u(a)=NaN; p(a)=NaN; t(a)=NaN; c(a)=NaN; uhf(a)=NaN; vhf(a)=NaN; mss(a)=NaN; tlt(a)=NaN;
       % replot to show edited high frequency series
       plot(jd-jd(1),vhf);
       figure(spikes_u)
       plot(jd-jd(1),uhf);
       
       % interpolate over NaNs to fill in gaps in time series if less than
       % 24 hour gaps
       % this replace_spike function replaces Jon Molina's code for spike 
       % replacement which just used low pass filtering which can be skewed
       % if have very large spikes or spikes near either end of the series
       u=replace_spike(u,24/sampling_interval);
       v=replace_spike(v,24/sampling_interval);
       p=replace_spike(p,24/sampling_interval);
       t=replace_spike(t,24/sampling_interval);
       c=replace_spike(c,24/sampling_interval);
       mss=replace_spike(mss,24/sampling_interval);
       tlt=replace_spike(tlt,24/sampling_interval);
       
%        % Jon Molina's original method for spike removal
%        % determine likely value at spike position first iteration using
%        % lowpass filter incuding NaNs       
%        uf=mfilter(u,1,1/nt,0);
%        vf=mfilter(v,1,1/nt,0);
%        pf=mfilter(p,1,1/nt,0);
%        tf=mfilter(t,1,1/nt,0);
%        cf=mfilter(c,1,1/nt,0);
%        mssf=mfilter(mss,1,1/nt,0);
%        tltf=mfilter(tlt,1,1/nt,0);
% 
%        % replace spike value
%        % NB auxillary variables also replaced with low-pass filtered value 
%        % if identified as spike in U or V varible.
%        u(ii)=uf(ii); v(ii)=vf(ii); p(ii)=pf(ii); t(ii)=tf(ii); c(ii)=cf(ii); mss(ii)=mssf(ii); tlt(ii)=tltf(ii);
%        %despiked


       % NEED TO REMOVE NANS BEFORE FILTERING. OTHERWISE REPLACES WHOLE SERIES WITH NaNs.
       abc=find(~isnan(u));

       % low pass product again
       uf=mfilter(u(abc),1,1/nt,0);
       vf=mfilter(v(abc),1,1/nt,0);
       pf=mfilter(p(abc),1,1/nt,0);
       tf=mfilter(t(abc),1,1/nt,0);
       cf=mfilter(c(abc),1,1/nt,0);
       mssf=mfilter(mss(abc),1,1/nt,0);
       tltf=mfilter(tlt(abc),1,1/nt,0);

       % 12h resolution
       tim=fix(jd(1)):0.5:fix(jd(end));
       vf1=interp1(jd(abc),vf,tim);
       uf1=interp1(jd(abc),uf,tim);
       pf1=interp1(jd(abc),pf,tim);
       tf1=interp1(jd(abc),tf,tim);
       cf1=interp1(jd(abc),cf,tim);
       mssf1=interp1(jd(abc),mssf,tim);
       tltf1=interp1(jd(abc),tltf,tim);

       figure(low_pass_fig)
       plot(tim-jd1,uf1); hold on;
       plot(tim-jd1,vf1,'-r'); 
       xlim([0 jd2-jd1]);
       set(gca,'xTickLabel',xticklabels);
       set(gca,'XTick',jdxticks-jd1);
       if plot_x_labels>0
            ylimits=get(gca,'ylim');
            y_label_pos=ylimits(1)-abs(ylimits(1)-ylimits(2))/10;
            for k=1:length(year_indexes)
                text((jd2-jd1)*(year_indexes(k)-1)/(length(xticklabels)-1),y_label_pos,num2str(xticks(year_indexes(k),1)),'FontSize',10);
            end
       end
       hold off;
       title(['RCM11 SN: ',sprintf('%3.3d',vecRCM11(i)), ' 12-hour resolution low-pass filtered currents'])
       legend({'u-component','v-component'})
       

       %replace NANs with -9999
       tf1(isnan(tf1))=-9999;
       pf1(isnan(pf1))=-9999;
       uf1(isnan(uf1))=-9999;
       vf1(isnan(vf1))=-9999;
       cf1(isnan(cf1))=-9999;
       mssf1(isnan(mssf1))=-9999;
       tltf1(isnan(tltf1))=-9999;
       
       % Save to rodb format
       gt=gregorian(tim);
       YY=gt(:,1); MM=gt(:,2); DD=gt(:,3); 
       if size(gt,2) == 6
           HH=hms2h(gt(:,4),gt(:,5),gt(:,6)); 
       else 
           HH= gt(:,4);
       end  
       fort = '%4d  %2d  %2d  %8.4f  %6.4f  %6.1f  %6.4f  %7.2f  %7.2f  %2.2f  %2.2f';
       cols = 'YY:MM:DD:HH:T:P:C:U:V:MSS:TLT';
       data=[YY MM DD HH tf1' pf1' cf1' uf1' vf1' mssf1' tltf1'];
       infovar=['Mooring:Latitude:Longitude:Columns:Start_Date:Start_Time:End_Date:End_Time:SerialNumber:' ...
           'WaterDepth:InstrDepth:MagDeviation'];
       rodbsave(outfile,infovar,fort,moor,lat,lon,cols,s_d,s_t,e_d,e_t,serialno,wd,depths(iiiRCM11(i),2),magdec,...
                data);
       
    
       j=j+1;
       disp('PAUSED - press any key to continue to next instrument')
       pause
    end
end



%%
% --------------------
% Read in NORTEK data if required.
% --------------------
if iiNOR>0 
    
    % loop to read one file at a time
    j=1;
    for i=1:length(vecNOR);
       serialno = vecNOR(i);
       disp('*************************************************************')
       disp(['Reading NORTEK - ',num2str(serialno)])
       disp('*************************************************************')

       if isunix
           infile = [procpath,moor,'/nor/',moor,'_',sprintf('%4.4d',vecNOR(i)),'.use'];
           outfile = ['./',moor,'_',sprintf('%4.4d',vecNOR(i)),'.edt'];
       elseif ispc
           infile = [procpath,moor,'\nor\',moor,'_',sprintf('%4.4d',vecNOR(i)),'.use'];
           outfile = ['.\',moor,'_',sprintf('%4.4d',vecNOR(i)),'.edt'];
       end

       % read data into vectors and then into structure array
       [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd');
       jd=julian(yy,mm,dd,hh);

       t(t==-9999)=NaN;
       p(p==-9999)=NaN;
       u(u==-9999)=NaN;
       v(v==-9999)=NaN;
       w(w==-9999)=NaN;
       hdg(hdg==-9999)=NaN;
       pit(pit==-9999)=NaN;
       rol(rol==-9999)=NaN;
       uss(uss==-9999)=NaN;
       vss(vss==-9999)=NaN;
       wss(wss==-9999)=NaN;
       ipow(ipow==-9999)=NaN;
       % don't bother with current speed (cs) and current direction (cd) as
       % can be calculated from u an v.
       
       %%
       % CORRECT FOR SOUND VELOCITY OF SEAWATER
       % no option to use as measured c for a Nortek
       % currently automatically uses measured temperature
       % ask for operator's input of fixed salinity value
       s_old = input('What was the fixed value of salinity used when the instrument was setup? ');
       use_c=input('What value to use for fixed salinity (in PSU) for speed of sound correction? ');
       
       t68=t*1.00024; % convert T90 temp values to T68 for use with CSIRO routines
       s(1:length(t),1)=use_c; %creates vector of same length as t filled with use_c value for fixed salinity.
       svel = sw_svel(s,t68,p);
       %Norteks and Sonteks usually set up to use measured temperature and
       %fixed salinity value of 35 but need to check for each file.
       s_old=s*0+s_old; % adjust array to be same length as t68 and p
       svel_old=sw_svel(s_old,t68,p); 
       corfac= (svel./svel_old);
       us = u .* corfac;
       vs = v .* corfac;
       % end of correcting for speed of sound
       %%
       
       % CORRECT FOR MAGNETIC DEVIATION GIVEN BY MAGDEC
       %use uvrot by Visbeck
       [u,v]=uvrot(us,vs,magdec);

       sampling_interval = 24*median(diff(jd)); % in hours
       nt=round(40/sampling_interval); % number of samples in 40 hours - used for 40-hour lowpass filtering
       
       % high pass data for spike identification
       uhf=sqrt(mfilter(u(~isnan(u)),1,0,1/nt).^2);
       vhf=sqrt(mfilter(v(~isnan(v)),1,0,1/nt).^2);
       ii=find(uhf>filfac*mean(uhf) | vhf>filfac*mean(vhf)); % line that identifies spikes - doesn't identify spikes in T and P but uses
                                                             % spikes found in U and V to position spikes for T and P.
       
       figure(spikes_u)
       plot(jd-jd(1),uhf,jd(ii)-jd(1),uhf(ii),'*r'); 
       
       title(['NORTEK SN: ',sprintf('%4.4d',vecNOR(i)), ' high pass filtered u-component'])
       % end of automatic spike identification for u
       
       % replace spikes with NaNs
       u(ii)=NaN;
       uhf(ii)=NaN;
       
       figure(spikes_v)
       plot(jd-jd(1),vhf,jd(ii)-jd(1),vhf(ii),'or');
       
       title(['NORTEK SN: ',sprintf('%4.4d',vecNOR(i)), ' high pass filtered v-component'])
       % end of spike identification
       % replace spikes with NaNs
       v(ii)=NaN;
       vhf(ii)=NaN;
       
       figure(spikes_u)
       disp('u-component despiking:')
       % manually replace spikes and bad data in u if required
       u=manual_spike_replace(u,jd,uhf,spikes_u);
       % replace corresponding values in other timeseries
       a=find(isnan(u));
       v(a)=NaN; p(a)=NaN; t(a)=NaN; w(a)=NaN; uhf(a)=NaN; vhf(a)=NaN; pit(a)=NaN;
       uss(a)=NaN; vss(a)=NaN; wss(a)=NaN; ipow(a)=NaN;
       % replot to show edited high frequency series
       plot(jd-jd(1),uhf);
       figure(spikes_v)
       plot(jd-jd(1),vhf);
       
       figure(spikes_v)
       disp('v-component despiking:')
       % manually replace spikes and bad data in v if required
       v=manual_spike_replace(v,jd,vhf,spikes_v);
       % replace corresponding values in other timeseries
       a=find(isnan(v));
       u(a)=NaN; p(a)=NaN; t(a)=NaN; c(a)=NaN; uhf(a)=NaN; vhf(a)=NaN; pit(a)=NaN;
       uss(a)=NaN; vss(a)=NaN; wss(a)=NaN; ipow(a)=NaN;
       % replot to show edited high frequency series
       plot(jd-jd(1),vhf);
       figure(spikes_u)
       plot(jd-jd(1),uhf);
       
       % interpolate over NaNs to fill in gaps in time series if less than
       % 24 hour gaps
       % this replace_spike function replaces Jon Molina's code for spike 
       % replacement which just used low pass filtering which can be skewed
       % if have very large spikes or spikes near either end of the series
       u=replace_spike(u,24/sampling_interval);
       v=replace_spike(v,24/sampling_interval);
       p=replace_spike(p,24/sampling_interval);
       t=replace_spike(t,24/sampling_interval);
       w=replace_spike(w,24/sampling_interval);
       pit=replace_spike(pit,24/sampling_interval);
       rol=replace_spike(rol,24/sampling_interval);
       uss=replace_spike(uss,24/sampling_interval);
       vss=replace_spike(vss,24/sampling_interval);
       wss=replace_spike(wss,24/sampling_interval);
       ipow=replace_spike(ipow,24/sampling_interval);
       
%        % Jon Molina's original method for spike removal
%        % determine likely value at spike position first iteration using
%        % lowpass filter incuding NaNs
%        uf=mfilter(u,1,1/nt,0);
%        vf=mfilter(v,1,1/nt,0);
%        pf=mfilter(p,1,1/nt,0);
%        tf=mfilter(t,1,1/nt,0);
%        wf=mfilter(w,1,1/nt,0);
%        pitf=mfilter(pit,1,1/nt,0);
%        rolf=mfilter(rol,1,1/nt,0);
%        ussf=mfilter(uss,1,1/nt,0);
%        vssf=mfilter(vss,1,1/nt,0);
%        wssf=mfilter(wss,1,1/nt,0);
%        ipowf=mfilter(ipow,1,1/nt,0);
%        
%        % replace spike value
%        % NB auxillary variables also replaced with low-pass filtered value 
%        % if identified as spike in U or V varible.
%        u(ii)=uf(ii); v(ii)=vf(ii); p(ii)=pf(ii); t(ii)=tf(ii); w(ii)=wf(ii); pit(ii)=pitf(ii); rol(ii)=rolf(ii);
%        uss(ii)=ussf(ii); vss(ii)=vssf(ii); wss(ii)=wssf(ii); ipow(ii)=ipowf(ii);
%        %despiked

       % NEED TO REMOVE NANS BEFORE FILTERING. OTHERWISE REPLACES WHOLE SERIES WITH NaNs.
       abc=find(~isnan(u));

       % low pass product again
       uf=mfilter(u(abc),1,1/nt,0);
       vf=mfilter(v(abc),1,1/nt,0);
       pf=mfilter(p(abc),1,1/nt,0);
       tf=mfilter(t(abc),1,1/nt,0);
       wf=mfilter(w(abc),1,1/nt,0);
       ussf=mfilter(uss(abc),1,1/nt,0);
       vssf=mfilter(vss(abc),1,1/nt,0);
       wssf=mfilter(wss(abc),1,1/nt,0);
       pitf=mfilter(pit(abc),1,1/nt,0);
       rolf=mfilter(rol(abc),1,1/nt,0);
       ipowf=mfilter(ipow(abc),1,1/nt,0);
       
       % 12h resolution
       tim=fix(jd(1)):0.5:fix(jd(end));
       vf1=interp1(jd(abc),vf,tim);
       uf1=interp1(jd(abc),uf,tim);
       pf1=interp1(jd(abc),pf,tim);
       tf1=interp1(jd(abc),tf,tim);
       wf1=interp1(jd(abc),wf,tim);
       ussf1=interp1(jd(abc),ussf,tim);
       vssf1=interp1(jd(abc),vssf,tim);
       wssf1=interp1(jd(abc),wssf,tim);
       pitf1=interp1(jd(abc),pitf,tim);
       rolf1=interp1(jd(abc),rolf,tim);
       ipowf1=interp1(jd(abc),ipowf,tim);


       figure(low_pass_fig)
       plot(tim-jd1,uf1); hold on;
       plot(tim-jd1,vf1,'-r');
       plot(tim-jd1,wf1,'-g');
       xlim([0 jd2-jd1]);
       set(gca,'xTickLabel',xticklabels);
       set(gca,'XTick',jdxticks-jd1);
       if plot_x_labels>0
            ylimits=get(gca,'ylim');
            y_label_pos=ylimits(1)-abs(ylimits(1)-ylimits(2))/10;
            for k=1:length(year_indexes)
                text((jd2-jd1)*(year_indexes(k)-1)/(length(xticklabels)-1),y_label_pos,num2str(xticks(year_indexes(k),1)),'FontSize',10);
            end
       end
       hold off;
       title(['NORTEK SN: ',sprintf('%4.4d',vecNOR(i)), ' 12-hour resolution low-pass filtered currents'])
       legend({'u-component','v-component','w-component'})
       

       %replace NANs with -9999
       tf1(isnan(tf1))=-9999;
       pf1(isnan(pf1))=-9999;
       uf1(isnan(uf1))=-9999;
       vf1(isnan(vf1))=-9999;
       wf1(isnan(wf1))=-9999;
       ussf1(isnan(ussf1))=-9999;
       vssf1(isnan(vssf1))=-9999;
       wssf1(isnan(wssf1))=-9999;
       pitf1(isnan(pitf1))=-9999;
       rolf1(isnan(rolf1))=-9999;
       ipowf1(isnan(ipowf1))=-9999;
       
       % Save to rodb format
       gt=gregorian(tim);
       YY=gt(:,1); MM=gt(:,2); DD=gt(:,3); 
       if size(gt,2) == 6
           HH=hms2h(gt(:,4),gt(:,5),gt(:,6)); 
       else 
           HH= gt(:,4);
       end  
       fort = '%4d  %2d  %2d  %8.4f  %6.4f  %6.1f  %7.2f  %7.2f  %7.2f  %2.0d  %2.0d  %2.0d  %3.1d  %3.1d  %3.1d';
       cols = 'YY:MM:DD:HH:T:P:U:V:W:USS:VSS:WSS:PIT:ROL:IPOW';
       data=[YY MM DD HH tf1' pf1' uf1' vf1' wf' ussf1' vssf1' wssf1' pitf1' rolf1' ipowf1'];
       infovar=['Mooring:Latitude:Longitude:Columns:Start_Date:Start_Time:End_Date:End_Time:SerialNumber:' ...
           'WaterDepth:InstrDepth:MagDeviation'];
       rodbsave(outfile,infovar,fort,moor,lat,lon,cols,s_d,s_t,e_d,e_t,serialno,wd,depths(iiiNOR(i),2),magdec,...
                data);
    
       j=j+1;
       disp('PAUSED - press any key to continue to next instrument')
       pause
       
    end
end



%%
% --------------------
% Read in SONTEK data if required.
% --------------------
if iiARG>0 
    
    % loop to read one file at a time
    j=1;
    for i=1:length(vecARG);
       serialno = vecARG(i);
       disp('*************************************************************')
       disp(['Reading SONTEK - ',num2str(serialno)])
       disp('*************************************************************')

       if isunix
           infile = [procpath,moor,'/arg/',moor,'_',sprintf('%3.3d',vecARG(i)),'.use'];
           outfile = ['./',moor,'_',sprintf('%3.3d',vecARG(i)),'.edt'];
           if exist(infile,'file')==0  % older Arg files had 4 digit serial number starting with zero in filename
               infile = [procpath,moor,'/arg/',moor,'_0',sprintf('%3.3d',vecARG(i)),'.use'];
           end
       elseif ispc
           infile = [procpath,moor,'\arg\',moor,'_',sprintf('%3.3d',vecARG(i)),'.use'];
           outfile = ['.\',moor,'_',sprtinf('%3.3d',vecARG(i)),'.edt'];
           if exist(infile,'file')==0  % older Arg files had 4 digit serial number starting with zero in filename
               infile = [procpath,moor,'\arg\',moor,'_0',sprtinf('%3.3d',vecARG(i)),'.use'];
           end
       end

       % read data into vectors and then into structure array
       [yy,mm,dd,hh,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow] = ...
           rodbload(infile,'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow');
       jd=julian(yy,mm,dd,hh);
       
       t(t==-9999)=NaN;
       p(p==-9999)=NaN;
       u(u==-9999)=NaN;
       v(v==-9999)=NaN;
       w(w==-9999)=NaN;
       pit(pit==-9999)=NaN;
       rol(rol==-9999)=NaN;
       uss(uss==-9999)=NaN;
       vss(vss==-9999)=NaN;
       wss(wss==-9999)=NaN;
       ipow(ipow==-9999)=NaN;
       
       %%
       % CORRECT FOR SOUND VELOCITY OF SEAWATER
       % no option to use as measured c for a Sontek
       % currently automatically uses measured temperature
       % ask for operator's input of fixed salinity value
       s_old = input('What was the fixed value of salinity used when the instrument was setup? ');
       use_c=input('What value to use for fixed salinity (in PSU) for speed of sound correction? ');
       
       t68=t*1.00024; % convert T90 temp values to T68 for use with CSIRO routines
       s(1:length(t),1)=use_c; %creates vector of same length as t filled with use_c value for fixed salinity.
       svel = sw_svel(s,t68,p);
       %Norteks and Sonteks usually set up to use measured temperature and
       %fixed salinity value of 35 but need to check for each file.
       s_old=s*0+s_old; % adjust array to be same length as t68 and p
       svel_old=sw_svel(s_old,t68,p); 
       corfac= (svel./svel_old);
       us = u .* corfac;
       vs = v .* corfac;
       % end of correcting for speed of sound
       %%

       % CORRECT FOR MAGNETIC DEVIATION GIVEN BY MAGDEC
       %use uvrot by Visbeck
       [u,v]=uvrot(us,vs,magdec);

       sampling_interval = 24*median(diff(jd)); % in hours
       nt=round(40/sampling_interval); % number of samples in 40 hours - used for 40-hour lowpass filtering
       
       % high pass data for spike identification
       uhf=sqrt(mfilter(u(~isnan(u)),1,0,1/nt).^2);
       vhf=sqrt(mfilter(v(~isnan(v)),1,0,1/nt).^2);
       ii=find(uhf>filfac*mean(uhf) | vhf>filfac*mean(vhf)); % line that identifies spikes - doesn't identify spikes in T and P but uses
                                                             % spikes found in U and V to position spikes for T and P.
       
       figure(spikes_u)
       plot(jd-jd(1),uhf,jd(ii)-jd(1),uhf(ii),'*r'); 
       
       title(['SONTEK SN: ',sprintf('%3.3d',vecARG(i)), ' high pass filtered u-component'])
       % end of automatic spike identification for u
       
       % replace spikes with NaNs
       u(ii)=NaN;
       uhf(ii)=NaN;
       
       figure(spikes_v)
       plot(jd-jd(1),vhf,jd(ii)-jd(1),vhf(ii),'or');
       
       title(['SONTEK SN: ',sprintf('%3.3d',vecARG(i)), ' high pass filtered v-component'])
       % end of automatic spike identification for v
       
       % replace spikes with NaNs
       v(ii)=NaN;
       vhf(ii)=NaN;
       
       figure(spikes_u)
       disp('u-component despiking:')
       % manually replace spikes and bad data in u if required
       u=manual_spike_replace(u,jd,uhf,spikes_u);
       % replace corresponding values in other timeseries
       a=find(isnan(u));
       v(a)=NaN; p(a)=NaN; t(a)=NaN; w(a)=NaN; uhf(a)=NaN; vhf(a)=NaN; pit(a)=NaN;
       rol(a)=NaN; uss(a)=NaN; vss(a)=NaN; wss(a)=NaN; ipow(a)=NaN;
       % replot to show edited high frequency series
       plot(jd-jd(1),uhf);
       figure(spikes_v)
       plot(jd-jd(1),vhf);
       
       figure(spikes_v)
       disp('v-component despiking:')
       % manually replace spikes and bad data in v if required
       v=manual_spike_replace(v,jd,vhf,spikes_v);
       % replace corresponding values in other timeseries
       a=find(isnan(v));
       u(a)=NaN; p(a)=NaN; t(a)=NaN; w(a)=NaN; uhf(a)=NaN; vhf(a)=NaN; pit(a)=NaN;
       rol(a)=NaN; uss(a)=NaN; vss(a)=NaN; wss(a)=NaN; ipow(a)=NaN;
       % replot to show edited high frequency series
       plot(jd-jd(1),vhf);
       figure(spikes_u)
       plot(jd-jd(1),uhf);
       
       % interpolate over NaNs to fill in gaps in time series if less than
       % 24 hour gaps
       % this replace_spike function replaces Jon Molina's code for spike 
       % replacement which just used low pass filtering which can be skewed
       % if have very large spikes or spikes near either end of the series
       u=replace_spike(u,24/sampling_interval);
       v=replace_spike(v,24/sampling_interval);
       p=replace_spike(p,24/sampling_interval);
       t=replace_spike(t,24/sampling_interval);
       w=replace_spike(w,24/sampling_interval);
       pit=replace_spike(pit,24/sampling_interval);
       rol=replace_spike(rol,24/sampling_interval);
       ipow=replace_spike(ipow,24/sampling_interval);
       uss=replace_spike(uss,24/sampling_interval);
       vss=replace_spike(vss,24/sampling_interval);
       wss=replace_spike(wss,24/sampling_interval);
      
       
%        Jon Molina's original method for spike replacement
%        % determine likely value at spike position first iteration using
%        % lowpass filter incuding NaNs
%        uf=mfilter(u,1,1/nt,0);
%        vf=mfilter(v,1,1/nt,0);
%        pf=mfilter(p,1,1/nt,0);
%        tf=mfilter(t,1,1/nt,0);
%        wf=mfilter(w,1,1/nt,0);
%        pitf=mfilter(pit,1,1/nt,0);
%        rolf=mfilter(rol,1,1/nt,0);
%        ussf=mfilter(uss,1,1/nt,0);
%        vssf=mfilter(vss,1,1/nt,0);
%        wssf=mfilter(wss,1,1/nt,0);
%        ipowf=mfilter(ipow,1,1/nt,0);
%        
%        
%        % replace spike value
%        % NB auxillary variables also replaced with low-pass filtered value 
%        % if identified as spike in U or V varible.
%        u(ii)=uf(ii); v(ii)=vf(ii); p(ii)=pf(ii); t(ii)=tf(ii); w(ii)=wf(ii); pit(ii)=pitf(ii); rol(ii)=rolf(ii);
%        uss(ii)=ussf(ii); vss(ii)=vssf(ii); wss(ii)=wssf(ii); ipow(ii)=ipowf(ii);
%        %despikedu(ii)=uf(ii);

       % NEED TO REMOVE NANS BEFORE FILTERING. OTHERWISE REPLACES WHOLE SERIES WITH NaNs.
       abc=find(~isnan(u));

       % low pass product again
       uf=mfilter(u(abc),1,1/nt,0);
       vf=mfilter(v(abc),1,1/nt,0);
       pf=mfilter(p(abc),1,1/nt,0);
       tf=mfilter(t(abc),1,1/nt,0);
       wf=mfilter(w(abc),1,1/nt,0);
       ussf=mfilter(uss(abc),1,1/nt,0);
       vssf=mfilter(vss(abc),1,1/nt,0);
       wssf=mfilter(wss(abc),1,1/nt,0);
       pitf=mfilter(pit(abc),1,1/nt,0);
       rolf=mfilter(rol(abc),1,1/nt,0);
       ipowf=mfilter(ipow(abc),1,1/nt,0);
       
       % 12h resolution
       tim=fix(jd(1)):0.5:fix(jd(end));
       vf1=interp1(jd(abc),vf,tim);
       uf1=interp1(jd(abc),uf,tim);
       pf1=interp1(jd(abc),pf,tim);
       tf1=interp1(jd(abc),tf,tim);
       wf1=interp1(jd(abc),wf,tim);
       ussf1=interp1(jd(abc),ussf,tim);
       vssf1=interp1(jd(abc),vssf,tim);
       wssf1=interp1(jd(abc),wssf,tim);
       pitf1=interp1(jd(abc),pitf,tim);
       rolf1=interp1(jd(abc),rolf,tim);
       ipowf1=interp1(jd(abc),ipowf,tim);


       figure(low_pass_fig)
       plot(tim-jd1,uf1); hold on;
       plot(tim-jd1,vf1,'-r');
       plot(tim-jd1,wf1,'-g');
       xlim([0 jd2-jd1]);
       set(gca,'xTickLabel',xticklabels);
       set(gca,'XTick',jdxticks-jd1);
       if plot_x_labels>0
            ylimits=get(gca,'ylim');
            y_label_pos=ylimits(1)-abs(ylimits(1)-ylimits(2))/10;
            for k=1:length(year_indexes)
                text((jd2-jd1)*(year_indexes(k)-1)/(length(xticklabels)-1),y_label_pos,num2str(xticks(year_indexes(k),1)),'FontSize',10);
            end
       end
       hold off;
       title(['SONTEK SN: ',sprintf('%3.3d',vecARG(i)), ' 12-hour resolution low-pass filtered currents'])
       legend({'u-component','v-component','w-component'})
       

       %replace NANs with -9999
       tf1(isnan(tf1))=-9999;
       pf1(isnan(pf1))=-9999;
       uf1(isnan(uf1))=-9999;
       vf1(isnan(vf1))=-9999;
       wf1(isnan(wf1))=-9999;
       ussf1(isnan(ussf1))=-9999;
       vssf1(isnan(vssf1))=-9999;
       wssf1(isnan(wssf1))=-9999;
       pitf1(isnan(pitf1))=-9999;
       rolf1(isnan(rolf1))=-9999;
       ipowf1(isnan(ipowf1))=-9999;
       
       % Save to rodb format
       gt=gregorian(tim);
       YY=gt(:,1); MM=gt(:,2); DD=gt(:,3); 
       if size(gt,2) == 6
           HH=hms2h(gt(:,4),gt(:,5),gt(:,6)); 
       else 
           HH= gt(:,4);
       end  
       fort = '%4d  %2d  %2d  %8.4f  %6.4f  %6.1f  %7.2f  %7.2f  %7.2f  %2.0d  %2.0d  %2.0d  %3.1d  %3.1d  %3.1d';
       cols = 'YY:MM:DD:HH:T:P:U:V:W:USS:VSS:WSS:PIT:ROL:IPOW';
       data=[YY MM DD HH tf1' pf1' uf1' vf1' wf' ussf1' vssf1' wssf1' pitf1' rolf1' ipowf1'];
       infovar=['Mooring:Latitude:Longitude:Columns:Start_Date:Start_Time:End_Date:End_Time:SerialNumber:' ...
           'WaterDepth:InstrDepth:MagDeviation'];
       rodbsave(outfile,infovar,fort,moor,lat,lon,cols,s_d,s_t,e_d,e_t,serialno,wd,depths(iiiARG(i),2),magdec,...
                data);
    
       j=j+1;
       disp('PAUSED - press any key to continue to next instrument')
       pause
    end
end

%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for mr=1:length(ifile);
%     fidi=[ifile{mr}];
%  
%     % make time to julian;
%     time=julian(yy,mm,dd,h);
%     dt=(time(2)-time(1))*24    % sampling interval
%     nt=round(40/dt);            % 40 hr lowpass
% 
%     % CORRECT FOR SOUND VELOCITY OF SEAWATER
%     % FIRST CALC Sound velocity as fxn of T P S
%  
%     S(1:length(t),1)=35; 
%     svel = sw_svel(S,t,p);
%     svel_old=1500;
%     corfac= (svel/svel_old);
% 
%     u = uraw .* corfac;
%     v = vraw .* corfac;
% 
%    
% 
% 
%    
%     
% end

function [ur,vr]=uvrot(u,v,rot)
% function [ur,vr]=uvrot(u,v,rot)
% rotate velocities
% by angle rot in deg
% M. Visbeck
rot = -rot*pi/180;
cr=cos(rot);
sr=sin(rot);
ur=u.*cr-v.*sr;
vr=u.*sr+v.*cr;
return
end

function [y,b]=mfilter(x,f0,f1,f2,n,ss)
% function [y,b]=mfilter(x,f0,f1,f2,n,ss);
%
% filter data vector x by FIR1 filter
% Hamming window second order used
% x timeseries
% f0 sample frequency
%
% f1 , f2 filter frequency
% f1=0 f2=highpass
% f1=lowpass f2=0
% f1 < f2 bandpass
% f1 > f2 bandstopp
%
%  if you like
% n filter elements 
% ss delta for subsample
%
%M. Visbeck 29.07.91
if ~exist('n','var'), n=f0/max(f1,f2); end
if ~exist('ss','var'), ss=1; end

%lp
if f2==0, Wn=f1/f0*2; ftype=''; 
%hp
elseif f1==0, Wn=f2/f0*2; ftype='high'; 
%bp
elseif f2>f1, Wn=[f1 f2]/f0*2; ftype=''; 
%bs
elseif f1>f2, Wn=[f2 f1]/f0*2; ftype='stop'; 
end

b=fir1(n,Wn,ftype);

[lc,lt]=size(x);

if ( lc==1 || lt==1 ) 
 y=filtfilt(b,1,x);
 y=y(1:ss:lt*lc);
else
 y=x;
 for i=1:lc
  y(i,:)=filtfilt(b,1,x(i,:));
 end
 y=y(:,1:ss:lt);
end
end

function [series,hf_series] = manual_spike_replace(series,jd_time,hf_series,fig_handle)
    % function to work with data brushing tool to manually edit spikes and
    % replace them with NaNs
    finished=1;
    while finished==1
        finished=input('Further manual editing required? 0=no, 1=yes: ');
        if finished==1
            disp('Manually select data using data brushing tool and create variable "spike" from bad data.')
            brush(fig_handle)
            pause
            for ij=1:length(spike)
                a=find(jd_time-jd_time(1)==spike(ij,1));
                series(a)=NaN;
                hf_series(a)=NaN;
            end
            figure(fig_handle)
            hold off
            plot(jd_time-jd_time(1),hf_series); 
            %clear spike
        end
    end
end

function [cleaned_series] = replace_spike(series,sampling_rate)
    % function to linearly interpolate across NaNs if less than 24hr gaps
    % if 24 hours or greater gaps then continue to replace interpolated 
    % value with NaNs once more
    % series = timeseries to replace NaNs in
    % sampling_rate = number of samples in 24 hours
    index=find(~isnan(series));
    
    % now find periods when consecutive NaNs greater than 24 hours
    abc=zeros(size(series));
    for i=1:length(series)+1-sampling_rate
        isnan_string='isnan(series(i';
        eval_string='isnan(series(i)) ';
        for j=1:sampling_rate-1
            eval_string=[eval_string ' & ' isnan_string '+' num2str(j) ')) '];
        end
        eval(['if(find(' eval_string '));   abc(i)=1;   end'])
    end
    nan_index24=find(abc>0);
    dummy_time=1:1:length(series);
    
    cleaned_series=interp1(dummy_time(index),series(index),dummy_time);
    cleaned_series(nan_index24:nan_index24+sampling_rate-1)=NaN;
    
end