% function nortek_raw2use_02(moor,varargin)
%
% basic preprocessing for Nortek Aquadopp data
% 
% required inputs: moor - mooring name e.g 'eb3_1_200406'
%
% optional inputs: procpath - path to proc directory if not using standard
%                             paths
%                  outpath - path to output processed files to - otherwise
%                            outputs to directory function run from
%                  plot_interval - matrix of start and end dates for plot
%                                  e.g. [2004 02 01 00; 2005 06 01 00]
%                                  dates are:- yyyy mm dd hh
%
% features
%      1. eliminate launching and recovery period
%      2. save data to rodb file
%
% uses timeaxis.m, auto_filt.m, julian.m, rodbload.m, rodbsave.m

% 15/01/2007 - DR modified from argocat_raw2use_003
% 20/12/2012 - DR added hdg, pit, rol and sig strength plot output
% 06/01/2016 - Loic Houpert, SAMS, add option to detect the matlab version use
% to prevent crash in the suptitle function (Since R2014n, graphical functions
% changed)
%            - Loic Houpert, add png export for the saved plots

function nortek_raw2use_02(moor,varargin)

global MOORPROC_G

if nargin==0
    help nortek_raw2use_02;
    return
end

% check for optional arguments
if nargin==1
    pd = moor_inoutpaths('nor',moor);
else
    pd = varargin{1};
end

a = find(strcmp('plot_interval', varargin));
if a>0 && ~isempty(varargin{a+1})
    plot_interval=reshape(varargin{a+1},4,2)';
else
    plot_interval=[];
end

operator = MOORPROC_G.operator;
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(pd.infofile,'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

vec=find((id==368|id==370)); % Nortek id number

sn=sn(vec);
z = z(vec);

fid_stat= fopen(pd.stage2log,'a');
fprintf(fid_stat,['Processing steps taken by ' mfilename ':\n']);
fprintf(fid_stat,'  1. eliminate lauching and recovery period\n');
fprintf(fid_stat,'  2. save data to rdb file\n');
fprintf(fid_stat,'\n Operated by:%s on %s\n',operator,datestr(clock)); 
fprintf(fid_stat,['        Norteks in Mooring ',moor,'\n\n\n']);

dummy    = -9999;

% Determine plot_interval if not input to function
if isempty(plot_interval)
    plot_interval = [s_d(1) s_d(2)-1 s_d(3) 0; e_d(1) e_d(2)+1 e_d(3) 0];
    if plot_interval(1,2)==0
        plot_interval(1,2)=12; plot_interval(1,1)=plot_interval(1,1)-1;
    end
    if plot_interval(2,2)==13
        plot_interval(2,2)=1; plot_interval(2,1)=plot_interval(2,1)+1;
    end
end


%-----------------------------------------
% --- preprocessing loop -------------------
% ----------------------------------------

jd_s  = julian(s_d(1),s_d(2),s_d(3),s_t(1)+s_t(2)/60);  % start time
jd_e  = julian(e_d(1),e_d(2),e_d(3),e_t(1)+e_t(2)/60);  % end time

%all_z=[combo_z individual_z];
%all_sn=[combo_sn individual_sn];

for proc = 1 : length(vec)
    columns = 'YY:MM:DD:HH:T:P:U:V:W:HDG:PIT:ROL:USS:VSS:WSS:IPOW:CS:CD';
    indep  = z(proc);
    infile  = fullfile(pd.stage1path,sprintf(pd.stage1form,sn(proc)));

    if exist(infile,'file')==0
        disp(['infile: ' infile ' does not exist.'])

    else
        outfile = fullfile(pd.stage2path,sprintf(pd.stage2form,sn(proc)));
        fprintf(fid_stat,'Serialnumber %d \n',sn(proc));
        fprintf(fid_stat,'Infile %s \n',infile);
        fprintf(fid_stat,'Outfile %s \n',outfile);

        [YY,MM,DD,HH,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd] = ...
            rodbload(infile,[columns]);

        %------------------------------------------
        %----- cut off launching and recovery period
        %------------------------------------------
        disp('cut off launching and recovery period')

        jd  = julian(YY,MM,DD,HH);
        %gregorian(jd);
        %         jd_e
        %         gregorian(jd_e)
        %         jd_s
        %         gregorian(jd_s)

        ii  = find(jd <= jd_e & jd >= jd_s );

        YY=YY(ii);MM=MM(ii);DD=DD(ii);HH=HH(ii);
        t=t(ii);p=p(ii);
        u=u(ii);v=v(ii);w=w(ii);
        hdg=hdg(ii);pit=pit(ii);rol=rol(ii);
        uss=uss(ii);vss=vss(ii);wss=wss(ii);
        ipow=ipow(ii); cs=cs(ii); cd=cd(ii);
        jd  = jd(ii);

        cycles     = length(ii);
        Y = [YY(1) MM(1) DD(1)];
        Start_Time = HH(1);
        End_Date = [YY(cycles) MM(cycles) DD(cycles)];
        End_Time = HH(cycles);


        %------------------------------------------
        %---- fill time gaps  with dummy
        %------------------------------------------
        disp(' fill time gaps  with dummy')

        djd = diff(jd);           % time step
        sr  = median(djd);        % sampling interval
        ii  = djd > 1.5*sr;  % find gaps
        gap = round(djd(ii)/sr)-1;
        addt= [];



        %-----------------------------------------------------
        %  write output to logfile ---------------------------
        %-----------------------------------------------------

        fprintf(fid_stat,'Operation interval: %s  to  %s\n', ...
            datestr(gregorian(jd(1))),datestr(gregorian(jd(end)) ));
        fprintf(fid_stat,'\n');

        %-----------------------------------
        %--- write data to rodb format -----
        %-----------------------------------

        disp(['writing data to ',outfile])
        fort =['%4.4d %2.2d %2.2d  %6.4f  %4.2f %7.3f  %4.1f %4.1f %4.1f  '...
            '%4.1f %4.1f %4.1f  %2d %2d %2d  %2.1f  %4.1f %5.2f'];

        rodbsave(outfile,...
            'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
            fort,...
            lat,lon,columns,s_d,Start_Time,sn(proc),mr,wd,indep,End_Date,End_Time,...
            [ YY,MM,DD,HH,t,p,u,v,w,hdg,pit,rol,uss,vss,wss,ipow,cs,cd]);


        % aurelie: i have replace Start_Date by s_d and it works !!!!

        %%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%
        %who
        jd0 = julian(-1,1,1,24);
        jd1 = julian(plot_interval(1,:))-jd0;
        jd2 = julian(plot_interval(2,:))-jd0;
        sampling_rate = 1/median(diff(jd));

        STR = ['Temperature [deg C]   ';
            'Pressure [dbar]       ';
            'Zonal Velocity [cm/s] ';
            'Merid. Velocity [cm/s]';
            'Vert. Velocity [cm/s] '];
        STR2 = {'Heading [deg]';
            'angle [deg]';
            'sig. str. [counts]'};
        VAR1= ['t';'p';'u';'v';'w'];
        VAR2= {'hdg';'pit';'rol';'uss';'vss';'wss'};
        panels=5;
        panels2=3;

        figure;
        for sub = 1 : 5
            eval(['var1 = ',VAR1(sub),';'])
            var2=[];
            var3=[];
            ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR(sub,:),sub,[jd1 jd2],'n',panels);

        end

        subplot(5,1,1)
        title(['Nortek s/n: ',num2str(sn(proc)), ...
            '; Target Depth: ',num2str(indep)])

        orient tall

        %eval(['print -depsc2 -tiff ',outfile,'.eps'])
        eval(['print -dpng ',outfile,'.png'])

        figure;

        for sub = 1 : 5
            eval(['var1 = ',VAR1(sub),';'])
            var2=[];
            var3=[];
            ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR(sub,:),sub,[jd1 jd2],'y',panels);
        end

        subplot(5,1,1)
        title(['Nortek s/n: ',num2str(sn(proc)), ...
            '; Target Depth: ',num2str(indep)])
        orient tall
        %       eval(['print -depsc ',outfile,'.filtered.eps'])
        eval(['print -dpng ',outfile,'.filtered.png'])

        % plot of diagnostics info
        figure;
        subplot(3,1,1)
        eval(['var1 = ',VAR2{1},';'])
        var2=[];
        var3=[];
        ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR2{1},1,[jd1 jd2],'n',panels2);
        legend({'hdg'})

        subplot(3,1,2)
        eval(['var1 = ',VAR2{2},';'])
        eval(['var2 = ',VAR2{3},';'])
        var3=[];
        ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR2{2},2,[jd1 jd2],'n',panels2);
        legend({'pit','rol'})

        subplot(3,1,3)
        eval(['var1 = ',VAR2{4},';'])
        eval(['var2 = ',VAR2{5},';'])
        eval(['var3 = ',VAR2{6},';'])

        ok = plot_timeseries(jd,var1,var2,var3,sampling_rate,STR2{3},3,[jd1 jd2],'n',panels2);
        legend({'beam1','beam2','beam3'})


        suptitle(['Nortek s/n: ',num2str(sn(proc)), ...
            '; Target Depth: ',num2str(indep)])

        orient tall
        % eval(['print -depsc2 -tiff ',outfile,'_diagnostics.eps'])
        eval(['print -dpng ',outfile,'_diagnostics.png'])
    end % if exist(infile)==0
end  % for proc=1:length(combo_sn)+length(individual_sn) loop

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = plot_timeseries(jd,var1,var2,var3,sr,str,sub,jdlim,filt,panels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot time series


jd0 = julian(-1,1,1,24);

i1    = find(~isnan(var1) & var1~=0);

if strcmp(filt,'y')
    var1  = auto_filt(var1(i1),sr,1/2,'low',4);
else
    var1  = var1(i1);
end

if ~isempty(var2)
    i2    = find(~isnan(var2) & var2~=0);
    if strcmp(filt,'y')
        var2  = auto_filt(var2(i2),sr,1/2,'low',4);
    else
        var2  = var2(i2);
    end
end

if ~isempty(var3)
    i3    = find(~isnan(var3) & var3~=0);
    if strcmp(filt,'y')
        var3  = auto_filt(var3(i3),sr,1/2,'low',4);
    else
        var3  = var3(i3);
    end
end

subplot(panels,1,sub);

plot(jd(i1)-jd0,var1)
hold on

if ~isempty(var2)
    plot(jd(i2)-jd0,var2,'r')
end

if ~isempty(var3)
    plot(jd(i3)-jd0,var3,'g')
end

ylabel(str)
grid on
xlim([jdlim])
datetick('x',12,'keeplimits')

ok=1;

end
