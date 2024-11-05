function insitu_cal_osnap2(calp)
% compare Mooring sensors to CTD at bottle-stops for post-calibration
%
% uses read_rosfile.m, ctd_impact.m, mc_concorr.m, num_legend.m, pload.m (+subroutines)
%
% This version can be used for MicroCAT / Sontek Argonaut data
%
%  Input
%    - simultaneous Mooring sensor and CTD data (+ corresponding info.dat and mooring information)
%
%
%  Output
%
%  dt,dc,dp                   -- T / C / P  diff. MicroCAT - CTD at the bottlestops
%  dt_mcdep,dc_mcdep,dp_mcdep -- same as above, but interp. onto sensor deployment depths
%  dt_av,dc_av          -- same quantity, but aver. over pressure range (average_interval)
%

% Kanzow, adapted from 'microcat_insitu_cal.m'.
%         'microcat_insitu_cal.m' is obsolete and will not be seviced any more
%         03.01.06  Sontek / Argonaut data option added
% Loic Houpert,  - adapted the function read_rosfile to be able to read time in julian days
% 19.10.2015,    - adapted the function read_rosfile and ctd_impact to deal with cnv and ros files
%                 with a problem in the  header variable #start_time
%                - add reference to structure p_insitucal in which the parameters
%                 are defined (see osnap/users/loh/mcatpostcruisecalib/postcalib_process)
%               - add a parameter mc_cunit
%
% Lewis Drysdale, - bug fix with Julian date conversion of microcvat data
%     06.12.2022,     see issues in GitHub https://github.com/ScotMarPhys/m_moorproc_toolbox/issues/17
%
% Yvonne Firing   - updated for same issue
%     20.01.2022,     (as well as a Matlab version incompatibility in inputs to juliandate)
%


global MOORPROC_G
warning off
calp.jd_mdn = 1721058.5; %offset from julian (as calculated by matlab's juliandate) and matlab datenum
lat = 57; % set latitude for calculating pressure to set mc's to deployment depths
% ------ output / input files and directories ---------------

id_z_sn = all_inst_table(calp.sensor_id(1),-999,-999);
sensortype = char(id_z_sn.instl(1));
if isempty(sensortype)
    disp('Sensor type unknown - return')
    return
end
if calp.sensor_id(1) == 366
    calp.cond_threshold = 10; % dbar
    calp.impact_var     = 'p';
    calp.c_interval     = [-2 2]; % cond. dev.plot range
    calp.t_interval     = [-2 2]; % temp. dev.plot range
    calp.p_interval    = [-35  35];   % pres. dev. plot range
end

output_name = sprintf('%s_insitu_cal_%s_%s_cast%d_cond%d',sensortype,calp.cruise,calp.depl_period,calp.cast,calp.sensorselect);

% Data directories (and information about file/variable formats)
pd = moor_inoutpaths('cal_coef',calp.cast);


%-------------------------------------
% 1. --- load data  ------------------
%-------------------------------------

% ----MicroCAT ID -----------
[~,typ,instr,sn1] = rodbload(pd.info_file,'z:instrument:serialnumber:deployment');
m = (typ>=calp.sensor_id(1) & typ<=calp.sensor_id(end));
instr = instr(m);
sn1 = sn1(m);

% --- deployment depths --------
ein = fullfile(pd.coef_dir,[calp.depl_period '_deploymentdepths.dat']);
[dep,typ,ssn] = rodbload(ein,'z:instrument:serialnumber');
if isempty(dep)
    error(['!! MAKE sure the file ' ein ' exists or edit the filepath !!' ])
    return
end
% %make sure they are in the same order
% [~,ia,ib] = intersect(ssn,sn1,'stable');
% if length(ia)<length(ssn)
%     warning('some instruments in deployment file not in info.dat file')
% end
% dep = dep(ia); typ = typ(ia); ssn = ssn(ia);
% instr = instr(ib); clear sn1
% % now ssn, dep, typ, instr all match; saves looping later


% ---- CTD and Bottle data ----------------
fprintf(1,'\n Loading CTD  and bottle data ...\n')
[ctd, bottle] = load_ctdcaldip(pd, calp.cast, calp.sensorselect, calp.cnv_time_correction,calp.jd_mdn);


% --- MicroCAT and Seaphox -------
fprintf(1,'\n loading MicroCAT data ... \n\n')
ninst = length(instr);
for mct = 1 : ninst
    fname = sprintf('%s%4.4d%s',pd.mc_file,instr(mct),pd.mc_ext);
    [yy,mm,dd,hh,p,t,c,s]  = rodbload(fname,'YY:MM:DD:HH:P:T:C:S');
    lyy                = length(yy);
    mc.datnum(1:lyy,mct) = datenum(yy,mm,dd,hh,0,0);
    mc.T(1:lyy,mct)        = t;
    mc.C(1:lyy,mct)        = c;
    mc.P(1:lyy,mct)        = p;
    mc.S(1:lyy,mct)        = s;
end
mc.datnum(mc.datnum==0) = NaN;
m = isnan(mc.datnum);
mc.T(m) = NaN; mc.C(m) = NaN; mc.P(m) = NaN; mc.S(m) = NaN;
% check which variables have been measured by sensor
calp.pstat = double(sum(~isnan(mc.P))>0);
calp.cstat = double(sum(~isnan(mc.C))>0);
calp.tstat = double(sum(~isnan(mc.T))>0);
% Seaphox outputs salinity but not conductivity so back calculate
% conductivity if salinity is a variable
if calp.sensor_id==375
    if sum(~isnan(mc.C))==0 && sum(~isnan(mc.S))>0
        mc.C = gsw_C_from_SP(mc.S,mc.T,mc.P);
    end
end
missdat = sum(~isnan(mc.C))==0;
if sum(missdat)
    warning('NO DATA FOR INSTRUMENTS:')
    fprintf(1,'%d\n',instr(missdat))
    instm1 = instr(missdat); instm1 = instm1(1);
    error(['MAKE SURE THAT FILES cast%d_%d.raw (etc.) exist ...' ...
        'in cal_dip/microcat/cast%d/ \nor remove these instruments ...' ...
        'from moor/proc_calib/%s/cal_dip/cast%dinfo.dat', ...
        calp.cast,instm1,calp.cast,MOORPROC_G.cruise,calp.cast]);
end
% MicroCAT units
if strcmp(pd.mc_cunit,'S/m')
    mc.C = mc.C*10;
end
if strcmp(pd.ctdcnv_cunit,'S/m')
    ctd_cscale = 1/10; %***
else
    ctd_cscale = 1;
end

% --------------------------------------------------------------------------
% 2. ---- water impact times: determine time offsets between ctd and mc ----
% --------------------------------------------------------------------------
% CTD water impact time
[wit_ctd,~] = ctd_impact(pd.ctdcnv_file,calp.impact_var,calp.cond_threshold*ctd_cscale,MOORPROC_G.cruise);
wit_ctd(4)  = wit_ctd(4) - calp.cnv_time_correction*24;
if size(wit_ctd,2)==4
    wit_ctd = [wit_ctd 0 0];
end
wit_ctd_mdn = datenum(wit_ctd(1:6));
% wit_mc_mdn = NaN+mc.C(:,1); % produces array of data sample size not
% instrument size
wit_mc_mdn = nan(1,ninst);
if strcmp(calp.impact_var,'c')
    m = mc.C>calp.cond_threshold;
elseif strcmp(calp.impact_var,'p')
    m = mc.P>calp.cond_threshold; %***
end
for inst = 1 : ninst
    ii = find(m(:,inst),1,'first');
    if isempty(ii)
        warning('no good C or P data for inst %d',instr(inst))
    else
        wit_mc_mdn(inst)  = mc.datnum(ii,inst);
    end
end
% time difference between start of instrument and ocean surface impact [s]
start_mc_mdn  = mc.datnum(1,:);
dwit_mc   = (wit_mc_mdn - start_mc_mdn)*24*3600;
% time offset  between CTD and MC: needed to compare bottle stop values
m = (dwit_mc > 60);   % only consider mc with impact time > 60 for others
% may have started after surface impact
impact_offset =  wit_ctd_mdn - wit_mc_mdn; % individual impact time offsets ctd - mc
if sum(m)
    offset    = wit_ctd_mdn -  median(wit_mc_mdn(m));   % decimal days
else
    offset = 0;
end
[ohms(1),ohms(2),ohms(3)] = s2hms(offset*24*3600);
[oh,om,os] = s2hms(impact_offset*86400);
fprintf(1,'\n Time offset of MicroCATs rel. to CTD:\n\n')
fprintf(1,'  ID    HH    MM   SS\n')
fprintf(1,'  %d: %d  %d  %d\n',[instr';oh;om;round(os)]);


%-------------------------------------------------------
% 3. ---- extract data from bottle stops ---------------
%-------------------------------------------------------
[mc, ctd] = get_bottlestop_data(bottle, ctd, mc, pd, calp, ninst);


% --------------------------------------------------------------------
% ----------------- CALCULATE DEVIATIONS MC - CTD  -------------------
% --------------------------------------------------------------------

% -------- allocate sensors to deployment depths ----
mcdep = NaN+zeros(1,ninst); mcdep2 = mcdep;
for ii=1:length(instr)
    m = find(typ>=calp.sensor_id(1) & typ<=calp.sensor_id(end) & ssn == instr(ii));
    mcdep2(ii) = dep(m);
    mcdep(ii) = gsw_p_from_z(-mcdep2(ii),lat);
end

% --- compute dc, dt, dp

bottle.p0av    =  ctd.pav'*ones(1,ninst);
mc.dp          = mc.pav - (bottle.p0av) ;
bottle.t0av    = ctd.tav'*ones(1,ninst);
mc.dt          = mc.tav - (bottle.t0av);

pproblem    = find(rms(mc.dp) > calp.dp_tol);
kompx       = find(isnan(mc.P(1,:)));   % index of sensors without own pressure sensor
cavc        = mc_concorr(mc.cav(:,kompx),bottle.p0av(:,kompx),0);
mc.cav(:,kompx)= cavc;      % insert corrected C
cav_problem = mc_concorr(mc.cav(:,pproblem),bottle.p0av(:,pproblem),mc.pav(:,pproblem));
bottle.c0av    = ctd.cav'*ones(1,ninst);
mc.dc          = mc.cav - (bottle.c0av );
dc_pproblem  = cav_problem - bottle.c0av(:,pproblem);
average_interval = calp.average_interval(1):20:calp.average_interval(2);

% compute interpolated and averaged versions of dt,dc,dp
pp_count   = 1; % pressure problem index
MP = max(ctd.press) +100;
oknan = find(~isnan(bottle.p0av(:,1)));
p0 = [MP; bottle.p0av(oknan,1)];
if isempty(calp.p_interval)
    calp.p_interval = [0 ceil(round(MP/100)*100)];
end

dt_mcdep = diag(interp1(p0, [mc.dt(1,:); mc.dt(oknan,:)], mcdep(:)))';
dc_mcdep = diag(interp1(p0, [mc.dc(1,:); mc.dc(oknan,:)], mcdep(:)))';
dp_mcdep = diag(interp1(p0, [mc.dp(1,:); mc.dp(oknan,:)], mcdep(:)))';
ti = interp1(p0, [mc.dt(1,:); mc.dt(oknan,:)], average_interval(:));
ci = interp1(p0, [mc.dc(1,:); mc.dc(oknan,:)], average_interval(:));
dt_av = mean(ti, 'omitnan'); dt_sd = std(ti, 'omitnan');
dc_av = mean(ci, 'omitnan'); dc_sd = std(ci, 'omitnan');

dp_mcdep_ext = NaN+dp_mcdep;
dc_av_pproblem = dp_mcdep_ext;
for i = 1:ninst
    if isnan(dp_mcdep(i)) % extrapolate if deployment pressure > max. pressure of cast
        pol             = polyfit([bottle.p0av(:,1)'],[mc.dp(:,i)'],3);
        dp_mcdep_ext(i) = polyval(pol,mcdep(i));
    end

    if ismember(i,pproblem) % if MC press. bad, use ctd press. to correct conduct.
        dc_av_pproblem(i) = mean(interp1([MP bottle.p0av(oknan,1)'],...
            [dc_pproblem(1,pp_count) dc_pproblem(oknan,pp_count)'],average_interval),'omitnan');
        pp_count = pp_count + 1; %***???
    end

end

% ---- INTERNAL SUBROUTINES ----------------------

% ------ save ----------------------------------
%identify numbr of stops (copied out of get_bottlestop_data.m subroutine)
bst2 = [1; find(diff(bottle.p)<-calp.bottlestop_dpmin)+1];
bot_start = bottle.datnum(bst2);
nstop = length(bot_start);

instr_id = repmat(instr(:)',nstop,1);

% -------------------------------------------------
% ------ GRAPHICS --------------------------------
% -------------------------------------------------

% plot temp_deviations

col  = ['brgkmcybrgkmcybrgkmcybrgkmcy'];
mrk  = ['ddddddd+++++++xxxxxxxsssssss'];
lin  = ['----------------------------'];

% ----------  Graphics  ----------------------

if sum(~isnan(mc.dp))
    sub = 3;
    figure(2);subplot(sub,1,1);hold off;subplot(sub,1,2);hold off;
    subplot(sub,1,3);hold off;hold off;
else
    sub = 2;
    figure(2);subplot(sub,1,1);hold off;subplot(sub,1,2);hold off;hold off
end


%plot values at all and highlight those at MicroCAT deployment depths
for i = 1 : ninst
    subplot(1,sub,1)
    plot(bottle.p0av(:,i),mc.dt(:,i),'color',col(i),'marker',mrk(i),'linestyle',lin(i)); hold on
    subplot(1,sub,2)
    plot(bottle.p0av(:,i),mc.dc(:,i),'color',col(i),'marker',mrk(i),'linestyle',lin(i)); hold on
    subplot(1,sub,3)
    plot(bottle.p0av(:,i),mc.dp(:,i),'color',col(i),'marker',mrk(i),'linestyle',lin(i)); hold on
end

subplot(1,sub,1)
num_legend(instr','''best''',5)  %legend
%num_legend(instr','southeast',5)  %legend

for i = 1:ninst
    subplot(1,sub,1)
    plot(mcdep(i),dt_mcdep(i),'linestyle','none','marker',mrk(i),'markersize',11,'color',col(i))
    subplot(1,sub,2)
    plot(mcdep(i),dc_mcdep(i),'linestyle','none','marker',mrk(i),'markersize',11,'color',col(i))
    subplot(1,sub,3)
    if ~isnan(dp_mcdep)
        plot(mcdep(i),dp_mcdep(i),'linestyle','none','marker',mrk(i),'markersize',11,'color',col(i))
    else
        plot(mcdep(i),dp_mcdep_ext(i),'linestyle','none','marker',mrk(i),'markersize',11,'color',col(i))
    end
end


subplot(1,sub,1)
grid on
set(gca,'Fontsize',12,'xlim',calp.p_interval,'ylim',calp.t_interval)
xlabel('pressure [dbar]')
ylabel(['temp. diff. [C]'])
title([sprintf('Deviations %s-CTD //%s',upper(sensortype),date) ])
set(gca,'color',[.7 .7 .7])
subplot(1,sub,2)
grid on
set(gca,'Fontsize',12,'xlim',calp.p_interval,'ylim',calp.c_interval)
xlabel('pressure [dbar]')
ylabel('cond. diff. [mS/cm]')
title(['Cruise: ',calp.cruise,'  Cast: ',num2str(calp.cast)])
set(gca,'color',[.7 .7 .7])

subplot(1,sub,3)
grid on
%   set(gca,'Fontsize',12,'xlim',p_interval,'ylim',dp_interval)
set(gca,'Fontsize',12,'ylim',calp.dp_interval)
xlabel('pressure [dbar]')
ylabel('pres. diff. [dbar]')
set(gca,'color',[.7 .7 .7])


set(figure(2),'Paperunits','centimeters','Paperposition',[0 0 29 21])
set(0,'currentfigure',figure(2))

figname = fullfile(pd.mc_dir,sprintf('cast%d',calp.cast),output_name);
print([figname '.png'],'-dpng')

set(0,'currentfigure',figure(2))
orient landscape
print([figname '.ps'],'-dpsc')


if ~isempty(pproblem)
    figure(4);clf;hold on

    for i = 1 : length(pproblem)

        plot(bottle.p0av(:,i),dc_pproblem(:,i),mrk(:,pproblem(i)))
        hold on
        %ylim(c_interval)
        xlabel('pressure [dbar]')
        ylabel('cond. dev. [mS/cm]')
        xlim(calp.p_interval)
        num_legend(instr(pproblem),'''best''',5)
        grid on
        title('MicroCATs pressure problems: cond. diff. MC-CTD (corrected by CTD pressure)')
        set(gca,'FontSize',12)
    end

end
% -- TEXT OUTPUT ---------
outfile=fullfile(pd.mc_dir,['cast',num2str(calp.cast)],[output_name '.txt']);
disp(['Printing cal coefs to ' outfile])
fid = fopen(outfile,'w');
fprintf(fid,['%s',calp.cruise,'  cast:',num2str(calp.cast),' MicroCAT - CTD // processing date: ',date,' \n'],'%');

fprintf(fid,'%sID   CAST   DEPTH RANGE  DEPTH   dt  dt_sd     dc  dc_sd       dp    |  dp_ext | dc_pproblem\n','%');
fprintf(fid,'%sID   CAST    T/C [m]     P [m]   [K]   [K] [mS/cm]   [ms/cm]  [dbar] |  [dbar] |   [mS/cm] \n','%');

caldata = [instr repmat([calp.cast average_interval([1 end])],ninst,1) ...
    mcdep2' dt_av' dt_sd' dc_av' dc_sd' dp_mcdep' dp_mcdep_ext' dc_av_pproblem'];

fprintf(fid,' %4.4d  %2.2d    %4.4d-%4.4d  %4.4d    %5.4f %5.4f   %5.4f  %5.4f   %3.1f  |  %3.1f  | %5.4f\n',caldata');

fprintf(fid,['\n At nominal depth of instrument \n' ],'%');
fprintf(fid,'%sID   CAST   DEPTH   dt   dc     dp    |  dp_ext \n','%');
fprintf(fid,'%sID   CAST   P [m]   [K]  [mS/cm]  [dbar] |  [dbar] \n','%');

caldata2     = ...
    [instr repmat(calp.cast,ninst,1) ...
    mcdep2' dt_mcdep' dc_mcdep' dp_mcdep' dp_mcdep_ext' ];

fprintf(fid,' %4.4d  %2.2d   %4.4d    %5.4f    %5.4f  %3.1f  |  %3.1f  \n',caldata2');
fclose(fid);

