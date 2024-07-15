function offset = wit_toffset(pd, cast, inst, C, P, cond_threshold, impact_var, cnv_time_correction)

% --------------------------------------------------------------------------
% 2. ---- water impact times: determine time offsets between ctd and mc ----
% --------------------------------------------------------------------------
% CTD impact

global MOORPROC_G

% water impact time ctd
[wit_ctd,~] = ctd_impact(pd.ctd_cnvfile,impact_var,cond_threshold,MOORPROC_G.cruise); 
wit_ctd(4)  = wit_ctd(4) - cnv_time_correction*24;

if size(wit_ctd,2) == 6
  wit_ctd_mdn = datenum(wit_ctd(1:6));
elseif  size(wit_ctd,2) == 4
  wit_ctd_mdn = datenum(wit_ctd,0,0);
end

%MicroCAT impact
ninst = size(C,2);

wit_mc_mdn = nan(1,ninst);
for inst = 1 : ninst
    if strcmp(calp.impact_var,'c')
        ii = find(C(:,inst) > cond_threshold);           % water impact mc
    elseif strcmp(calp.impact_var,'p')
        ii = find(P(:,inst) > cond_threshold);           % water impact mc
    end
    if sum(~isnan(C(:,inst)))==0
        disp(' ');
        error(['NO DATA FOR INSTRUMENT ' num2str(instr(inst)) ...
            '. MAKE SURE THAT A FILE calp.cast' num2str(cast) '_' num2str(instr(inst)) ...
            '.raw EXISTS IN cal_dip/microcat/calp.cast' num2str(cast) ...
            '/  OR REMOVE THE INSTRUMENT ' num2str(instr(inst))  ' FROM THE FILE ' ....
            pd.info_file])
    end
    if ~isempty(ii)
        ii = ii(1);
        wit_mc_mdn(inst)  = datnum(ii,inst);
    else
        warning('no good C or P data for inst %d',instr(inst))
    end
end

% time difference between start of instrument and ocean surface impact [s]
dwit_mc   = (wit_mc_mdn - start_mc_mdn)*24*3600;

% time offset  between CTD and MC: needed to compare bottle stop values
ii        = find(dwit_mc > 60);   % only consider mc with impact time > 60 for others
                                  % may have started after surface impact
impact_offset =  wit_ctd_mdn - wit_mc_mdn; % individual impact time offsets ctd - mc

if ~isempty(ii)
  offset    = wit_ctd_mdn -  median(wit_mc_mdn(ii));         % odecimal days
  [ohms(1),ohms(2),ohms(3)] = s2hms(offset*24*3600);
else
  offset = 0;
  [ohms(1),ohms(2),ohms(3)] = s2hms(offset*24*3600);  
end

[oh,om,os] = s2hms(impact_offset*86400);
fprintf(1,'\n Time offset of MicroCATs rel. to CTD:\n\n') 
fprintf(1,'  ID    HH    MM   SS\n')
fprintf(1,'  %d: %d  %d  %d\n',[instr';oh;om;round(os)]);
