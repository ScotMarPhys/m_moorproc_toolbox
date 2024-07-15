function [d, varargout] = load_ctddata(pd, cast, sensorselec)
%
% d = load_ctddata(pd, cast, sensorselect);
% [d, bottle] = load_ctddata(pd, cast, sensorselect);
%
% load CTD data from a cruise for microcat calibration
% optionally load bottle data

global MOORPROC_G
mg = MOORPROC_G;

switch pd.ctdformat

    %%%%%%% AOML %%%%%%%
    case 'aoml'
        load(pd.ctd_file)
        disp('aoml')
        switch mg.cruise
            case {'kn221-02','kn221-03'}
                cnv_cor_save=sbe;
            case 'pe400'
                cnv_cor_save=ctd_cal;
        end
        for i = 1: length(cnv_cor_save)
            if strcmp('pe400',mg.cruise)
            	strf = regexp(cnv_cor_save(i).file,'_','split');
            	if ~isempty(strfind(strf{2},num2str(cast)))
                    d.cond  = cnv_cor_save(i).C; % better sensor
                    d.temp  = cnv_cor_save(i).T;
                    d.press = cnv_cor_save(i).P;
                    d.time =  cnv_cor_save(i).JD0 + julian([2015 1 1 0 0 0]) - 1;
                    d.datnum = d.time-jd_mdn; d = rmfield(d,'time');
            	end
            else
                if cnv_cor_save(i).station == cast
                    if strcmp('kn221-02',mg.cruise) || strcmp('kn221-03',mg.cruise)
                        d.cond  = cnv_cor_save(i).c1S; % better sensor
                        d.temp  = cnv_cor_save(i).t190C;
                        d.press = cnv_cor_save(i).prDM;

                        ctd_time_ori = julian([cnv_cor_save(i).gtime(1:3) 0]);
                        jul_day_frac = cnv_cor_save(i).timeJ - floor(cnv_cor_save(i).timeJ);
                        d.time = jul_day_frac + ctd_time_ori;
                        d.datnum = d.time-jd_mdn; d = rmfield(d,'time');
                        ctd_time_ori = ctd_time_ori-jd_mdn;
                    else
                        d.cond  = cnv_cor_save(i).conductivity;
                        d.temp  = cnv_cor_save(i).temperature;
                        d.press = cnv_cor_save(i).pressure;
                        %if strcmp(ctd_cunit,'S/m')
                        d.cond = d.cond;%*10;

                        ctd_time_ori = datenum(cnv_cor_save(i).gtime(1:6));
                        d.datnum = cnv_cor_save(i).elap_time_sec/86400 + ctd_time_ori;
                    end
                end
            end
        end

        clear cnv_cor_save
        pack

    %%%%%%% PSTAR %%%%%%%
    case 'pstar'
    pd.ctd_file
    [d, h]  = pload(pd.ctd_file,'press temp cond time','silent');

    year   = floor(h.iymd/10000);
    month  = floor((h.iymd - year*10000)/100);
    day    = h.iymd -year*10000 - month*100;
    hour   = h.ihms;

    ctd_time_ori = datenum(year+h.icent, month, day,hour,0,0);
    d.datnum       = d.time/86400 + ctd_time_ori; % ctd time in datenum
    d = rmfield(d,'time');

    %%%%%%% MSTAR %%%%%%%
    case 'mstar'
        pd.ctd_file

        switch mg.cruise
            case {'d344','d359'}
                dd  = netcdf(pd.ctd_file); %,'press temp cond time','silent');
                d = struct('press',dd{'press'}(:),'temp',dd{'temp'}(:),'cond',dd{'cond'}(:),'time',dd{'time'}(:));
                year   = dd.data_time_origin(1); %floor(h.iymd/10000);
                month  = dd.data_time_origin(2); %floor((h.iymd - year*10000)/100);
                day    = dd.data_time_origin(3); %h.iymd -year*10000 - month*100;
                hour   = dd.data_time_origin(4) + (dd.data_time_origin(5)+(dd.data_time_origin(6)/60))/60 ; %h.ihms;

                ctd_time_ori = datenum(year,month,day,hour,0,0);
                d.datnum       = d.time/86400 + ctd_time_ori; % ctd time in datenum
                d = rmfield(d,'time');
                % if strcmp(mg.cruise,'jc064')
            case {'dy120','ar304','dy078','dy053','pe399'}
                if sensorselec==1
                    [d, h]=mload(pd.ctd_file,'time','press','temp1','cond1',' ','q');
                    d.cond = d.cond1;
                    d.temp = d.temp1;
                elseif sensorselec == 2
                    [d, h]=mload(pd.ctd_file,'time','press','temp2','cond2',' ','q');
                    d.cond = d.cond2;
                    d.temp = d.temp2;
                end
                d.datnum=datenum(h.data_time_origin)+d.time/86400;
                d = rmfield(d,'time');
            case 'jc238'
                [d, h] = mload(pd.ctd_file,'time','press','temp','cond',' ','q');
                d.datnum=datenum(h.data_time_origin)+d.time/86400;
                d = rmfield(d,'time');
            otherwise
                if MOORPRC_G.YEAR>=2024
                    [d, h] = mload(pd.ctd_file,'time','press','temp','cond',' ','q');
                    d.datnum = m_commontime(d, 'time', h, 'datenum');
                    d = rmfield(d,'time');
                else
                    dd  = netcdf(pd.ctd_file); %,'press temp cond time','silent');
                    d = struct('press',ncread(ctd_file,'press'),...
                        'temp',ncread(ctd_file,'temp'),...
                        'cond',ncread(ctd_file,'cond'),...
                        'time',ncread(ctd_file,'time'));
                    d.temp=d.temp(:);
                    d.cond=d.cond(:);
                    d.press=d.press(:);
                    d.time=d.time(:);
                    ctd_time_ori = datenum(dd.data_time_origin);
                    d.datnum       = d.time/86400 + ctd_time_ori; % ctd time in datenum
                    d = rmfield(d,'time');
                end
        end

end

if nargout>1

    if contains(pd.bottle_file,'.btl')
        bottle = read_botfile(pd.bottle_file);
        for jjj=1:length(bottle)
            bottle(jjj).p = bottle(jjj).pav;
        end
        bottle.datnum = [datenum(bottle.yy,bottle.mm,bottle.dd,bottle.hh)]';
    elseif contains(pd.bottle_file,'.ros')
        bottle            = read_rosfile(pd.bottle_file);
        bottle.datnum = bottle.jd - jd_mdn;
    else
        bottle = [];
        disp('PATH TO BOTTLE FILES NOT DEFINED')
    end
    
    varargout{1} = bottle;

end
