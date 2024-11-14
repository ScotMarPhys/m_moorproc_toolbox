function [ctd, varargout] = load_ctddata(pd, cast, sensorselec, cnv_time_correction,jd_mdn)
%
% ctd = load_ctddata(pd, cast, sensorselect);
% [ctd, bottle] = load_ctddata(pd, cast, sensorselect);
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
                    ctd.cond  = cnv_cor_save(i).C; % better sensor
                    ctd.temp  = cnv_cor_save(i).T;
                    ctd.press = cnv_cor_save(i).P;
                    ctd.time =  cnv_cor_save(i).JD0 + julian([2015 1 1 0 0 0]) - 1;
                    ctd.datnum = ctd.time-jd_mdn; ctd = rmfield(ctd,'time');
            	end
            else
                if cnv_cor_save(i).station == cast
                    if strcmp('kn221-02',mg.cruise) || strcmp('kn221-03',mg.cruise)
                        ctd.cond  = cnv_cor_save(i).c1S; % better sensor
                        ctd.temp  = cnv_cor_save(i).t190C;
                        ctd.press = cnv_cor_save(i).prDM;

                        ctd_time_ori = julian([cnv_cor_save(i).gtime(1:3) 0]);
                        jul_day_frac = cnv_cor_save(i).timeJ - floor(cnv_cor_save(i).timeJ);
                        ctd.time = jul_day_frac + ctd_time_ori;
                        ctd.datnum = ctd.time-jd_mdn; ctd = rmfield(ctd,'time');
                        ctd_time_ori = ctd_time_ori-jd_mdn;
                    else
                        ctd.cond  = cnv_cor_save(i).conductivity;
                        ctd.temp  = cnv_cor_save(i).temperature;
                        ctd.press = cnv_cor_save(i).pressure;
                        %if strcmp(ctd_cunit,'S/m')
                        ctd.cond = ctd.cond;%*10;

                        ctd_time_ori = datenum(cnv_cor_save(i).gtime(1:6));
                        ctd.datnum = cnv_cor_save(i).elap_time_sec/86400 + ctd_time_ori;
                    end
                end
            end
        end

        clear cnv_cor_save
        pack

    %%%%%%% PSTAR %%%%%%%
    case 'pstar'
    pd.ctd_file
    [ctd, h]  = pload(pd.ctd_file,'press temp cond time','silent');

    year   = floor(h.iymd/10000);
    month  = floor((h.iymd - year*10000)/100);
    day    = h.iymd -year*10000 - month*100;
    hour   = h.ihms;

    ctd_time_ori = datenum(year+h.icent, month, day,hour,0,0);
    ctd.datnum       = ctd.time/86400 + ctd_time_ori; % ctd time in datenum
    ctd = rmfield(ctd,'time');

    %%%%%%% MSTAR %%%%%%%
    case 'mstar'
        ctd_file=pd.ctd1hz_file

        switch mg.cruise
            case {'d344','d359'}
                dd  = netcdf(pd.ctd1hz_file); %,'press temp cond time','silent');
                d = struct('press',dd{'press'}(:),'temp',dd{'temp'}(:),'cond',dd{'cond'}(:),'time',dd{'time'}(:));
                year   = dd.data_time_origin(1); %floor(h.iymd/10000);
                month  = dd.data_time_origin(2); %floor((h.iymd - year*10000)/100);
                day    = dd.data_time_origin(3); %h.iymd -year*10000 - month*100;
                hour   = dd.data_time_origin(4) + (dd.data_time_origin(5)+(dd.data_time_origin(6)/60))/60 ; %h.ihms;

                ctd_time_ori = datenum(year,month,day,hour,0,0);
                ctd.datnum       = ctd.time/86400 + ctd_time_ori; % ctd time in datenum
                ctd = rmfield(ctd,'time');
                % if strcmp(mg.cruise,'jc064')
            case {'dy120','ar304','dy078','dy053','pe399'}
                if sensorselec==1
                    [ctd, h]=mload(pd.ctd1hz_file,'time','press','temp1','cond1',' ','q');
                    ctd.cond = ctd.cond1;
                    ctd.temp = ctd.temp1;
                elseif sensorselec == 2
                    [ctd, h]=mload(pd.ctd1hz_file,'time','press','temp2','cond2',' ','q');
                    ctd.cond = ctd.cond2;
                    ctd.temp = ctd.temp2;
                end
                ctd.datnum=datenum(h.data_time_origin)+ctd.time/86400;
                ctd = rmfield(ctd,'time');
            case 'jc238'
                [ctd, h] = mload(pd.ctd1hz_file,'time','press','temp','cond',' ','q');
                ctd.datnum=datenum(h.data_time_origin)+ctd.time/86400;
                ctd = rmfield(ctd,'time');
            otherwise
                if MOORPROC_G.YEAR>=2024
                    [ctd, h] = mload(pd.ctd1hz_file,'time','press','temp','cond',' ','q');
                    ctd.datnum = m_commontime(ctd, 'time', h, 'datenum');
                    ctd = rmfield(ctd,'time');
                else
                    dd  = netcdf(pd.ctd1hz_file); %,'press temp cond time','silent');
                    d = struct('press',ncread(ctd_file,'press'),...
                        'temp',ncread(ctd_file,'temp'),...
                        'cond',ncread(ctd_file,'cond'),...
                        'time',ncread(ctd_file,'time'));
                    ctd.temp=ctd.temp(:);
                    ctd.cond=ctd.cond(:);
                    ctd.press=ctd.press(:);
                    ctd.time=ctd.time(:);
                    ctd_time_ori = datenum(dd.data_time_origin);
                    ctd.datnum       = ctd.time/86400 + ctd_time_ori; % ctd time in datenum
                    ctd = rmfield(ctd,'time');
                end
        end

end

if nargout>1

    if isfield(pd, 'bottle_file')
        if isfile(pd.bottle_file)
            if contains(pd.bottle_file,'.btl')
                bottle = read_botfile(pd.bottle_file);
                for jjj=1:length(bottle)
                    bottle(jjj).p = bottle(jjj).pav;
                end
                bottle.datnum = [datenum(bottle.yy,bottle.mm,bottle.dd,bottle.hh)]';
            elseif contains(pd.bottle_file,'.ros')
                bottle            = read_rosfile(pd.bottle_file);
                bottle.datnum = bottle.jd - jd_mdn;
            elseif contains(pd.bottle_file,'fir_')
                [db, hb] = mload(pd.bottle_file,'/');
                %bottle = read_blfile(pd.bottle_file); %***create?
                bottle.datnum = m_commontime(db, 'utime', hb, 'datenum');
                bottle.p = db.upress;
            end
            bottle.firing='BottlesFired';
        else
            btlcont=input('No bottle file found. Proceed using CTD stop time only? (Y/N) ','s')
            if strcmpi(btlcont,'y')
                bottle.datnum =[];
                bottle.start_time=[];
                bottle.firing='NoBottleFired';
            else
                return
            end
        end

        bottle.datnum         = bottle.datnum - cnv_time_correction; %***what about for cnv?
        if isfield(bottle,'start_time')
        bottle.start_time = bottle.start_time - cnv_time_correction;
        end
        varargout{1} = bottle;
    else
        disp('PATH TO BOTTLE FILES NOT DEFINED')
        varargout{1} = [];
    end

end
