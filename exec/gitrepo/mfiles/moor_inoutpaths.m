function pd = moor_inoutpaths(datatype,loc)
% predirs = moor_inoutpaths(datatype,loc)
%
% predirs = moor_inoutpaths('microcat_cal_dip',cast)
% predirs = moor_inoutpaths('microcat',moor)
% predirs = moor_inoutpaths('nor',moor)
% predirs = moor_inoutpaths('bpr',moor)
% predirs = moor_inoutpaths('adcp',moor)
% predirs = moor_inoutpaths('cal_coef',cast)
%
% set input and output subdirectories (relative to data and reports
% top-level directories defined by moor_setup and held in MOORPROC_G) and
% filename prefixes to be accessed by successive stages of processing
%
% inputs are datatype and loc, where loc is the CTD cast number or the
% mooring deployment
%   e.g. 
%   moor_inoutpaths('microcat_cal_dip', 3)
%   moor_inoutpaths('microcat_cal_dip', '3')
%   moor_inoutpaths('nortek', 'ebh3_15_2022')

% if you add more cases for older cruises that don't follow the same
% conventions, please put them farther down, below the second 
% "switch datatype" line (search for mg.YEAR<2022) 


global MOORPROC_G
mg = MOORPROC_G;

switch datatype

    case 'microcat_cal_dip'
        cast = loc;
        if isnumeric(cast)
            castn = cast;
            cast = num2str(cast);
        else
            castn = str2double(cast);
        end
        pd.rawpath = fullfile(mg.moordatadir, 'raw', mg.cruise, 'microcat_cal_dip', ['cast' cast]);
        pd.infofile = fullfile(mg.moordatadir, 'proc_calib', mg.cruise, 'cal_dip', ['cast' cast 'info.dat']);
        pd.stage1path = fullfile(mg.moordatadir, 'proc_calib', mg.cruise, 'cal_dip', 'microcat', ['cast' cast]);
        pd.stage1form = ['cast' cast '_%4.4d.raw'];
        pd.stage1log = fullfile(pd.stage1path,'microcat2rodb.log');
        pd.stage1fig = fullfile(mg.reportdir,'figs','caldip');
        pd.ctdfile = fullfile(mg.ctddir,sprintf('ctd_%s_%03d_psal.nc',mg.cruise_ctd,castn));
        pd.stage2path = fullfile(mg.moordatadir, 'proc_calib'); %***
        pd.stage2fig = fullfile(mg.reportdir,'figs','caldip',['microcat_check_cast_' cast '_plot']);
        pd.stage2log = fullfile(mg.reportdir,'stats',['microcat_check' cast '.log']);

    case 'microcat'
        moor = loc;
        pd.rawpath = fullfile(mg.moordatadir, 'raw', mg.cruise, 'microcat');
        pd.infofile = fullfile(mg.moordatadir, 'proc', moor, [moor 'info.dat']);
        pd.stage1path = fullfile(mg.moordatadir, 'proc', moor, 'microcat');
        pd.stage1form = [moor '_%4.4d.raw'];
        pd.stage1log = fullfile(pd.stage1path,'stage1_log');
        pd.stage2path = fullfile(mg.moordatadir, 'proc', moor, 'microcat');
        pd.stage2form = [moor '_%4.4d.use'];
        pd.stage2log = fullfile(pd.stage2path, ['stage2_log_' moor,'.log']);
        pd.stage2figpath = fullfile(mg.reportdir, 'figs');
        pd.stage3path = fullfile(mg.moordatadir, 'proc', moor, 'microcat');
        pd.stage3form = [moor '_%0.3d.microcat'];

    case {'nor','nortek'}
        datatype = 'nor';
        moor = loc;
        pd.rawpath = fullfile(mg.moordatadir, 'raw', mg.cruise, datatype);
        if ~exist(pd.rawpath,'dir')
            pd.rawpath = fullfile(mg.moordatadir, 'raw', mg.cruise, 'nortek');
        end
        pd.infofile = fullfile(mg.moordatadir, 'proc', moor, [moor 'info.dat']);
        pd.listfile = fullfile(pd.rawpath, [moor '_filenames.txt']);
        pd.stage1path = fullfile(mg.moordatadir, 'proc', moor, datatype);
        pd.stage1log = fullfile(pd.stage1path, [moor '_Nortek_stage1.log']);
        pd.stage1form = [moor '_%d.raw'];
        pd.stage2path = fullfile(mg.moordatadir, 'proc', moor, datatype);
        pd.stage2log = fullfile(pd.stage2path, ['stage2_log_Nortek_' moor]);
        pd.stage2form = [moor '_%d.use'];
        pd.stage3path = fullfile(mg.moordatadir, 'proc', moor, datatype);
        pd.stage3form = [moor '_%d.edt'];
        pd.stage3log = [moor '_%d.edt.log'];
        pd.stage3formh = [moor '_%d.highedt'];
        pd.stage3logh = [moor '_%d.highedt.log'];   
        pd.stage3forml = [moor '_%d.lowedt'];
        pd.stage3logl = [moor '_%d.lowedt.log'];

    case 'bpr'
        moor = loc;
        pd.rawpath = fullfile(mg.moordatadir, 'raw', mg.cruise, 'seagauge');
        pd.rawform = [moor '_%4.4d_data.tid'];
        pd.infofile = fullfile(mg.moordatadir, 'proc', moor, [moor 'info.dat']);
        pd.offsetfile = fullfile(mg.moordatadir, 'raw', mg.cruise, 'clock_offset.dat');
        pd.stage1path = fullfile(mg.moordatadir, 'proc', moor, 'seagauge');
        pd.stage1log = fullfile(pd.stage1path, [moor '_seaguard_stage1.log']);
        pd.stage1form = [moor '_%5.5d.raw'];

    case {'adcp' 'adp'}
        disp(datatype)
        moor = loc;
        pd.rawpath = fullfile(mg.moordatadir, 'raw', mg.cruise, datatype);
        pd.infofile = fullfile(mg.moordatadir, 'proc', moor, [moor 'info.dat']);
        pd.stage1path = fullfile(mg.moordatadir, 'proc', moor, datatype);
        pd.stage1log = fullfile(pd.stage1path, [moor '_ADCP_stage1.log']);
        pd.stage1form = [moor '_%d_bin%02.f.raw'];
        pd.stage2inform = [moor '_%d_bin'];
        pd.stage2path = fullfile(mg.moordatadir, 'proc', moor, datatype);
        pd.stage2log = fullfile(pd.stage2path, [moor '_ADCP_stage2.log']);
        pd.stage2form = [moor '_%d_bin%02.f.use'];
        pd.stage3path = fullfile(mg.moordatadir, 'proc', moor, datatype);
        pd.stage3form = [moor '_%d_bin%02.f.edt'];
        pd.stage3log = [moor '_%d_bin%02.f.edt.log'];
        pd.stage3formh = [moor '_%d_bin%02.f.highedt'];
        pd.stage3logh = [moor '_%d_bin%02.f.highedt.log'];   
        pd.stage3forml = [moor '_%d_bin%02.f.lowedt'];
        pd.stage3logl = [moor '_%d_bin%02.f.lowedt.log'];

    case 'cal_coef' %***ctd path setting should be its own thing, same for microcat_cal_dip and for cal_coef?
        cast = loc;
        pd.datadir = mg.moordatadir;
        fp = fileparts(which('moor_setup'));
        pd.coef_dir = fullfile(fp,'metadata','cal_coef',lower(mg.project));
        pd.mc_dir = fullfile(pd.datadir, 'proc_calib', mg.cruise, 'cal_dip','microcat');
        pd.ctd1hz_file = fullfile(mg.ctddir,sprintf('ctd_%s_%3.3d_psal.nc',mg.cruise,cast));
        if ~exist(pd.ctd1hz_file,'file')
            pd.ctd1hz_file = fullfile(mg.ctddir,sprintf('ctd_%s_%3.3d_1hz.nc',mg.cruise,cast));
        end
        pd.ctdcnv_file = fullfile(mg.ctddir,'ASCII_FILES',sprintf('%s_CTD_%3.3d.cnv',upper(mg.cruise),cast));
        if ~exist(pd.ctdcnv_file,'file')
            pd.ctdcnv_file = fullfile(mg.ctddir,'ASCII_FILES',sprintf('%s_CTD%3.3d.cnv',upper(mg.cruise),cast));
            if ~exist(pd.ctdcnv_file,'file')
                pd.ctdcnv_file = fullfile(mg.ctddir,'ASCII_FILES',sprintf('%s_%3.3d.cnv',upper(mg.cruise),cast));
            end
        end
        pd.bottle_file = fullfile(mg.ctddir,'ASCII_FILES',sprintf('%s_CTD_%3.3d.ros',upper(mg.cruise),cast));
        if ~exist(pd.bottle_file,'file')
            pd.bottle_file = fullfile(mg.ctddir,'ASCII_FILES',sprintf('%s_CTD%3.3d.ros',upper(mg.cruise),cast));
            if ~exist(pd.bottle_file,'file')
                pd.bottle_file = fullfile(mg.ctddir,'ASCII_FILES',sprintf('%s_%3.3d.ros',upper(mg.cruise),cast));
            end
        end
        pd.mc_file = fullfile(pd.mc_dir,sprintf('cast%d',cast),sprintf('cast%d_',cast));
        pd.info_file = fullfile(pd.datadir,'proc_calib',mg.cruise,'cal_dip',sprintf('cast%dinfo.dat',cast));
        pd.ctdformat = 'mstar';
        pd.ctdcnv_cunit = 'S/m'; %***
        pd.ctd1hz_cunit = 'mS/cm';
        pd.mc_cunit = 'mS/cm';
        pd.mc_ext = '.raw';

    case 'mcgrid'
        pd.hydrodir = fullfile(mg.moordatadir,'proc','hydro_grid');
        pd.grdatdir = fullfile(mg.moordatadir,'proc','hydro_grid_merged');
        pd.figdir = fullfile(mg.moordatadir,'Figures');
        pd.mooringpath = fullfile(mg.moordatadir,'proc');

    case 'oceansites'
        moor = loc;
        pd.infofile = fullfile(mg.moordatadir,'proc',moor,[moor 'info.dat']);
        pd.ncpre = fullfile(mg.moordatadir,'oceansites_format',moor);


end

%now modify defaults if necessary
% (keep above defaults more readable)
switch datatype
    case 'cal_coef'
        defs = {'dy053' 'dy078' 'ar304' 'dy120' 'jc238' 'dy181'};
        if ~ismember(mg.cruise,defs)
        %cruise, ctdformat, ctdcnv_cunit, ctd1hz_cunit, mc_cunit, mc_ext
        a = {'kn221-02' 'aoml' 'S/m' 'S/m' 'S/m' '.raw';...
            'kn221-03' 'aoml' 'S/m' 'S/m' 'mS/cm' 'raw';...
            'pe399' 'mstar' 'mS/cm' 'mS/cm' 'mS/cm' '.raw';...
            };
        m = strcmp(mg.cruise,a(:,1));
        [pd.ctdformat,pd.ctdcnv_cunit,pd.ctd1hz_cunit,pd.mc_cunit,pd.mc_ext] = deal(a{m,2:end});
        end
        switch mg.cruise %some are different
            case 'kn221-02'
                pd.ctd1hz_file = fullfile(mg.ctddir,'..','CTD_FINAL','CTD','binavg_1Hz_mat','ctd_cal_kn211_alsta_1Hz.mat');
                pd.ctdcnv_file = fullfile(mg.ctddir,'..','CTD_FINAL','CTD','os1407_ctd_binavg_1Hz',sprintf('OS1407_%.3.3d.cnv',cast));
                pd.bottle_file = fullfile(mg.ctddir,'..','ctd_uncalib','datcnv',sprintf('OS1407_%3.3d.ros',cast));
            case 'kn221-03'
                pd.ctd1hz_file = fullfile(mg.ctddir,'..','1hz_bin_averaged_mat','ctd_cal_kn221-03_1Hz.mat');
                pd.ctdcnv_file = fullfile(mg.ctddir,'..','1hz_bin_averaged',sprintf('kn221-03_%3.3d.cnv',cast));
                pd.bottle_file = fullfile(mg.ctddir,'..','ros_files',sprintf('kn221-03_%3.3d.ros',cast));
            case 'pe400'
                pd.ctd1hz_file = fullfile(mg.ctddir,'..','1hz_bin_averaged_mat','PE400_1Hz.mat');
                pd.ctdcnv_file = fullfile(mg.ctddir,'..','1hz_bin_averaged',sprintf('PE400_%d.cnv',cast));
                pd.bottle_file = fullfile(mg.ctddir,'..','ros_files',sprintf('PE400_%d.ros',cast));
            case 'pe399'
                btldir = fullfile(pd.btldir,'CTD_Processed_by_Sven',sprintf('CTD-%2.2d',cast));
                dd=dir([btldir '/*.ros']);
                pd.bottle_file = fullfile(btldir,dd.name);%['PE399_',sprintf('%3.3d',calp.cast),'.ros']; % .ros
                pd.ctdcnv_file = fullfile(mg.ctddir,'ASCII_FILES',sprintf('PE399_%3.3d_align_ctm.cnv',cast));
            case 'dy078'
                pd.bottle_file = fullfile(mg.ctddir,'ASCII_FILES',sprintf('CTD%3.3d.ros',cast));
                pd.ctd_cnvfile = fullfile(mg.ctddir,'ASCII_FILES',sprintf('CTD%3.3d_align_actm.cnv',cast));
            case 'ar304'
                pd.ctd_1hzfile = fullfile(mg.ctddir,sprintf('ctd_%s_%3.3d_raw.nc',mg.cruise,cast));
                pd.bottle_file = fullfile(mg.ctddir,'ASCII_FILES',sprintf('ar30-04%3.3d.ros',cast));
                pd.ctd_cnvfile = fullfile(mg.ctddir,'ASCII_FILES',sprintf('ar30-04%3.3d.cnv',cast));
            case 'dy181'
                pd.bottle_file = fullfile(mg.ctddir,sprintf('fir_%s_%03d.nc',mg.cruise,cast));
        end
end
