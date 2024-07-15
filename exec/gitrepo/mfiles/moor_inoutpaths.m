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
            castn = str2num(cast);
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
        pd.stage2log = fullfile(mg.reportdir, 'stats', ['stage2_log_' moor]);
        pd.stage2figpath = fullfile(mg.reportdir, 'figs');

    case {'nor'}
        pd.rawpath = fullfile(mg.moordatadir, 'raw', mg.cruise, datatype);
        pd.infofile = fullfile(mg.moordatadir, 'proc', moor, [moor 'info.dat']);
        pd.listfile = fullfile(pd.rawpath, [moor '_filenames.txt']);
        pd.stage1path = fullfile(mg.moordatadir, 'proc', moor, datatype);
        pd.stage1log = fullfile(pd.stage1path, [moor '_Nortek_stage1.log']);
        pd.stage1form = [moor '_%d.raw'];
        pd.stage2path = fullfile(mg.moordatadir, 'proc', moor, datatype);
        pd.stage2log = fullfile(mg.reportdir, 'stats', ['stage2_log_Nortek_' moor]);
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

    case 'adcp'
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
        pd.datadir = mg.moordatadir;
        fp = fileparts(which('moor_setup'));
        pd.coef_dir = fullfile(fp,'metadata','cal_coef',lower(mg.project));
        pd.mc_dir = fullfile(pd.datadir, 'proc_calib', mg.cruise, 'cal_dip','microcat');
        pd.info_dir = fullfile(pd.datadir, 'proc_calib', mg.cruise, 'cal_dip');
        pd.ctd_dir = fullfile(mg.ctddir);
        pd.ctdraw_dir = fullfile(mg.ctddir,'ASCII_FILES'); %.cnv location
        pd.btldir = fullfile(mg.ctddir,'ASCII_FILES');
        pd.ctdfile = sprintf('ctd_%s_%3.3d_psal.nc');
        pd.ctd_cnvfile = sprintf('%s_CTD_%3.3d.cnv',upper(mg.cruise),loc);
        pd.bottle_file = sprintf('%s_CTD_%3.3d.ros',upper(mg.cruise),loc);
        pd.mc_file = fullfile(sprintf('cast%d',loc),sprintf('cast%d_',loc));
        pd.info_file = sprintf('cast%dinfo.dat',loc);

end

%now modify defaults if necessary (keep above defaults more readable)
if mg.YEAR<2022

    switch datatype
        
        case 'cal_coef'
        switch mg.cruise %some are different
            case 'kn221-02'
                pd.ctd_dir = fullfile(mg.ctddir,'..','CTD_FINAL','CTD','binavg_1Hz_mat');
                pd.ctdraw_dir = fullfile(mg.ctddir,'..','CTD_FINAL','CTD','os1407_ctd_binavg_1Hz');
                pd.btldir = fullfile(mg.ctddir,'..','ctd_uncalib','datcnv');
                pd.ctdfile = 'ctd_cal_kn211_alsta_1Hz.mat';
                pd.bottle_file = sprintf('OS1407_%3.3d.ros',loc);
                pd.ctd_cnvfile = sprintf('OS1407_%3.3d.cnv',loc);
            case 'kn221-03'
                pd.ctd_dir = fullfile(mg.ctddir,'..','1hz_bin_averaged_mat');
                pd.ctdraw_dir = fullfile(mg.ctddir,'..','1hz_bin_averaged');
                pd.btldir = fullfile(mg.ctddir,'..','ros_files');
                pd.ctd_file    = 'ctd_cal_kn221-03_1Hz.mat' ; % ['ab1104_ctd_caldip_1Hz.mat'];
                pd.bottle_file = ['kn221-03_',sprintf('%3.3d',calp.cast),'.ros'];
                pd.ctd_cnvfile = ['kn221-03_',sprintf('%3.3d',calp.cast),'.cnv'];
            case 'pe400'
                pd.ctd_dir = fullfile(mg.ctddir,'..','1hz_bin_averaged_mat');
                pd.ctdraw_dir = fullfile(mg.ctddir,'..','1hz_bin_averaged');
                pd.btldir = fullfile(mg.ctddir,'..','ros_files');
                pd.ctd_file    = 'PE400_1Hz.mat' ; % ['ab1104_ctd_caldip_1Hz.mat'];
                pd.bottle_file = ['PE400_',num2str(calp.cast),'.ros'];
                pd.ctd_cnvfile = ['PE400_',num2str(calp.cast),'.cnv'];
            case 'pe399'
                pd.btldir = fullfile(pd.btldir,'CTD_Processed_by_Sven',sprintf('CTD-%2.2d',loc));
                pd.ctd_file    = ['ctd_pe399_',sprintf('%3.3d',calp.cast),'_1hz.nc'];
                dd=dir([pd.btldir '/*.ros']);
                pd.bottle_file = dd.name;%['PE399_',sprintf('%3.3d',calp.cast),'.ros']; % .ros
                pd.ctd_cnvfile = ['PE399_',sprintf('%3.3d',calp.cast),'_align_ctm.cnv'];
            case 'dy053'
                pd.ctd_file = sprintf('ctd_%s_%3.3d_1hz.nc',mg.cruise,loc);
                pd.bottle_file = sprintf('%s_%3.3d.ros',upper(mg.cruise),loc);
                pd.ctd_cnvfile = sprintf('%s_%3.3d.cnv',upper(mg.cruise),loc);
            case 'dy078'
                pd.ctd_file = sprintf('ctd_%s_%3.3d_1hz.nc',mg.cruise,loc);
                pd.bottle_file = sprintf('CTD%3.3d.ros',loc);
                pd.ctd_cnvfile = sprintf('CTD%3.3d_align_actm.cnv',loc);
            case 'ar304'
                pd.ctd_file = sprintf('ctd_%s_%3.3d_raw.nc',mg.cruise,loc);
                pd.bottle_file = sprintf('ar30-04%3.3d.ros',loc);
                pd.ctd_cnvfile = sprintf('ar30-04%3.3d.cnv',loc);
            case 'dy120'
                pd.bottle_file = sprintf('%s_%3.3d.ros',upper(mg.cruise),loc);
                pd.ctd_cnvfile = sprintf('%s_%3.3d.cnv',upper(mg.cruise),loc);
        end

    end

end