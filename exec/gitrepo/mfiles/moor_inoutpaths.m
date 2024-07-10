function pd = moor_inoutpaths(datatype,loc)
% predirs = moor_inoutpaths(datatype,loc)
%
% predirs = moor_inoutpaths('microcat_cal_dip',cast)
% predirs = moor_inoutpaths('microcat',moor)
% predirs = moor_inoutpaths('nor',moor)
% predirs = moor_inoutpaths('bpr',moor)
% predirs = moor_inoutpaths('adcp',moor)
%
% set input and output directories and filename prefixes to be accessed by
% successive stages of processing for each type of data given by datatype: 
%   microcat_cal_dip (2nd argument is the CTD cast number as a scalar or a string, e.g. 5 or '5')
%   microcat or nor (2nd argument is the mooring e.g. 'ebh3_15_2022')

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
    case {'nor','nortek'}
        moor = loc;
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

end