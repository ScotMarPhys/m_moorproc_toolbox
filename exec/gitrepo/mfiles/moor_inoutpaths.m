function pd = moor_inoutpaths(datatype,plat)
% predirs = moor_inoutpaths(datatype,plat)
%
% predirs = moor_inoutpaths('microcat_cal_dip',cast)
% predirs = moor_inoutpaths('microcat',moor)
%
% set input and output directories and filename prefixes to be accessed by
% successive stages of processing for each type of data given by datatype: 
%   microcat_cal_dip (2nd argument is the CTD cast number as a scalar or string, e.g. 1 or '1')
%   microcat or nor (2nd argument is the mooring e.g. 'ebh3_15_2022')

global MOORPROC_G
mg = MOORPROC_G;

if isnumeric(plat)
    plat = num2str(plat);
end

switch datatype
    case 'microcat_cal_dip'
        cast = plat;
        pd.infofile = fullfile(mg.moordatadir, 'raw', mg.cruise, 'microcat_cal_dip', ['cast' cast 'info.dat']);
        pd.rawpath = fullfile(mg.moordatadir, 'raw', mg.cruise, 'microcat_cal_dip', ['cast' cast]);
        pd.stage1path = fullfile(mg.moordatadir, 'proc_calib', mg.cruise, 'cal_dip', 'microcat', ['cast' cast]);

        stage2path = fullfile(mg.moordatadir, 'proc_calib', mg.)