function infile = mc_raw_filenames(sn, pd)
%try different possible filename forms until we find the raw microcat data
%file (for caldip or deployment files)

infiles = {
    sprintf('%4.4dcal2.asc',sn);
    sprintf('%4.4dcal.asc',sn);
    sprintf('%3.3dcal.asc',sn);
    sprintf('%4.4dCAL.asc',sn);
    sprintf('%3.3dCAL.asc',sn);
    sprintf('cal%4.4d.asc',sn)
    sprintf('%4.4d_cal_dip2.asc',sn)
    sprintf('%4.4d_cal_dip_data2.asc',sn)
    sprintf('%4.4d_test.asc',sn)
    sprintf('%4.4d_cal_dip.asc',sn)
    sprintf('%4.4d_cal_dip_data.asc',sn)
    sprintf('%4.4d__cal_dip_data.asc',sn)
    sprintf('%4.4d_cal_dip_data.cnv',sn)
    sprintf('%4.4d_Cal_Dip_Data.cnv',sn)
    sprintf('%4.4d%srec2.asc',sn)
    sprintf('%4.4d%rec.asc',sn)
    sprintf('%3.3d%rec.asc',sn)
    sprintf('%4.4d%REC.asc',sn)
    sprintf('%3.3d%REC.asc',sn)
    fullfile('data',sprintf('%4.4d_2.asc',sn))
    fullfile('data',sprintf('%4.4d.asc',sn))
    sprintf('%4.4d%s',sn,'.asc')
    sprintf('%4.4d_DATA.asc',sn)
    sprintf('%4.4d_data.cnv',sn)
    sprintf('%4.4d_Data.cnv',sn)
    sprintf('%4.4d_data.asc',sn)
    };

for n = 1:length(infiles)
    infile = fullfile(pd.rawpath,infiles{n});
    if exist(infile,'file')
        datfileinfo = dir(infile);
        if datfileinfo.bytes>0
            %found it
            break
        end
    end
end
infile = fullfile(pd.rawpath,infiles{n});
