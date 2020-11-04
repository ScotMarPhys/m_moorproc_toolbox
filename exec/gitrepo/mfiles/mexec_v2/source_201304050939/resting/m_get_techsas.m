function ncfile = m_get_techsas(techsas_in,ncfile,dataname)

MEXEC_A.Mprog = 'm_get_techsas';
m_proghd;


% % techsas_in.name = '/Users/bak/Documents/aafrom_macbookpro_for_jamaica/rvs_nc_data_2/20080629-090012-gppat-GPPAT.att';
% % ncfile.name = 'techsas_out';
% % dataname = 'techsas_test';


m_global

ncfile.name = m_add_nc(ncfile.name);

ncfile = m_openot(ncfile);
nc_attput(ncfile.name,nc_global,'dataname',dataname); %set the dataname

nc_attput(ncfile.name,nc_global,'platform_type',MEXEC_G.PLATFORM_TYPE); %eg 'ship'
nc_attput(ncfile.name,nc_global,'platform_identifier',MEXEC_G.PLATFORM_IDENTIFIER); %eg 'James_Cook'
nc_attput(ncfile.name,nc_global,'platform_number',MEXEC_G.PLATFORM_NUMBER); %eg 'Cruise 31'



techsas_names = m_unpack_varnames(techsas_in);

for k = 1:length(techsas_names)
    clear techsas_data techsas_units v
    techsas_data = nc_varget(techsas_in.name,techsas_names{k});
    techsas_units = nc_attget(techsas_in.name,techsas_names{k},'units');
    v.name = techsas_names{k}; v.data = techsas_data; v.units = techsas_units;
    m_write_variable(ncfile,v);
end

m_set_data_time_origin(ncfile,1899,12,30,0,0,0)

nowstring = datestr(now,31);
m_add_comment(ncfile,'This mstar file created from techsas file');
m_add_comment(ncfile,techsas_in.name);
m_add_comment(ncfile,['at ' nowstring]);


m_finis(ncfile);

h = m_read_header(ncfile);
m_print_header(h);
return