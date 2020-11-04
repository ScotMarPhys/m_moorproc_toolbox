function m_create_file(ncfile)

% Create an empty mstar file and create 'header' data in NetCDF global attributes
m_common

% torg = datenum(MEXEC_G.MSTAR_TIME_ORIGIN); %default time origin

nc_create_empty(ncfile.name,'nc_clobber');

% % % % % nc_attput(ncfile.name,nc_global,'mstar_string','mstar_version_1.0'); %Always make the first 5 characters of this string identical to 'mstar'
% % % % % nc_attput(ncfile.name,nc_global,'openflag','W'); % set to W if file is open to write. Otherwise R.
% % % % % nc_attput(ncfile.name,nc_global,'date_file_updated',[0 0 0 0 0 0]); % This is the time of file update, stored as a recognisable 6 element vector
% % % % % nc_attput(ncfile.name,nc_global,'mstar_time_origin',MEXEC_G.MSTAR_TIME_ORIGIN); % This is the reference time for mstar data time origin and file update, stored as a recognisable 6 element vector
% % % % % nc_attput(ncfile.name,nc_global,'data_time_origin',MEXEC_G.MSTAR_TIME_ORIGIN); % This is the reference time for data, stored as a recognisable 6 element vector; usually time will be measured as decimal days or seconds since this time
% % % % % nc_attput(ncfile.name,nc_global,'time_convention','date_file_updated and data_time_origin are 6-element vectors, as commonly used in matlab date handling: [yyyy mo dd hh mm ss]');
% % % % % nc_attput(ncfile.name,nc_global,'version',0);
% % % % % nc_attput(ncfile.name,nc_global,'platform_type',' '); %eg 'ship'
% % % % % nc_attput(ncfile.name,nc_global,'platform_identifier',' '); %eg 'James_Cook'
% % % % % nc_attput(ncfile.name,nc_global,'platform_number',' '); %eg 'Cruise 31'
% % % % % nc_attput(ncfile.name,nc_global,'instrument_identifier','none_specified'); %eg 'CTD' or 'Current meter plus serial number'
% % % % % nc_attput(ncfile.name,nc_global,'recording_interval','none_specified'); %eg '1 Hz'
% % % % % nc_attput(ncfile.name,nc_global,'water_depth_metres',0); %eg 4000
% % % % % nc_attput(ncfile.name,nc_global,'instrument_depth_metres',0); %eg 3995; relevent for current meters
% % % % % nc_attput(ncfile.name,nc_global,'latitude',0); % decimal degrees; relevant for moorings or CTD stations
% % % % % nc_attput(ncfile.name,nc_global,'longitude',0); % decimal degrees; relevant for moorings or CTD stations
% % % % % nc_attput(ncfile.name,nc_global,'mstar_site',MEXEC_G.SITE); % identifier of computer where file was created
% % % % % % nc_attput(ncfile.name,nc_global,'data_time_origin',torg); % time origin. Will be used for converted pstar files, even if not used for anything else.
% % % % % comment_delimiter_string = MEXEC_G.COMMENT_DELIMITER_STRING;
% % % % % nc_attput(ncfile.name,nc_global,'comment_delimiter_string',comment_delimiter_string); % comment delimiter string
% % % % % nc_attput(ncfile.name,nc_global,'comment',comment_delimiter_string); % comments are free text containing any other useful information

% % hatt = m_default_attributes;
% % hatt_names = fieldnames(hatt);
% % 
% % for k = 1:length(hatt_names)
% %     attnam = hatt_names{k};
% %     cmd = ['attval = hatt.' attnam ';'];
% %     eval(cmd);
% %     nc_attput(ncfile.name,nc_global,attnam,attval);
% % end
hdef = m_default_attributes;

m_write_header(ncfile,hdef);

m_create_padvar(ncfile);

nc_padheader(ncfile.name,7752);

m_update_filedate(ncfile);

% m_create_data_time_origin(ncfile);

return

% % pstar header vars:


% % % blank6 = char(fread(fid,6,'*uchar'))';
% % % h.prefil = char(fread(fid,8,'*uchar'))';
% % % h.postfl = char(fread(fid,8,'*uchar'))';
% % % h.noflds = fread(fid,1,'uint32');
% % % h.norecs = fread(fid,1,'uint32');
% % % h.nrows = fread(fid,1,'uint32');
% % % h.nplane = fread(fid,1,'uint32');
% % % h.icent = fread(fid,1,'uint32');
% % % h.iymd = fread(fid,1,'uint32');
% % % h.ihms = fread(fid,1,'uint32');
% % % h.recint = char(fread(fid,16,'*uchar'))';
% % % 
% % % if lb == r5 %read them as real*5 and unpack; skip the MEXEC.status byte at read time. Then skip 2 blanks.
% % %     alat6 = char(fread(fid,5,'5*uchar',1))';
% % %     blank2 = char(fread(fid,2,'*uchar'))';
% % %     alat = punpack(alat6);
% % %     along6 = char(fread(fid,5,'5*uchar',1))';
% % %     blank2 = char(fread(fid,2,'5*uchar',1))';
% % %     along = punpack(along6);
% % %     dpthi6 = char(fread(fid,5,'5*uchar',1))';
% % %     blank2 = char(fread(fid,2,'*uchar'))';
% % %     dpthi = punpack(dpthi6);
% % %     dpthw6 = char(fread(fid,5,'5*uchar',1))';
% % %     blank2 = char(fread(fid,2,'*uchar'))';
% % %     dpthw = punpack(dpthw6);
% % % else %lb == r8; read them as real*8
% % %     alat = fread(fid,1,'double');
% % %     along = fread(fid,1,'double');
% % %     dpthi = fread(fid,1,'double');
% % %     dpthw = fread(fid,1,'double');
% % % end
% % % lonm = abs(60*(along-lond));
% % % disp(['Data Name :  ' h.datnam ' ' h.site ' ' h.vers]);
% % % disp(['Platform :   ' h.platyp ' ' h.platnam ' ' h.pltnum]);
% % % disp(['Instrument :' h.instmt '   dpthi ' sprintf('%8.2f',h.dpthi) '   dpthw ' sprintf('%8.2f',h.dpthw)]);
% % % disp(['Fields :    ' sprintf('%3d',h.noflds) '     Data Cycles ' sprintf('%8d',h.norecs) '     Rows ' sprintf('%4d',h.nrows) '   Planes ' sprintf('%4d',h.nplane)]);
% % % disp(['Position (lat lon) : '  sprintf('%10.5f',alat) '  ' sprintf('%10.5f',along)]);
% % % disp(['Position (lat lon) : '  sprintf('%4d %06.3f',latd,latm) ' ' sprintf('%4d %06.3f',lond,lonm)]);
% % % disp(['Time origin : ' sprintf('%02d',h.icent/100) '/' sprintf('%06d',h.iymd) '/' sprintf('%06d',h.ihms)]);
% % % disp('************************************************************************');
% % % for k = 1:h.noflds
% % %     disp(['*' sprintf('%3d',k) '*' h.fldnam{k} '*' h.fldunt{k} '* ' sprintf('%15.3f',h.alrlim(k)) ' * ' sprintf('%15.3f',h.uprlim(k)) ' * ' sprintf('%10.3f',h.absent(k)) ' *']);
% % % end
% % % disp('************************************************************************');
% % % for k = 1:12
% % %     thiscoment = h.coment(k,:);
% % %     blank72 = '                                                                        ';
% % %     if strcmp(blank72,thiscoment) == 0
% % %         disp(thiscoment);
% % %     end
% % % end
% % % 
% % % 
% % % fclose(fid);
