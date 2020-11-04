function pstar_to_mstar

% load pstar file into mstar file


m_common
m_margslocal
m_varargs

MEXEC_A.Mprog = 'pstar_to_mstar';
m_proghd;



pstar_fn = m_getpstarfilename;
mstar_fn = m_getfilename;

ncfile.name = mstar_fn;
ncfile = m_openot(ncfile); %check it is not an open mstar file


[d h] = pload(pstar_fn,'');



h2.dataname = m_remove_spaces(h.datnam); %pstar datanames are padded to length 8 with spaces
h2.version = m_pstarvers_to_mstarvers(h.vers);
h2.openflag = h.opwrit;
h2.platform_type = m_remove_outside_spaces(h.platyp);
h2.platform_identifier =  m_remove_outside_spaces(h.platnam);
h2.platform_number =  m_remove_outside_spaces(h.pltnum);
h2.instrument_identifier =  m_remove_outside_spaces(h.instmt);
h2.recording_interval =  m_remove_outside_spaces(h.recint);
h2.latitude = h.alat;
h2.longitude = h.along;
h2.water_depth_metres = h.dpthw;
h2.instrument_depth_metres = h.dpthi;
h2.mstar_site = ['pexec_' h.site];

cccc = h.icent;
ymd = h.iymd+0.1; %add 0.1 to be sure floor/rounding works as expected.
yy = floor(ymd/10000);
modd = ymd-10000*yy;
mo = floor(modd/100);
dd = floor(modd-100*mo);
yyyy = cccc+yy;

hms = h.ihms+0.1;
hh = floor(hms/10000);
mmss = hms - 10000*hh;
mm = floor(mmss/100);
ss = floor(mmss-100*mm);

% h2.data_time_origin = datenum(yyyy,mo,dd,hh,mm,ss)-datenum(1950,1,1,0,0,0);



ncfile.name = mstar_fn;

m_create_file(ncfile);

% mod by BAK 14 apr 2015, after DAS found some di279 data that had
% all-blank fields which crashed the nc_attput. The m_create_file sets up
% valid default values, whcih do not need to be overwritten if there are no
% genuine values to replace them.
if ~isempty(h2.dataname); nc_attput(ncfile.name,nc_global,'dataname',h2.dataname); end
nc_attput(ncfile.name,nc_global,'version',h2.version);
if ~isempty(h2.platform_type); nc_attput(ncfile.name,nc_global,'platform_type',h2.platform_type); end
if ~isempty(h2.platform_identifier); nc_attput(ncfile.name,nc_global,'platform_identifier',h2.platform_identifier); end
if ~isempty(h2.platform_number); nc_attput(ncfile.name,nc_global,'platform_number',h2.platform_number); end
if ~isempty(h2.instrument_identifier); nc_attput(ncfile.name,nc_global,'instrument_identifier',h2.instrument_identifier); end
if ~isempty(h2.recording_interval); nc_attput(ncfile.name,nc_global,'recording_interval',h2.recording_interval); end
nc_attput(ncfile.name,nc_global,'water_depth_metres',h2.water_depth_metres);
nc_attput(ncfile.name,nc_global,'instrument_depth_metres',h2.instrument_depth_metres);
nc_attput(ncfile.name,nc_global,'latitude',h2.latitude);
nc_attput(ncfile.name,nc_global,'longitude',h2.longitude);
m_set_data_time_origin(ncfile,yyyy,mo,dd,hh,mm,ss)

% nc_varput(ncfile.name,'data_time_origin',h2.data_time_origin); m_uprlwr(ncfile,'data_time_origin');
nc_attput(ncfile.name,nc_global,'mstar_site',h2.mstar_site);

for k = 1:size(h.coment,1)
    com = h.coment(k,:);
    com(strfind(com,char(0))) = ' ';
    if length(strfind(com,' ')) == 72; continue; end
    m_add_comment(ncfile,com);
end

%Now add variables

pst_fieldnames = fieldnames(d);
for k = 1:length(pst_fieldnames) %number of pstar variables
    clear v;
    v.name = pst_fieldnames{k};
    v.units = m_remove_outside_spaces(h.fldunt{k}); % pstar units are padded with spaces
    eval(['v.data = d.' pst_fieldnames{k} ';']);
    v.type = 'double';    
    disp(['writing mstar variable ' v.name]);
    m_write_variable(ncfile,v);
end

% m_finis(ncfile,0);

% m_verson(ncfile,version_increment); %advance the version

m_pstar_filedate(ncfile,pstar_fn); % set the file update variable; use unix date of pstar file being loaded
nowstring = datestr(now,31);
m_add_comment(ncfile,' ');
m_add_comment(ncfile,'This mstar file created from pstar file');
m_add_comment(ncfile,pstar_fn);
m_add_comment(ncfile,['at ' nowstring]);
nc_attput(ncfile.name,nc_global,'openflag','R'); %set the open/writing attribute

h = m_read_header(ncfile);
m_print_header(h)

hist = h;
hist.filename = ncfile.name;
MEXEC_A.Mhistory_ot{1} = hist;
% fake the input file details so that write_history works
histin = h;
histin.filename = pstar_fn;
histin.dataname = [];
histin.version = [];
histin.mstar_site = [];
MEXEC_A.Mhistory_in{1} = histin;
m_write_history;

return

% %
% %
% % h =
% %
% %      opwrit: 'R'
% %      rawdat: 'P'
% %      pipefl: ' '
% %      archiv: 'N'
% %       magic: 1.3450e+09
% %      prefil: '        '
% %      postfl: '        '
% %      noflds: 14
% %      norecs: 384
% %       nrows: 0
% %      nplane: 0
% %       icent: 2000
% %        iymd: 61208
% %        ihms: 85235
% %      fldnam: {1x128 cell}
% %      fldunt: {1x128 cell}
% %      alrlim: [1x128 double]
% %      uprlim: [1x128 double]
% %      absent: [1x128 double]
% %        site: 'jc'
% %
% %
% % nc_attput(ncfile.name,nc_global,'mstar_string','mstar_version_1.0'); %Always make the first 5 characters of this string identical to 'mstar'
% % nc_attput(ncfile.name,nc_global,'openflag','W'); %set the open/writing attribute immediately
% % nc_attput(ncfile.name,nc_global,'dataname',' ');
% % nc_attput(ncfile.name,nc_global,'version',0);
% % nc_attput(ncfile.name,nc_global,'platform_type',' '); %eg 'ship'
% % nc_attput(ncfile.name,nc_global,'platform_identifier',' '); %eg 'James_Cook'
% % nc_attput(ncfile.name,nc_global,'platform_identifier',' '); %eg 'Cruise 31'
% % nc_attput(ncfile.name,nc_global,'instrument_identifier',' '); %eg 'CTD' or 'Current meter plus serial number'
% % nc_attput(ncfile.name,nc_global,'recording_interval',' '); %eg '1 Hz'
% % nc_attput(ncfile.name,nc_global,'water_depth_metres',0); %eg 4000
% % nc_attput(ncfile.name,nc_global,'instrument_depth_metres',0); %eg 3995; relevent for current meters
% % nc_attput(ncfile.name,nc_global,'latitude',0); % decimal degrees; relevant for moorings or CTD stations
% % nc_attput(ncfile.name,nc_global,'longitude',0); % decimal degrees; relevant for moorings or CTD stations
% % nc_attput(ncfile.name,nc_global,'mstar_site',' '); % identifier of computer where file was created
% % comment_delimiter_string = ' \n ';
% % nc_attput(ncfile.name,nc_global,'comment_delimiter_string',comment_delimiter_string); % comment delimiter string
% % nc_attput(ncfile.name,nc_global,'comment',comment_delimiter_string); % comments are free text containing any other useful information
% % % pstar also has
% % % time origin
% %
% % % all files should have dimensions nrows and ncols, but we can't
% % % pre-define these because we don't know the number of rows and cols, and
% % % dimensions can't have their value changed (at least I think not...)
% %
% % m_update_filedate(ncfile);
% %
% % nc_padheader(ncfile.name,7752);
