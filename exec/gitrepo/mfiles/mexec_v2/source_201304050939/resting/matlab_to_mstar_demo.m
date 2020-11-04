% testdata;
torg = datenum(1950,1,1,0,0,0)

t = [1:100]+now-torg;
t = [1:100];
t = t(:)';
t2 = repmat(t,[10 1]);

p = 1:10;
p = p(:);
p = repmat(p,[1,100]);

pt = t;
pt = p(:,1)*t;

u = 2*pt;
v= pt+cos(pt);


% end of making test data

MEXEC_A.Mprog = 'mat2mstar demo';

m_proghd;

ncfile.name = './f1';
dataname = 'bak_test_mat2mstar_demo';

ncfile = m_openot(ncfile);
nc_attput(ncfile.name,nc_global,'dataname',dataname); %set the dataname


var.name = 'time';
var.data = t;
var.units = 'days';
var.type = 'double';


% % % need to consider possibility that if writing to an existign variable, units aren't needed
% % % need to consider possibility that if writing to an existign variable, units aren't needed
% % % need to consider possibility that if writing to an existign variable, units aren't needed
% % % need to consider possibility that if writing to an existign variable, units aren't needed
% % % need to consider possibility that if writing to an existign variable, units aren't needed
% % % need to consider possibility that if writing to an existign variable, units aren't needed
% % % need to consider possibility that if writing to an existign variable, units aren't needed
% % % need to consider possibility that if writing to an existign variable, units aren't needed
% % % need to consider possibility that if writing to an existign variable, units aren't needed
% % % need to consider possibility that if writing to an existign variable, units aren't needed

%maybe need some routines like

% % % put_units
% % % put_data (check dimensions OK)
% % % get_data

m_write_variable(ncfile,var);

var4.name = 'p';
var4.data = p;
var4.units = 'metres';
% var2.type = 'double';   %use default type which is double
m_write_variable(ncfile,var4);


var2.name = 'u';
var2.data = u;
var2.units = 'cm/s';
% var2.type = 'double';   %use default type which is double
m_write_variable(ncfile,var2);



var3.name = 'v';
var3.data = v;
var3.units = 'furlongs';
% var2.type = 'double';   %use default type which is double
m_write_variable(ncfile,var3);


% m_write_variable(ncfile,t,'t','days','double');
% m_write_variable(ncfile,u,'u','cm/s'); %type 'double' is default
% m_write_variable(ncfile,v,'v','m/s');

m_add_comment(ncfile,'test mstar file');
m_add_comment(ncfile,'next comment');
m_add_comment(ncfile,'');
m_add_comment(ncfile,'last was null');



m_finis(ncfile);


metadata = nc_info(ncfile.name); %refresh metadata
ncfile.metadata = metadata;

h = m_read_header(ncfile);
m_print_header(h);
