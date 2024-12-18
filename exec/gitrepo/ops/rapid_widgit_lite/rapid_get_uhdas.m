% function rapid_get_uhdas

ncf_nb = '/mnt/uhdas_data/JC192/proc/os150nb/contour/os150nb.nc';
ncf_bb = '/mnt/uhdas_data/JC192/proc/os150bb/contour/os150bb.nc';

fl_nb = dir(ncf_nb);
fl_bb = dir(ncf_bb);

if fl_bb.datenum > fl_nb.datenum
	ncf = ncf_bb;
else
	ncf = ncf_nb;
end

ncf_inf = ncinfo(ncf)

nt = ncf_inf.Dimensions(1).Length;
nz = ncf_inf.Dimensions(2).Length;
nv = 10;

tme = ncread(ncf,'time',nt-nv,nv+1);
zz = ncread(ncf,'depth',[1 nt-nv],[nz nv+1]);
uv = ncread(ncf,'u',[1 nt-nv],[nz nv+1]);
vv = ncread(ncf,'v',[1 nt-nv],[nz nv+1]);

uM = ncreadatt(ncf,'u','missing_value');
vM = ncreadatt(ncf,'v','missing_value');

uv(uv == uM) = NaN;
vv(vv == vM) = NaN;