function m_print_header(h)

% return

latd = fix(h.latitude);
latm = abs(60*(h.latitude-latd));
lond = fix(h.longitude);
lonm = abs(60*(h.longitude-lond));
disp(['Data Name :  ' h.dataname ' <version> ' sprintf('%d',h.version) ' <site> ' h.mstar_site]);
disp(['Platform :   ' h.platform_type ' | ' h.platform_identifier ' | ' h.platform_number]);
disp(['Instrument : ' h.instrument_identifier '   dpthi ' sprintf('%8.2f',h.instrument_depth_metres) '   dpthw ' sprintf('%8.2f',h.water_depth_metres)]);
% disp(['Fields :    ' sprintf('%3d',h.noflds) '     Data Cycles ' sprintf('%8d',h.norecs) '     Rows ' sprintf('%4d',h.nrows) '   Planes ' sprintf('%4d',h.nplane)]);
disp(['Position (lat lon) : '  sprintf('%10.5f',h.latitude) '  ' sprintf('%10.5f',h.longitude)]);
disp(['Position (lat lon) : '  sprintf('%4d %06.3f',latd,latm) ' ' sprintf('%4d %06.3f',lond,lonm)]);
% disp(['Time origin : ' sprintf('%02d',h.icent/100) '/' sprintf('%06d',h.iymd) '/' sprintf('%06d',h.ihms)]);
disp(['Data time origin : ' h.data_time_origin_string]);
disp(['Fields :    ' sprintf('%3d',h.noflds)]);% '     Data Cycles ' sprintf('%8d',h.norecs) '     Rows ' sprintf('%4d',h.nrows) '   Planes ' sprintf('%4d',h.nplane)]);
disp('Dimension sets:');
disp('set  nrows      ncols      norecs');

for k = 1:h.numdimsets
    rn = h.rowname{k};
    suffix = rn(6:end);
    disp([sprintf('%-4s %-10d %-10d %-10d',[suffix ':'],h.rowlength(k),h.collength(k),h.rowlength(k)*h.collength(k))])
end
disp('************************************************************************************');
disp('*   *name      *units   *dims*      min     *      max     *     nabs * absval     *');
disp('************************************************************************************');
for k = 1:h.noflds
    lenfn = 10;    
    fn = h.fldnam{k};
    if length(fn) > lenfn
        fn = [fn(1:lenfn-1) '@'];
    end
       lenfu = 8;    
    fu = h.fldunt{k};
    if length(fu) > lenfu
        fu = [fu(1:lenfu-1) '@'];
    end

    lwrform = '%12.3f';
    if abs(h.alrlim(k)) < 0.1; lwrform = '%12.3e'; end
    if abs(h.alrlim(k)) == 0; lwrform = '%10.1f  '; end
    if abs(h.alrlim(k)) > 9999999.9; lwrform = '%12.3e'; end
    uprform = '%12.3f';
    if abs(h.uprlim(k)) < 0.1; uprform = '%12.3e'; end
    if abs(h.uprlim(k)) == 0; uprform = '%10.1f  '; end
    if abs(h.uprlim(k)) > 9999999.9; uprform = '%12.3e'; end
    
%     disp(['*' sprintf('%3d',k) '*' sprintf('%-10s',fn) '*' sprintf('%-8s',fu) '*' sprintf('%8d',h.dimrows(k)) '*' sprintf('%8d',h.dimcols(k)) '* ' sprintf(lwrform,h.alrlim(k)) ' * ' sprintf(uprform,h.uprlim(k)) ' * ' sprintf('%6d',h.num_absent(k)) ' * ' sprintf('%10.3f',h.absent(k)) ' *']);
    disp(['*' sprintf('%3d',k) '*' sprintf('%-10s',fn) '*' sprintf('%-8s',fu) '*' sprintf('%3s ',h.dimsset{k}) '* ' sprintf(lwrform,h.alrlim(k)) ' * ' sprintf(uprform,h.uprlim(k)) ' *' sprintf('%9d',h.num_absent(k)) ' * ' sprintf('%10.3f',h.absent(k)) ' *']);
end
disp('************************************************************************************');


c = h.comment;
delim = h.comment_delimiter_string;
h.comment_delimiter_string = delim;
delimindex = strfind(c,delim);
ncoms = length(delimindex);
for k = 2:ncoms %if there are no genuine comments ncoms will be 1 and this loop won't be executed
    disp(['comment: ' sprintf('%s',c(delimindex(k-1)+length(delim):delimindex(k)-1))]);
end



disp(['File last updated : ' h.last_update_string ]);
disp(' ');




return
