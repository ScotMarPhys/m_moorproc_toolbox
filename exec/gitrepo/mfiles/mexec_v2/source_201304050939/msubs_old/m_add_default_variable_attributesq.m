function m_add_default_variable_attributes(ncfile,var_name)

% Add the default attributes for a variable

vinfo = nc_getvarinfo(ncfile.name,var_name);
global toco
toc-toco
disp('dc1'); toco = toc;


if vinfo.Nctype == 6 %this is a double variable
    attlist = {
        'long_name'         ' '
        'units'             ' '
        'min_value'          -99999
        'max_value'          -99999
        '_FillValue'         -99999
        'missing_value'      -99999
        'number_fillvalue'   0
        };

end
toc-toco
disp('dc2'); toco = toc;

if vinfo.Nctype == 2 %this is a char variable
    attlist = {
        'long_name'         ' '
        'units'             ' '
        '_FillValue'        ' '
        'missing_value'     ' '
        'number_fillvalue'   0
        };

end
toc-toco
disp('dc3'); toco = toc;

for k = 1:size(attlist,1)
    nc_attput(ncfile.name,var_name,attlist{k,1},attlist{k,2});
    toc-toco
disp('dc5'); num2str(k)
!ls -ld gash2.nc
toco = toc;
end
toc-toco
disp('dc4'); toco = toc;
