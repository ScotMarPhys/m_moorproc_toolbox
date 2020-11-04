function m_copy_variable(ncfile_in,vname,ncfile_ot,newname,indexrows,indexcols)

% Copy a variable from input file to output file, with optional new name in
% output file
% 
% Optional copy of subset of vars
% 
% if it is recognised as a time variable, adjust for data_time_origin

m_common
global toco
clear v

v.name = newname;
vdata = nc_varget(ncfile_in.name,vname);

% % % % %  % no time adjustment in quick copy
% % % % % hin = m_read_header(ncfile_in);
% % % % % hot = m_read_header(ncfile_ot);
% % % % % torg1 = hin.data_time_origin;
% % % % % torg2 = hot.data_time_origin;
% % % % % tdif = torg1-torg2;
% % % % % if max(abs(tdif)) > 0
% % % % %     if m_isvartime(vname)
% % % % %         % this is a time variable name; adjust for data time origin
% % % % %         vdata = m_adjtime(vname,vdata,hin,hot);
% % % % %     end
% % % % % end

if nargin == 4
    v.data = vdata;
else
    v.data = vdata(indexrows,indexcols);
end
toc-toco
disp('c4'); toco=toc;
m_write_variableq(ncfile_ot,v,'nodata'); %write the variable information into the header but not the data
toc-toco
disp('c5'); toco=toc;

% next copy the attributes
vinfo = nc_getvarinfo(ncfile_in.name,vname);
va = vinfo.Attribute;
toc-toco
disp('c6'); toco=toc;

for k2 = 1:length(va)
    vanam = va(k2).Name;
    vaval = va(k2).Value;
    nc_attputq(ncfile_ot.name,newname,vanam,vaval);
    toc-toco
disp('c7'); toco=toc;
num2str(k2)
end

% now write the data, using the attributes already saved in the output file
% this provides the opportunity to change attributes if required, eg fillvalue
toc-toco
disp('c1'); toco = toc;
nc_varput(ncfile_ot.name,newname,v.data);
toc-toco
disp('c2'); toco = toc;
% m_uprlwr(ncfile_ot,newname); % not strictly needed if straight copy
m_uprlwr(ncfile_ot,newname,v.data); % not strictly needed if straight copy
toc-toco
disp('c3'); toco = toc;

return