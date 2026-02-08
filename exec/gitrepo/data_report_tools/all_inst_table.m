function id_z_sn = all_inst_table(id, z, sn, varargin)
% function id_z_sn = all_inst_table(id, z, sn, varargin)
%
% search for different types of instrument by id number, then add
% subdirectories and abbreviations as well as variables for reports
% 
% output is table that can be looped through, or have rows extracted
%
% z and sn are preserved but not parsed; therefore, to simply look up the
% instrument corresponding to a sensor id code, do 
% id_z_sn = all_inst_table(id, -99, -99); 

%initialise output table
id = id(:); z = z(:); sn = sn(:);
id_z_sn = table(id,z,sn);
id_z_sn.inst = repmat({' '},size(id));
id_z_sn.instl = id_z_sn.inst;

% get lookup table between id and instrument type. each column is a cell
% array
l = inst_id_lookup;

%loop through rows of l and populate corresponding lines in output table
%using values in l.id
for ino = 1:length(l.inst)
    ii = find(ismember(id, l.id{ino}));
    id_z_sn.inst(ii) = l.inst(ino);
    id_z_sn.instl(ii) = l.instl(ino);
    id_z_sn.vars_st(ii) = l.vars_st(ino);
    id_z_sn.vars_cs(ii) = l.vars_cs(ino);
end

%subdirectory (is this used? it is, but should be replaced with moor_inoutpaths, as should suffix)***
id_z_sn.dirs = lower(id_z_sn.inst);
id_z_sn.dirs(strcmp(id_z_sn.inst,'MC')) = {'microcat'};
id_z_sn.dirs(strcmp(id_z_sn.inst,'ODOMC')) = {'microcat'};
id_z_sn.dirs(strcmp(id_z_sn.inst,'RCM11')) = {'rcm'};
%filename suffix
id_z_sn.suf = repmat({''},length(id),1);
id_z_sn.suf(strcmp(id_z_sn.inst,'DVS')) = {'_bin2'};
