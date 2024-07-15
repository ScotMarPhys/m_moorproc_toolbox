function id_z_sn = all_inst_table(id, z, sn)
% function id_z_sn = all_inst_table(id, z, sn)
%
% search for different types of instrument by id number, then add
% subdirectories and abbreviations as well as variables for reports
% 
% output is table that can be looped through, or have rows extracted
%
% z and sn are preserved but not parsed; therefore, to simply look up the
% instrument corresponding to a sensor id code, do 
% id_z_sn = all_inst_table(id, -99, -99); 

%initialise table
id = id(:); z = z(:); sn = sn(:);
id_z_sn = table(id,z,sn);

% find the index number of Microcats
%iiMC = find(id == 337 | id == 334);
iiMC = find(id>=332 & id<=337 & id~=335);
% find the index number of ODO Microcats
iiODOMC = find(id == 335);
% find the index number of RBRs
iiRBR = find(id == 330);
% find the index number of Idronauts
iiIDR = find(id == 339);
% find the index number of S4s
iiS4 = find(id == 302);
% and find index number of RCM11s
iiRCM11 = find(id == 310);
% and find index number of Sontek Argonauts
iiARG = find(id == 366);
% and find index number of Nortek Aquadopps
iiNOR = find(id == 368 | id==370);
% and find index number of TRDI DVSs
iiDVS = find(id == 369);
% and find index number of both Seabird and Ixsea BPRs
iiBPR = find(id == 465 | id==470);
% SeaGuard (BPR)
iiSG = find(id == 301);
% SeaPhox
iiSP = find(id == 375);

%instrument abbreviation
id_z_sn.inst = cell(length(id),1); id_z_sn.instl = id_z_sn.inst;
id_z_sn.inst(iiMC) = {'MC'}; id_z_sn.instl(iiMC) = 'microcat';
id_z_sn.inst(iiODOMC) = {'ODOMC'}; id_z_sn.instl(iiMC) = 'microcat';
id_z_sn.inst(iiRBR) = {'RBR'}; id_z_sn.instl(iiMC) = 'rbr';
id_z_sn.inst(iiIDR) = {'IDR'}; id_z_sn.instl(iiMC) = 'idr';
id_z_sn.inst(iiS4) = {'S4'}; id_z_sn.instl(iiMC) = 's4';
id_z_sn.inst(iiRCM11) = {'RCM11'}; id_z_sn.instl(iiMC) = 'rcm11';
id_z_sn.inst(iiARG) = {'ARG'}; id_z_sn.instl(iiMC) = 'argonaut';
id_z_sn.inst(iiNOR) = {'NOR'}; id_z_sn.instl(iiMC) = 'nor';
id_z_sn.inst(iiDVS) = {'DVS'}; id_z_sn.instl(iiMC) = 'dvs';
id_z_sn.inst(iiBPR) = {'BPR'}; id_z_sn.instl(iiMC) = 'bpr';
id_z_sn.inst(iiSG) = {'SG'}; id_z_sn.instl(iiMC) = 'seaguard';
id_z_sn.insts(iiSP) = {'SP'}; id_z_sn.instl(iiMC) = 'seaphox';

%variables (not all, just those used by stats_table and *_overlay.m)
id_z_sn.vars = repmat({'t:c:p'},length(id),1);
id_z_sn.vars(iiODOMC) = {'t:c:p:o2'};
id_z_sn.vars([iiS4;iiRCM11;iiSG]) = {'t:c:p:u:v'};
id_z_sn.vars(iiARG) = {'t:p:u:v'};
id_z_sn.vars(iiNOR) = {'t:p:u:v:w'};
id_z_sn.vars(iiBPR) = {'t:p'};
id_z_sn.vars(iiDVS) = {'t:u:v'};
id_z_sn.vars(iiSP) = {'t:p:o2:s:ph:phv'};

%subdirectory (is this used?)***
id_z_sn.dirs = lower(id_z_sn.inst);
id_z_sn.dirs(iiMC) = {'microcat'};
id_z_sn.dirs(iiODOMC) = {'microcat'};
id_z_sn.dirs(iiRCM11) = {'rcm'};

%filename suffix
id_z_sn.suf = repmat({''},length(id),1);
id_z_sn.suf(iiDVS) = {'_bin2'};
