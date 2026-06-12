function lookup = inst_id_lookup
%
% outputs table connecting instrument abbreviations (inst), instrument long
% names (instl), and instrument IDs (id) as used in info.dat files, as well
% as list of variables used in stats_table and *_overlay.m (vars_st) and
% list of variables used in currents_stacked (vars_cs)


%%%% add inst, id, instl entries as elements in structure

n = 1;
%Microcats
lookup(n).inst = 'MC';
lookup(n).id = [332 333 334 336 337];
lookup(n).instl = 'microcat';

n=n+1;
% ODO Microcats
lookup(n).inst = 'ODOMC';
lookup(n).id = 335;
lookup(n).instl = 'microcat';

n=n+1;
% RBRs
lookup(n).inst = 'RBR';
lookup(n).id = 330;
lookup(n).instl = 'rbr';

n=n+1;
% Idronauts
lookup(n).inst = 'IDR';
lookup(n).id = 339;
lookup(n).instl = 'idr';

n=n+1;
% S4s
lookup(n).inst = 'S4';
lookup(n).id = 302;
lookup(n).instl = 's4';

n=n+1;
% RCM11s
lookup(n).inst = 'RCM11';
lookup(n).id = 310;
lookup(n).instl = 'rcm11';

n=n+1;
% Sontek Argonauts
lookup(n).inst = 'ARG';
lookup(n).id = 366;
lookup(n).instl = 'arg';

n=n+1;
% Nortek Aquadopps
lookup(n).inst = 'NOR';
lookup(n).id = [368 370];
lookup(n).instl = 'nor';

n=n+1;
% TRDI DVSs
lookup(n).inst = 'DVS';
lookup(n).id = 369;
lookup(n).instl = 'dvs';

n=n+1;
% both Seabird and Ixsea BPRs
lookup(n).inst = 'BPR';
lookup(n).id = [465 470];
lookup(n).instl = 'bpr';

n=n+1;
% SeaGuard (BPR)
lookup(n).inst = 'SG';
lookup(n).id = 301;
lookup(n).instl = 'seaguard';

% SeaPhox
n=n+1;
lookup(n).inst = 'SP';
lookup(n).id = 375;
lookup(n).instl = 'seaphox';


%%%% convert to table
lookup = struct2table(lookup);


%%%% add variable lists

lookup.vars_st = repmat({'t:c:p'},length(lookup.inst),1);
lookup.vars_st(strcmp(lookup.inst,'ODOMC')) = {'t:c:p:o2'};
lookup.vars_st(ismember(lookup.inst,{'S4','RCM11','SG'})) = {'t:c:p:u:v'};
lookup.vars_st(strcmp(lookup.inst,'ARG')) = {'t:p:u:v'};
lookup.vars_st(strcmp(lookup.inst,'NOR')) = {'t:p:u:v:w'};
lookup.vars_st(strcmp(lookup.inst,'BPR')) = {'t:p'};
lookup.vars_st(strcmp(lookup.inst,'DVS')) = {'t:u:v'};
lookup.vars_st(strcmp(lookup.inst,'SP')) = {'t:p:o2:s:ph:phv'};

lookup.vars_cs = repmat({' '},length(lookup.inst),1);
lookup.vars_cs(strcmp(lookup.inst,'S4')) = {'yy:mm:dd:hh:u:v:t:c:p:hdg'};
lookup.vars_cs(strcmp(lookup.inst,'RCM11')) = {'yy:mm:dd:hh:ref:u:v:t:c:p:tlt:mss'};
lookup.vars_cs(strcmp(lookup.inst,'ARG')) = {'yy:mm:dd:hh:t:tcat:p:pcat:c:u:v:w:hdg:pit:rol:usd:vsd:wsd:uss:vss:wss:hdgsd:pitsd:rolsd:ipow'};
lookup.vars_cs(strcmp(lookup.inst,'NOR')) = {'yy:mm:dd:hh:t:p:u:v:w:hdg:pit:rol:uss:vss:wss:ipow:cs:cd'};
