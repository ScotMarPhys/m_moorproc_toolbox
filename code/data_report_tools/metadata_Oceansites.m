% for one or more moorings in moor (char or cell array), 
% write mooring metadata to OceanSITES-format .csv file
% metadata come from info.dat files and include position and deployment and
%   recovery dates and ships
%
% 19/04/23 - YLF
%

datadir = fullfile(pathosnap,'data','moor','proc');
outdir = fullfile(pathosnap,'data','moor','export_meta');
outfile = fullfile(outdir, ['osnap_' mcruise '_oceanSITES.csv']);

wm = 'overwrite'; 
%wm = 'append';

if ~iscell(moor)
    moor = {moor};
end

clear t0

% instrument make/model and numbers in rodb files
insts = {'SEABIRD_SBE37' 337;
    'RBR' 330;
    'IDRONAUT' 339;
    'S4' 302;
    'RCM11' 310;
    'SONTEK_ARGONAUT' 366;
    'NORTEK_AQUADOP' [368 370];
    'AADI_SEAGUARD' 301;
    'ADCP' [319:328]};
insts_all = []; instn_all = {};
insts_all = [337 330 339 302 310 366 368 370 301 319:328];
instn_all = cell(size(insts_all));
for no = 1:size(insts,1)
    instn_all = [instn_all; repmat(inst{no,1},length(inst{no,2}),1)];
    insts_all = [insts_all; inst{no,2}'];
end

for mno = 1:length(moor)
    
    clear t0
    t0.REF = 'OSNAP';
    % Load vectors of mooring information
    % id instrument id, sn serial number, z nominal depth of each instrument
    % s_t, e_t, s_d, e_d start and end times and dates
    % lat lon mooring position, wd corrected water depth (m)
    % mr mooring name
    [id,sn,z,s_t,s_d,e_t,e_d,t0.lat,t0.lon,wd,mr]  =  rodbload(infofile,...
        'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');
    
    t0.dn_start = datenum([s_d' s_t' 0]);
    t0.dn_end = datenum([e_d' e_t' 0]);
    
    t = struct2table(t0);
    t = repmat(struct2table(t0,length(sn),1));
    
    %sort by depth
    a = [id z];
    a = sortrows(a,2);
    
    vecMC = sn(iiMC);
    zMC = z(iiMC);
    
    depths(:,1) = id([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG;iiADCP]);
    depths(:,2) = z([iiMC;iiRBR;iiIDR;iiS4;iiRCM11;iiARG;iiNOR;iiSG;iiADCP]);
    depths=sortrows(depths,2);
    iiiMC=find(depths(:,1)==337);
    iiiRBR=find(depths(:,1)==330);
    iiiIDR=find(depths(:,1)==339);
    iiiS4=find(depths(:,1)==302);
    iiiRCM11=find(depths(:,1)==310);
    iiiARG=find(depths(:,1)==366 | depths(:,1)==366337);
    iiiNOR=find(depths(:,1)==368|depths(:,1)==370);
    iiiSG=find(depths(:,1)==301);
    iiiADCP=find(depths(:,1)>=319 & depths(:,1)<=328);
    iii=[iiiS4;iiiRCM11;iiiARG;iiiNOR;iiiMC;iiiRBR;iiiIDR;iiiSG;iiiADCP];
    
    %add code to
    
    writetable(t, outfile, 'WriteMode', wm)
end
