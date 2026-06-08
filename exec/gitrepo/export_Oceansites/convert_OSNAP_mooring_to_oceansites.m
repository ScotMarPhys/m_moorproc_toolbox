function convert_OSNAP_mooring_to_oceansites(depyear,institution,varargin)
% convert_OSNAP_mooring_to_oceansites(depyear, institution,varargin)
% convert OSNAP_mooring_to_oceansites('07_2022', 'SAMS')
% convert_OSNAP_mooring_to_oceansites('07_2022', 'SAMS', {'rtwb2'})
%
% for each of three instrument types (microcat, nortek, adcp), were
% relevant for the moorings to be processed, set metadata specific to
% project, then specific to institution, then specific to
% mooring/sensor/processor 
%
% then call rodb_to_oceansitesnetcdf

%which moorings to work on
%defaults
switch institution
    case 'SAMS'
        moorlist = {'rteb1','rtwb1','rtwb2'};
    case 'NOC'
        moorlist = {'ib3','ib4','ib5'};
end
%optional input argument
if nargin>2
    moorlist = varargin{1};
end

% metadata that doesn't change (much)
clear m
m.project = 'OSNAP';
m.array  ='OSNAP';
m.network = 'OSNAP';
m.history = 'Delayed time processed quality controlled';
switch institution
    case 'SAMS'
        moorlist = {'rteb1','rtwb1','rtwb2'};
        m.principal_investigator = ''; 
        m.principal_investigator_email = ''; 
        m.principal_investigator_url = 'http://www.sams.ac.uk'; 
        m.institution = 'Scottish Association for Marine Science, SAMS';
    case 'NOC'
        moorlist = {'ib3','ib4','ib5'};
        m.principal_investigator = 'Penny Holliday'; 
        m.principal_investigator_email = 'penny.holliday@noc.ac.uk'; 
        m.principal_investigator_url = 'https://www.noc.ac.uk';
        m.institution = 'National Oceanography Centre, NOC';
end

%now copy-past these fields to moorinfo.mcat etc.
insts = {'mcat','nortek','adcp'};
for no = 1:length(insts)
    moorinfo.(insts{no}) = m;
end

%% MCAT
moorinfo.mcat.processing_level = 'calibration using pre- and post- deployment calibration casts; data manually reviewed';
moorinfo.mcat.QC_indicator = 'excellent';
moorinfo.mcat.proclvl             = 'calibrated;good data';
moorinfo.mcat.sensor_model        = 'SBE37-SMP Microcat';
moorinfo.mcat.sensor_manufacturer = 'Sea-Bird Electronics Inc';
moorinfo.mcat.sensor_ref          = 'https://www.seabird.com/asset-get.download.jsa?id=66339519258';
moorinfo.mcat.sensor_mount        = 'mounted_on_mooring_line';
moorinfo.mcat.sens_orientation    = 'vertical';
moorinfo.mcat.presuncertainty     = ' ';
moorinfo.mcat.presaccuracy        =  '+/- 0.1% of full scale range';
moorinfo.mcat.presresolution      =  '0.002 % of full scale range';
moorinfo.mcat.tempuncertainty     = ' ';
moorinfo.mcat.tempaccuracy        =  '+/- 0.002';
moorinfo.mcat.tempresolution      =  '0.0001';
moorinfo.mcat.tempprecision       =  ' ';
moorinfo.mcat.conduncertainty     = ' ';
moorinfo.mcat.condaccuracy        =  '+/- 0.003';
moorinfo.mcat.condresolution      =  '0.0001';
moorinfo.mcat.condprecision       =  ' ';
moorinfo.mcat.saluncertainty      = ' ';
moorinfo.mcat.salaccuracy         =  '+/- 0.002';
moorinfo.mcat.salresolution       =  '0.0001';
moorinfo.mcat.salprecision        =  ' ';
%% NORTEK
moorinfo.nortek.processing_level = 'speed of sound and magnetic deviation corrections; data manually reviewed';
moorinfo.nortek.QC_indicator = 'excellent';
moorinfo.nortek.proclvl             = 'calibrated;good data';
moorinfo.nortek.sensor_model        = 'Aquadopp current meter';
moorinfo.nortek.sensor_manufacturer = 'Nortek AS';
moorinfo.nortek.sensor_ref          = 'https://www.nortekgroup.com/products/aquadopp2-500m/pdf';
moorinfo.nortek.presaccuracy         = '0.5%';
moorinfo.nortek.presresolution       = '0.005% of full scale';
moorinfo.nortek.sensor_mount        = 'mounted_on_mooring_line';
moorinfo.nortek.sens_orientation    = 'downward';
moorinfo.nortek.velaccuracy         = '1% of measured value +/-0.5 cm/s';
moorinfo.nortek.velresolution       = 'N/A';
%% ADCP
moorinfo.adcp.processing_level = 'speed of sound and magnetic deviation corrections; data manually reviewed';
moorinfo.adcp.QC_indicator = 'excellent';

%metadata based on institution/PI, plus which moorings to process


for ijk= 1:length(moorlist)
    moor = [moorlist{ijk} '_' depyear];
    switch moorlist{ijk}

        case 'rtadcp1'
            moorinfo.adcp.contributor_name = 'Lewis Drysdale; ';	% CHANGE
            moorinfo.adcp.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.adcp.contributor_email = 'lewis.drysdale@sams.ac.uk';	% CHANGE
            moorinfo.adcp.proclvl             = 'calibrated;good data';
            moorinfo.adcp.sensor_model        = 'Workhorse Long Ranger 75 kHz';
            moorinfo.adcp.sensor_manufacturer = 'Teledyne RDI';
            moorinfo.adcp.sensor_ref          = 'http://rdinstruments.com/__documents/Brochures/long_ranger_datasheet_lr.pdf';
            moorinfo.adcp.velaccuracy         = '+- 1% of the water velocity relative to ADCP +- 5mm/s';
            moorinfo.adcp.velresolution       = '1mm/s';
            moorinfo.adcp.sensor_mount        = 'benthos_frame';
            moorinfo.adcp.sens_orientation    = 'upward';
            moorinfo.adcp.coordsyst           = 'earth';
            moorinfo.adcp.dist1stbin          =  24.35;
            moorinfo.adcp.binsize             =	16;
            moorinfo.adcp.nbbins              =	60;
            moorinfo.adcp.pingperens          =	30;
            moorinfo.adcp.timepinggroup       =	'120s';

        case {'rteb1','rtwb1','rtwb2'}

            % Microcat
            moorinfo.mcat.contributor_name = 'Lewis Drysdale; Sam Jones ';	% CHANGE
            moorinfo.mcat.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.mcat.contributor_email = ' lewis.drysdale@sams.ac.uk; sam.jones@sams.ac.uk';	% CHANGE

            % Nortek
            moorinfo.nortek.contributor_name = ' Lewis Drysdale; Sam Jones;';	% CHANGE
            moorinfo.nortek.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.nortek.contributor_email = 'lewis.drysdale@sams.ac.uk; sam.jones@sams.ac.uk; ';	% CHANGE

        case {'ib3','ib4','ib5'}

            % Microcat
            moorinfo.mcat.contributor_name = 'Tiago Segabinazzi Dotto; Yvonne Firing';	% CHANGE
            moorinfo.mcat.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.mcat.contributor_email = ' tiago.dotto@noc.ac.uk; yvonne.firing@noc.ac.uk';	% CHANGE

            % Nortek
            moorinfo.nortek.contributor_name = ' Tiago Segabinazzi Dotto; Darren Rayner';	% CHANGE
            moorinfo.nortek.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.nortek.contributor_email = 'tiago.dotto@noc.ac.uk; darren.rayner@noc.ac.uk';	% CHANGE

            % ADCP
            moorinfo.adcp.contributor_name = ' Tiago Segabinazzi Dotto; Darren Rayner';	% CHANGE
            moorinfo.adcp.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.adcp.contributor_email = 'tiago.dotto@noc.ac.uk; darren.rayner@noc.ac.uk';	% CHANGE
            moorinfo.adcp.proclvl             = 'calibrated;good data';
            moorinfo.adcp.sensor_model        = 'Workhorse 300 kHz';
            moorinfo.adcp.sensor_manufacturer = 'Teledyne RDI';	
            moorinfo.adcp.sensor_ref          = 'http://rdinstruments.com/__documents/Brochures/sentinel_datasheet_lr.pdf';
            moorinfo.adcp.velaccuracy         = '+- 1% of the water velocity relative to ADCP +- 5mm/s';
            moorinfo.adcp.velresolution       = '1mm/s';
            moorinfo.adcp.sensor_mount        = 'mounted_on_mooring_line';
            moorinfo.adcp.sens_orientation    = 'downward';
            moorinfo.adcp.coordsyst           = 'earth';
            moorinfo.adcp.dist1stbin          =  6.2; 
            moorinfo.adcp.binsize             =	4;
            moorinfo.adcp.nbbins              =	28;
            moorinfo.adcp.pingperens          =	42;
            moorinfo.adcp.timepinggroup       =	'85s';   


    end


    rodb_to_oceansitesnetcdf(moor,moorinfo)

end

