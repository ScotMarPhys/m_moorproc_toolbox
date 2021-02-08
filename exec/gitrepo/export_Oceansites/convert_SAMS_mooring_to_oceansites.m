%clear
close all
%addpath(genpath('~/Dropbox/Work/function_MATLAB/Mooring_Processing_toolbox'))

%-------------------------------------
% Definition of the different path
%pathosnap = pwd;

%basedir  = '~/Dropbox/Work/Dataproc/Postdoc_OSNAP/OSNAP_mooring/'; %
basedir = [execdir 'gitrepo/export_Oceansites'];
% directory with the processed mooring data
% procpath = '/media/SAMS/m/Mar_Phys/OSNAP_mooring_data_processing/osnap/data/moor/proc/' % NEED TO BE UPDATE:[basedir 'backup_mdrive/proc/'];	
%procpath = '~/osnap/data/moor/proc/' % NEED TO BE UPDATE:[basedir 'backup_mdrive/proc/'];	
procpath = [pathosnap '/data/moor/proc/'];
% output directory for the netcdf
%outpath = [basedir 'OSNAP_oceanSITES_format_conversion_matlab_scripts_for_SAMSNOC/oceansites_format/'];
outpath = [basedir 'gitrepo/export_Oceansites/oceansites_format'];
%-------------------------------------
% Selection of the deployment year
% depyear ='01_2014';
% depyear ='02_2015';
% depyear = '03_2016';
% depyear = '04_2017';
depyear='05_2018';
%-------------------------------------
%% Selection of the mooring to process
%moorlist ={'nocm1','nocm2','nocm3','nocm4','nocm5'};
%moorlist ={'nocm2'};
moorlist = {'rteb1','rtwb1','rtwb2'};
%moorlist = {'rtadcp1'};%
 		
% eval(['cd ' basedir]);
cd(basedir)

for ijk= 1:length(moorlist)
    moor = [moorlist{ijk} '_' depyear];
    switch moorlist{ijk}
        case 'rtadcp1'
            moorinfo.adcp.principal_investigator = 'Stuart Cunningham';	% CHANGE
            moorinfo.adcp.principal_investigator_email = 'stuart.cunningham@sams.ac.uk';	% CHANGE
            moorinfo.adcp.principal_investigator_url = 'http://www.sams.ac.uk';	% CHANGE
            moorinfo.adcp.institution = 'Scottish Association for Marine Science, SAMS';	% CHANGE
            moorinfo.adcp.contributor_name = 'Stuart Cunningham; Loic Houpert';	% CHANGE
            moorinfo.adcp.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.adcp.contributor_email = 'stuart.cunningham@sams.ac.uk; loic.houpert@noc.ac.uk';	% CHANGE
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
            moorinfo.mcat.principal_investigator = 'Stuart Cunningham';	% CHANGE
            moorinfo.mcat.principal_investigator_email = 'stuart.cunningham@sams.ac.uk';	% CHANGE
            moorinfo.mcat.principal_investigator_url = 'http://www.sams.ac.uk';	% CHANGE
            moorinfo.mcat.institution = 'Scottish Association for Marine Science, SAMS';	% CHANGE
            moorinfo.mcat.contributor_name = 'Stuart Cunningham; Sam Jones; Loic Houpert; Lewis Drysdale';	% CHANGE
            moorinfo.mcat.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.mcat.contributor_email = 'stuart.cunningham@sams.ac.uk; sam.jones@sams.ac.uk; loic.houpert@noc.ac.uk; lewis.drysdale@sams.ac.uk';	% CHANGE      
            
       % Nortek
            moorinfo.nortek.principal_investigator = 'Stuart Cunningham';	% CHANGE
            moorinfo.nortek.principal_investigator_email = 'stuart.cunningham@sams.ac.uk';	% CHANGE
            moorinfo.nortek.principal_investigator_url = 'http://www.sams.ac.uk';	% CHANGE
            moorinfo.nortek.institution = 'Scottish Association for Marine Science, SAMS';	% CHANGE
            moorinfo.nortek.contributor_name = 'Stuart Cunningham; Sam Jones; Loic Houpert; Lewis Drysdale';	% CHANGE
            moorinfo.nortek.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.nortek.contributor_email = 'stuart.cunningham@sams.ac.uk;  sam.jones@sams.ac.uk; loic.houpert@noc.ac.uk; lewis.drysdale@sams.ac.uk';	% CHANGE      
            
        case {'nocm1','nocm2','nocm3','nocm4','nocm5'};
            
         %Microcat   
            moorinfo.mcat.principal_investigator = 'Penny Holliday';	% CHANGE
            moorinfo.mcat.principal_investigator_email = 'penny.holliday@noc.ac.uk';	% CHANGE
            moorinfo.mcat.principal_investigator_url = 'http://www.noc.ac.uk';	% CHANGE
            moorinfo.mcat.institution = 'National Oceanography Centre, NOC';	% CHANGE
            moorinfo.mcat.contributor_name = 'Penny Holliday; Loic Houpert';	% CHANGE
            moorinfo.mcat.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.mcat.contributor_email = 'penny.holliday@noc.ac.uk; loic.houpert@noc.ac.uk';	% CHANGE              
            
         %Nortek   
            moorinfo.nortek.principal_investigator = 'Penny Holliday';	% CHANGE
            moorinfo.nortek.principal_investigator_email = 'penny.holliday@noc.ac.uk';	% CHANGE
            moorinfo.nortek.principal_investigator_url = 'http://www.noc.ac.uk';	% CHANGE
            moorinfo.nortek.institution = 'National Oceanography Centre, NOC';	% CHANGE
            moorinfo.nortek.contributor_name = 'Penny Holliday; Loic Houpert';	% CHANGE
            moorinfo.nortek.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.nortek.contributor_email = 'penny.holliday@noc.ac.uk; loic.houpert@noc.ac.uk';	% CHANGE         
          
        %ADCP
            moorinfo.adcp.principal_investigator = 'Penny Holliday';	% CHANGE
            moorinfo.adcp.principal_investigator_email = 'penny.holliday@noc.ac.uk';	% CHANGE
            moorinfo.adcp.principal_investigator_url = 'http://www.noc.ac.uk';	% CHANGE
            moorinfo.adcp.institution = 'National Oceanography Centre, NOC';	% CHANGE
            moorinfo.adcp.contributor_name = 'Penny Holliday; Loic Houpert';	% CHANGE
            moorinfo.adcp.contributor_role = 'data processing and interpretation';	% CHANGE
            moorinfo.adcp.contributor_email = 'penny.holliday@noc.ac.uk; loic.houpert@noc.ac.uk';	% CHANGE
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
    
    
    % Nortek:       
    moorinfo.nortek.project = 'OSNAP';
    moorinfo.nortek.array = 'OSNAP';  
    moorinfo.nortek.network = 'OSNAP'; 
    moorinfo.nortek.history = 'Delayed time processed quality controlled';
    moorinfo.nortek.processing_level = 'speed of sound and magnetic deviation corrections; data manually reviewed';
    moorinfo.nortek.QC_indicator = 'excellent';
    moorinfo.nortek.proclvl             = 'calibrated;good data';
    moorinfo.nortek.sensor_model        = 'Aquadopp current meter';
    moorinfo.nortek.sensor_manufacturer = 'Nortek AS';	
    moorinfo.nortek.sensor_ref          = 'http://www.nortek-as.com/lib/brochures/Aquadopp%2006%20b.pdf';
    moorinfo.nortek.presaccuracy         = '0.5%';
    moorinfo.nortek.presresolution       = '0.005% of full scale';     
    moorinfo.nortek.sensor_mount        = 'mounted_on_mooring_line';
    moorinfo.nortek.sens_orientation    = 'downward';            
    moorinfo.nortek.velaccuracy         = '1% of measured value +/-0.5 cm/s';
    moorinfo.nortek.velresolution       = 'N/A';  

    % Microcat
    moorinfo.mcat.project = 'OSNAP';
    moorinfo.mcat.array = 'OSNAP';  
    moorinfo.mcat.network = 'OSNAP'; 
    moorinfo.mcat.history = 'Delayed time processed quality controlled';
    moorinfo.mcat.processing_level = 'calibration using pre- and post- deployment calibration casts; data manually reviewed';
    moorinfo.mcat.QC_indicator = 'excellent';
    moorinfo.mcat.proclvl             = 'calibrated;good data';
    moorinfo.mcat.sensor_model        = 'SBE37-SM Microcat';
    moorinfo.mcat.sensor_manufacturer = 'Sea-Bird Electronics Inc';	
    moorinfo.mcat.sensor_ref          = 'http://www.seabird.com/sbe37sm-microcat-ctd';
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
       
    
    %ADCP:
    moorinfo.adcp.project = 'OSNAP';
    moorinfo.adcp.array = 'OSNAP';  
    moorinfo.adcp.network = 'OSNAP'; 
    moorinfo.adcp.history = 'Delayed time processed quality controlled';
    moorinfo.adcp.processing_level = 'speed of sound and magnetic deviation corrections; data manually reviewed';
    moorinfo.adcp.QC_indicator = 'excellent';



rodb_to_oceansitesnetcdf(moor,procpath,moorinfo)
  	   
end

