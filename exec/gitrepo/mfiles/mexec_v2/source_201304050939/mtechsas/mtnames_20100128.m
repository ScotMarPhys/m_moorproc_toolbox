function matlist = mtnames
% function matlist = mtnames
%
% approximate triplets of mexec short names, rvs streams and techsas streams 
% on JC032
% If called with no output argumnets, list is printed to terminal.
%
% entries are
% mexec short name; rvs name; techsas name

% JC032. If you need to add lines, that is harmless. If you need a whole
% new set of correspondences, retain this list but comment it out, and add
% your new list.

m_common

matlist = {};

if ~isempty(strfind(MEXEC_G.PLATFORM_IDENTIFIER,'Cook'))
 matlist = {
             'adupos'                  ' '                  'ADUPOS-ADUPOS_JC1.gps'
            'smartsv'            'smartsv'                          'AML-AMLSV.SVP'
            'gravity'            'gravity'                  'AirSeaII-S84_JC1.grav'
             'ea600m'             'ea600m'                  'EA600-EA600_JC1.EA600'
            'gyropmv'            'gyropmv'                 'GyroJC-GYRO1_JC1.gyrJC'
             'gyro_s'             'gyro_s'                 'GyroJC-SGYRO_JC1.gyrJC'
              'winch'              'winch'                 'JCWinch-CLAM_JC1.winch'
          'surflight'            'surfmet'              'Light-JC-SM_JC1.SURFMETv2'
            'surfmet'            'surfmet'                'MET-JC-SM_JC1.SURFMETv2'
              'SBE45'              'SBE45'                          'SBE-SBE45.TSG'
            'surftsg'            'surfmet'               'Surf-JC-SM_JC1.SURFMETv2'
            'adu5pat'            'adu5pat'                    'gppat-GPPAT_JC1.att'
           'posmvpos'           'posmvpos'          'position-Applanix_GPS_JC1.gps'
             'dps116'             'dps116'               'position-DPS-116_JC1.gps'
             'seapos'                  ' '            'position-Seapath200_JC1.gps'
             'usbpos'                  ' '                  'position-usbl_JC1.gps'
       'satinfoposmv'                  ' '     'satelliteinfo-Applanix_GPS_JC1.gps'
         'satinfodps'                  ' '          'satelliteinfo-DPS-116_JC1.gps'
         'satinfosea'                  ' '       'satelliteinfo-Seapath200_JC1.gps'
         'satinfousb'                  ' '             'satelliteinfo-usbl_JC1.gps'
                'mag'                  ' '              'scalar_mag-SeaSpy_JC1.mag'
           'posmvtss'           'posmvtss'       'shipattitude-Aplanix_TSS_JC1.att'
             'attsea'                  ' '      'shipattitude-Seapath200AT_JC1.att'
           'attposmv'                  ' '   'shipattitude_aux-Aplanix_TSS_JC1.att'
          'attseaaux'                  ' '  'shipattitude_aux-Seapath200AT_JC1.att'
           'log_skip'           'log_skip'                 'vdvhw-log_skip_JC1.log'
            'log_chf'            'log_chf'                  'vmvbw-log_chf_JC1.log'
            };
end

if ~isempty(strfind(MEXEC_G.PLATFORM_IDENTIFIER,'Discovery'))
matlist = {           
    'gps_g12'       'gps_g12'  'ADUPOS-G12PAT.gps'
    'adupos'              ' '  'ADUPOS-PAPOS.gps'
    'adu5pat'       'gps_ash'  'gppat-GPPAT.att'
    'winch'           'winch'  'DWINCH-CLAM.DWINCH'
    'ea600m'          'ea500'  'PES-Simrad.PES'
    'surflight'     'surfmet'  'Light-SURFMET.SURFMETv2'
    'surfmet'       'surfmet'  'MET-SURFMET.SURFMETv2'
    'surftsg'       'surfmet'  'Surf-SURFMET.SURFMETv2'
    'SBE45'         'seabird'  'SBE45-SBE45.TSG'
    'gyro_s'           'gyro'  'gyro-GYRO.gyr'
    'log_chf'       'log_chf'  'DYLog-LOGCHF.DYLog'
    'gps4000'      'gps_4000'  'position-4000.gps'
    'satinfo4000'         ' '  'satelliteinfo-4000.gps'
    };
end

% Original names
%if ~isempty(strfind(MEXEC_G.PLATFORM_IDENTIFIER,'Discovery'))
%matlist = {           
%     'adu5pat'     ' '   'ADU2-ASH.gps'
%     'winch'       ' '   'DWINCH-CLAM.DWINCH'
%     'ea600m'      ' '   'PES-Simrad.PES'
%     'surfmet'     ' '   'SURFMET-Surfmet.met'
%     'usbl'        ' '   'USBL-USBL01.usbl'
%     'gyro_s'      ' '   'gyro-GYRO.gyr'
%     'log_chf'     ' '   'logchf-log.logchf'
%     'gps4000'     ' '   'position-4000.gps'
%     'satinfo4000' ' '   'satelliteinfo-4000.gps'
%    };
%end

if nargout > 0; return; end

fprintf(1,'\n%20s %20s %45s\n\n',['mexec short name'],['rvs stream name'],['techsas stream name']);

for kstream = 1:size(matlist,1)
fprintf(1,'%20s %20s %45s\n',['''' matlist{kstream,1} ''''],['''' matlist{kstream,2} ''''],['''' matlist{kstream,3} '''']);
end
