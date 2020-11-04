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

% list of Cook names significantly changed between JC032 and JC044

m_common

matlist = {};

if ~isempty(strfind(MEXEC_G.PLATFORM_IDENTIFIER,'Cook'))
 matlist = {
%             'adupos'                  ' '                  'ADUPOS-ADUPOS_JC1.gps'
             'adupos'                  ' '                  'PASHRPOS-ADUPOS_JC1.PASHR'
            'smartsv'            'smartsv'                          'AML-AMLSV.SVP'
%             'gravity'            'gravity'                  'AirSeaII-S84_JC1.grav'
            'gravity'            'gravity'                  'AirSeaII-S84_JC1.AirSeaII'
             'ea600m'             'ea600m'                  'EA600-EA600_JC1.EA600'
%            'gyropmv'            'gyropmv'                 'GyroJC-GYRO1_JC1.gyrJC'
            'gyropmv'            'gyropmv'                 'gyro-GYRO1_JC1.gyr'
%             'gyro_s'             'gyro_s'                 'GyroJC-SGYRO_JC1.gyrJC'
             'gyro_s'             'gyro_s'                 'gyro-SGYRO_JC1.gyr'
%              'winch'              'winch'                 'JCWinch-CLAM_JC1.winch'
              'winch'              'winch'                 'CLAM-CLAM_JC1.CLAM'
          'surflight'            'surfmet'              'Light-JC-SM_JC1.SURFMETv2'
            'surfmet'            'surfmet'                'MET-JC-SM_JC1.SURFMETv2'
%              'SBE45'              'SBE45'                          'SBE-SBE45.TSG'
              'SBE45'              'SBE45'                          'SBE45-SBE45_JC1.TSG'
            'surftsg'            'surfmet'               'Surf-JC-SM_JC1.SURFMETv2'
%            'adu5pat'            'adu5pat'                    'gppat-GPPAT_JC1.att'
            'adu5pat'            'adu5pat'                    'GPPAT-GPPAT_JC1.GPPAT'
           'posmvpos'           'posmvpos'          'position-Applanix_GPS_JC1.gps'
             'dps116'             'dps116'               'position-DPS-116_JC1.gps'
             'seapos'             'seapos'            'position-Seapath200_JC1.gps'
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
%            'log_skip'           'log_skip'                 'vdvhw-log_skip_JC1.log'
           'log_skip'           'log_skip'                 'VDVHW-log_skip_JC1.Log'
%            'log_chf'            'log_chf'                  'vmvbw-log_chf_JC1.log'
            'log_chf'            'log_chf'                  'EMLog-log_chf_JC1.EMLog'
            'em120'              '        '               'sb_depth-EM120_JC1.depth'
            'cnav'               'cnav'                              'cnav-CNAV.GPS'
            'satinfocnav'        'satinfocnav'              'satelliteinfo-CNAV.gps'
            
            };
end

if ~isempty(strfind(MEXEC_G.PLATFORM_IDENTIFIER,'Discovery'))
    % before di 368 13 jul 2011
% matlist = {           
%     'gps_g12'       'gps_g12'  'ADUPOS-G12PAT.gps'
%     'adupos'              ' '  'ADUPOS-PAPOS.gps'
%     'adu5pat'       'gps_ash'  'gppat-GPPAT.att'
%     'winch'           'winch'  'DWINCH-CLAM.DWINCH'
%     'ea600m'          'ea500'  'PES-Simrad.PES'
%     'surflight'     'surfmet'  'Light-SURFMET.SURFMETv2'
%     'surfmet'       'surfmet'  'MET-SURFMET.SURFMETv2'
%     'surftsg'       'surfmet'  'Surf-SURFMET.SURFMETv2'
%     'SBE45'         'seabird'  'SBE45-SBE45.TSG'
%     'gyro_s'           'gyro'  'gyro-GYRO.gyr'
%     'log_chf'       'log_chf'  'DYLog-LOGCHF.DYLog'
%     'gps4000'      'gps_4000'  'position-4000.gps'
%     'satinfo4000'         ' '  'satelliteinfo-4000.gps'
%     };
% end

%  di 368 jul 2011
matlist = {           
 'winch' ' '   'CLAM-CLAM.CLAM'
 'log_chf ' ' '   'EMLog-LOGCHF.EMLog'
  'adu5pat' ' '  'GPPAT-GPPAT.GPPAT'
  'surflight' ' '  'Light-SURFMET.SURFMETv2'
   'surfmet' ' ' 'MET-SURFMET.SURFMETv2'
 'adupos' ' '   'PASHRPOS-PAPOS.PASHR'
  'ea600m' ' '  'PES-Simrad_PT1.PES'
  'SBE45' ' '  'SBE45-SBE45.TSG'
 'surftsg' ' '   'Surf-SURFMET.SURFMETv2'
 'gyro_s' ' '   'gyro-GYRO.gyr'
 'gps4000' ' '   'position-4000.gps'
 'gpsfugro' ' '   's9200G2s-FUGRO.GPS'
 'satinfo4000' ' '   'satelliteinfo-4000.gps'
 'satinfofugro' ' '   'satelliteinfo-FUGRO.gps'
 'gps1' ' ' 's9200G2s-GPS1.GPS' % di368 gps splitter test
 'gps2' ' ' 's9200G2s-GPS2.GPS' % di368 gps splitter test
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
