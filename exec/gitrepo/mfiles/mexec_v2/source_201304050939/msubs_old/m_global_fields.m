% in mexec v2 which uses structures for variables required globally:
%
% in mexec v1, a variable was initialised to empty by declaring it global.
% in mexec v2 this doesnt happen, so some programs can crash by expecting a
% field to exist. Therefore ensure all fields exist in m_setup
%

MEXEC_G.SITE 
MEXEC_G.VERSION_FILE 
MEXEC_G.HISTORY_DIRECTORY
MEXEC_G.Mhousekeeping_version
MEXEC_G.PLATFORM_TYPE 
MEXEC_G.PLATFORM_IDENTIFIER 
MEXEC_G.PLATFORM_NUMBER
MEXEC_G.MSTAR_TIME_ORIGIN 
MEXEC_G.COMMENT_DELIMITER_STRING
MEXEC_G.MUSER
MEXEC_G.Mcarter
MEXEC_G.MEXEC_DATA_ROOT 
MEXEC_G.MEXEC_CWD 
MEXEC_G.MDIRLIST
MEXEC_G.PEXEC_DATA_ROOT 
MEXEC_G.PEXEC_CWD 
MEXEC_G.PDIRLIST
MEXEC_G.Mrsh_machine 
MEXEC_G.Mrsh_dataroot_local 
MEXEC_G.Mrsh_dataroot_remote
MEXEC_G.MSCRIPT_CRUISE_STRING 
MEXEC_G.Mtechsas_root 
MEXEC_G.Mtechsas_torg 
MEXEC_G.Mtechsas_default_navstream
MEXEC_G.MDEFAULT_DATA_TIME_ORIGIN
MEXEC_G.Muse_version_lockfile
MEXEC_G.Mshipdatasystem
MEXEC_G.Mscs_root 
MEXEC_G.Mscs_mat 
MEXEC_G.Mscs_sed 
MEXEC_G.Mscs_torg 
MEXEC_G.Mscs_default_navstream
MEXEC_G.Mship


MEXEC_A.MARGS_IN 
MEXEC_A.MARGS_OT 
MEXEC_A.MARGS_IN_OLD 
MEXEC_A.MARGS_IN_LOCAL 
MEXEC_A.MARGS_IN_LOCAL_OLD
MEXEC_A.Mhistory_in 
MEXEC_A.Mhistory_ot 
MEXEC_A.Mprog 
MEXEC_A.Mhistory_filename 
MEXEC_A.Mhistory_lastlines
MEXEC_A.Mhistory_skip
MEXEC_A.Mfidterm 
MEXEC_A.Mfider
MEXEC_A.MSAVE_VLIST
MEXEC_A.Mtimnames 
MEXEC_A.Mtimunits_days 
MEXEC_A.Mtimunits_seconds


MEXEC_P.Mhousekeeping_history