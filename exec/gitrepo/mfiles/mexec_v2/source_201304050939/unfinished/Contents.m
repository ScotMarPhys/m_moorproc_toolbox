% UNFINISHED
%
% Files
%   ._m_copy_variable       - 
%   ._m_verson              - 
%   m_adjtim_dataorg        - adjust a time variable from mstar_time_origin to data_time_origin
%   m_adjtim_mstarorg       - adjust a time variable from data_time_origin to mstar_time_origin
%   m_adjtime               - adjust a time variable from data_time_origin in h1 to data_time_origin in
%   m_assign_to_bins        - Find index boundaries of elements of z that lie in each bin
%   m_autolims              - return auto scaled lims for plotting 
%   m_cd                    - not function
%   m_check_nc_varname      - convert any characters found in var names in nc files that are illegal or troublesome in matlab
%   m_clear_global_args     - cleanest to clear global args, especially MEXEC_A.MARGS_IN before creating a
%   m_copy_variable         - Copy a variable from input file to output file, with optional new name in
%   m_copy_variable_test    - Copy a variable from input file to output file, with optional new name in
%   m_copy_variableq        - Copy a variable from input file to output file, with optional new name in
%   m_degmin_from_decdeg    - return degrees and decimal minutes from decimal degrees
%   m_ededit                - called from mplxyed
%   m_edfinddc              - 
%   m_editheader            - edit mstar header
%   m_edlist                - 
%   m_edplot                - This is all the work from mplotxy, but without the
%   m_edplot_old            - This is all the work from mplotxy, but without the
%   m_edrefresh             - 
%   m_edzoom                - 
%   m_exitifopen            - Check if an mstar file is already open for write and exit with error if
%   m_figure                - builtin('figure')
%   m_findvarnum            - seek exact match for vname in h.fldnam
%   m_finis                 - allow user to specify version increment, can be zero when loading pstar
%   m_flag_monotonic        - set a flag to be zero where a variable is not strictly monotonic
%   m_getfilename           - get the name of an mstar file. If there is an argument, unpack it; if
%   m_getinput              - Call with type = 's' for string or 'd' for double
%   m_getinput2             - Call with type = 's' for string or 'd' for double
%   m_getpstarfilename      - get the name of a pstar file. If there is an argument, unpack it; if
%   m_getsites              - 
%   m_gettechsasfilename    - get the name of a techsas file. If there is an argument, unpack it; if
%   m_getvlist              - getvlist: compare character string 'list' with mstar variable names in h.fldnam, and produce list of variable numbers
%   m_global_args           - 
%   m_global_paths          - 
%   m_index_to_rowcol       - convert k to row and col
%   m_interp                - Interpolate x,y onto xi.
%   m_ismstar               - Check if a file is an mstar NetCDF file
%   m_isunitdays            - determine if variable unit matches anything in the list of recognised
%   m_isunitsecs            - determine if variable unit matches anything in the list of recognised
%   m_isvartime             - determine if variable name matches anything in the list of recognised
%   m_matlab_to_mstar       - create mstar file containing vars that are the arguments of the function
%   m_median_despike        - function dataout = m_median_despike(data,s)
%   m_monit                 - presumably should monitor use of programs
%   m_niceticks             - 
%   m_numdims               - test if array is 1d or 2d
%   m_openin                - Check a file is a suitable input file
%   m_openio                - Check a file is a suitable input file, then
%   m_openot                - Create a new file for mstar output
%   m_proghd                - 
%   m_recolor               - jc032 script to allow changes of colormap
%   m_remove_outside_spaces - remove any outside spaces from a string
%   m_remove_spaces         - remove any spaces from a string
%   m_resolve_filename      - sort out whether a filename is already a structure ncfile.name or just the
%   m_rowcol_to_index       - convert k to row and col
%   m_step_fn               - Produces a step function.  =1 before t_switchover, =0 after. 
%   m_time                  - 
%   m_time_to_ymdhms        - convert time in days or seconds, as identified by h.fldunt
%   m_uprlwrq               - determine upper and lower limit of data and number of fillvalues
%   m_verson                - allow user to specify version file increment. Can specify as zero when
%   m_verson2               - allow user to specify version file increment. Can specify as zero when
%   m_verson_jc032          - allow user to specify version file increment. Can specify as zero when
%   m_wherenan              - makes var2=nan whenever var1=nan.
%   m_write_header          - write header to mstar file (global attributes)
%   m_write_history         - write processing steps in the calling program to a history file
%   maddvars                - add vars from one file to a second file
%   mapend                  - simple version, assumes all files have matching variables
%   mavmed                  - average data into bins
%   mavrge                  - average data into bins
%   maxmerc                 - function p2 = axmerc(p)
%   mcalc                   - perform a calculation g = f(x,y,z,....)
%   mcalcq                  - perform a calculation g = f(x,y,z,....)
%   mcalib                  - calibrate a variable, writing output back to the same variable
%   mcalib2                 - calibrate a variable, writing output back to the same variable
%   mcd                     - function newdir = mcd(varargin)
%   mcd_old                 - script to cd to the directory stored in global variable MEXEC_G.MEXEC_CWD
%   mchangetimeorigin       - change data time origin; adjust time data as required.
%   mcontr                  - function mcontr(cdfin,ncfile)
%   mcontr_old              - function mcontr(cdf,ncfile)
%   mcontr_old2             - function mcontr(cdf,ncfile)
%   mcontr_old3             - function mcontr(cdfin,ncfile)
%   mcontr_old4             - function mcontr(cdfin,ncfile)
%   mcontrnew               - function mcontr(cdfin,ncfile)
%   mcopya                  - copy selected vars and data cycles to a new output file
%   mcsetd                  - function mcsetd(varargin)
%   mdatpik                 - pick data cycles depending on whether control variable is inside or
%   medita                  - Edit data that lie outside a range to absent value
%   mgridp                  - simple version, assumes all files have matching variables
%   mgridp_flags            - simple version, assumes all files have matching variables
%   mgridp_old              - simple version, assumes all files have matching variables
%   mheadr                  - change header data
%   mhisto                  - plot histogram
%   mhistory                - type the history entry for the most program
%   mintrp                  - interpolate a set of variables to fill gaps
%   mlist                   - list data
%   mlistf                  - list data
%   mlisth                  - 
%   mload                   - load data and header contents of mstar NetCDF file into structure
%   mmerge                  - merge vars from a second file onto the first, using a control variable
%   mmerge_bim              - merge vars from a second file onto the first, using a control variable
%   mpaste                  - paste vars from a second file onto the first, optionally using a control variable
%   mplotxy                 - function mplotxy(pdf,ncfile)
%   mplxyed                 - function mplxyed
%   mplxyed_old             - function mplxyed
%   mpsetd                  - function mpsetd(varargin)
%   mreset                  - reset the write flag on an mstar file
%   msave                   - save matlab vars as mstar file
%   msbe_to_mstar           - load sbe ctd  file into mstar file
%   mskeleton_in1_ot1       - skeleton program with input and output to different files
%   mskeleton_inot          - skeleton program with input and output to same file
%   mskeleton_mcalc         - script to convert u,v to speed,direction
%   msort                   - sort data 
%   muvsd                   - script to convert u,v to speed,direction
%   pstar_to_mstar          - load pstar file into mstar file
%   techsas_to_mstar        - load techsas file into mstar file
