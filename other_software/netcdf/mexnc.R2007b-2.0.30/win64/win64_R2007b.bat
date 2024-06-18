@echo off
rem MSVC80OPTS.BAT
rem
rem    Compile and link options used for building MEX-files
rem    using the Microsoft Visual C++ compiler version 8.0
rem
rem    $Revision: 1.1.6.10 $  $Date: 2006/11/15 14:50:11 $
rem
rem ********************************************************************
rem General parameters
rem ********************************************************************
rem
set UNIDATA_INC=D:\sandbox\johne\Amake.3p-tmw\3p\install\unknown\win64\netcdf\include
set UNIDATA_LIB=D:\sandbox\johne\Amake.3p-tmw\3p\install\unknown\win64\netcdf\lib
set UNIDATA_LIBS=D:\sandbox\johne\Amake.3p-tmw\3p\install\unknown\win64\netcdf\lib\netcdf.lib

set MATLAB=%MATLAB%
set VS80COMNTOOLS=%VS80COMNTOOLS%
set VSINSTALLDIR=C:\Program Files (x86)\Microsoft Visual Studio 8
set VCINSTALLDIR=%VSINSTALLDIR%\VC
set PATH=%VCINSTALLDIR%\BIN\amd64;%VCINSTALLDIR%\PlatformSDK\bin;%VCINSTALLDIR%\PlatformSDK\bin\win64\amd64;%VSINSTALLDIR%\SDK\v2.0\bin;%VSINSTALLDIR%\Common7\Tools;%VSINSTALLDIR%\Common7\Tools\bin;%MATLAB_BIN%;%PATH%
set INCLUDE=%UNIDATA_INC%;%VCINSTALLDIR%\ATLMFC\INCLUDE;%VCINSTALLDIR%\INCLUDE;%VCINSTALLDIR%\PlatformSDK\INCLUDE;%VSINSTALLDIR%\SDK\v2.0\include;%INCLUDE%
set LIB=%VCINSTALLDIR%\LIB\amd64;%VCINSTALLDIR%\ATLMFC\LIB\amd64;%VCINSTALLDIR%\PlatformSDK\lib\amd64;%VSINSTALLDIR%\SDK\v2.0\lib\amd64;%MATLAB%\extern\lib\win64;%LIB%
set MW_TARGET_ARCH=win64

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=cl
rem set COMPFLAGS=-c -I"%UNIDATA_INC%" -Zp8 -GR -W3 -EHs -D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE -D_SECURE_SCL=0 -DMATLAB_MEX_FILE -nologo -DDLL_NETCDF
set COMPFLAGS=-c -ID:\sandbox\johne\Amake.3p-tmw\3p\install\unknown\win64\netcdf\include -Zp8 -GR -W3 -EHs -D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE -D_SECURE_SCL=0 -DMATLAB_MEX_FILE -nologo -DDLL_NETCDF
set OPTIMFLAGS=/MD -O2 -Oy- -DNDEBUG
set DEBUGFLAGS=/MD -Zi -Fd"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set NAME_OBJECT=/Fo

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LIBLOC=%MATLAB%\extern\lib\win64\microsoft
set LINKER=link
set LINKFLAGS=/dll /export:%ENTRYPOINT% /MAP /LIBPATH:"%UNIDATA_LIB%" %UNIDATA_LIBS% /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /implib:%LIB_NAME%.x /MACHINE:AMD64 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
set LINKOPTIMFLAGS=
set LINKDEBUGFLAGS=/DEBUG /PDB:"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=/out:"%OUTDIR%%MEX_NAME%%MEX_EXT%"
set RSP_FILE_INDICATOR=@

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=rc /fo "%OUTDIR%mexversion.res"
set RC_LINKER=

set POSTLINK_CMDS=del %LIB_NAME%.x
set POSTLINK_CMDS1=mt "-outputresource:%OUTDIR%%MEX_NAME%%MEX_EXT%;2" -manifest "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest"
set POSTLINK_CMDS2=del "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest"
set POSTLINK_CMDS3=del "%OUTDIR%%MEX_NAME%.map"
