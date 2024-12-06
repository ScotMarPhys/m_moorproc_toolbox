cd ..
mex -v -f win64/win64_R2007b.bat -output mexnc mexgateway.c netcdf2.c netcdf3.c common.c
cd win64
