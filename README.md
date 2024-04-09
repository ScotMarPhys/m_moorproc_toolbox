Matlab tools for working with instrument data from (OSNAP, RAPID) deep water moorings. The functions contained within this repository take data from raw to gridded and merged time series. Utilities for mooring positioning are also included. Based on the RAPID-MOC and RODB toolboxes, developed at IfM Kiel, NOC, and SAMS. 

To start, run moor_setup.m to point to the locations of the data (and output files) on your system. Run before processing or (suggested) call it in your MATLAB startup.m file, for example:
moor_setup('datadir','/data/pstar/projects/osnap','project','OSNAP')
or if run without input arguments it will prompt the user 

For documentation,  see Documents/Processing_documents/mooring_processing.pdf

For more info on the OSNAP data, see [website](https://scotmarphys.github.io/ScotMarPhys.OSNAP-Mooring-Processing.io/)

### Contact

Contribute by opening an issue.
