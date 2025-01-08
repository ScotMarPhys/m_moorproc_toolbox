% PROCESS_ADCPS is a script to process the adcp data.
close all
global MOORPROC_G
clearvars -except MEXEC_G MOORPROC_G
cruise   = MOORPROC_G.cruise;
operator = MOORPROC_G.operator;
moor = input('mooring deployment (e.g. ebh2_15_2022) to process:   ','s');
plot_interval=[]; %automatic based on available times

pd = moor_inoutpaths('s55',moor);

display('Starting stage 1 adcp2rodb_02')

adcp2rodb_02(moor,pd)


display('Starting stage 2 adcp_raw2use')

%adcp_raw2use_01(moor,pd,'plot_interval',plot_interval)
