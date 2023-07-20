% routine to make localised bathy files from etopo2 data for use with rapid_widgit_lite
% will probably update this to use swath for WB and possibly MAR if have
% coverage
moorings=   {'ebh4_5','ebh3', 'ebh2', 'ebh1', 'ebhi', 'eb1',  'nog',  'mar3', 'mar2', 'mar1', 'mar0', 'wb6',  'wb4',  'wbh2', 'wb2',  'wb1',  'wbadcp'};
mooring_lat=[27.8501, 27.805, 27.6147,27.1435,24.933, 23.7577,23.7710,23.807, 24.1775,24.2655,25.1047,26.4950,26.4782,26.4823,26.5136,26.5012,26.52492];
mooring_lon=[13.5408, 13.7468,14.2109,15.3777,21.2732,24.1582,41.0988,41.0984,49.763, 49.75  ,52.0102,70.5223,75.7028,76.6276,76.7379,76.8149,76.86803];

mooring_lat_upper
mooring_lat_lower


[lon lat depth]=textread('RapidAll_4396.xyz','%f %f %f');

lat=reshape(lat,2401,301);
lon=reshape(lon,2401,301);
depth=reshape(depth,2401,301);
lat=lat'; lon=lon'; depth=depth';
lat=flipud(lat); lon-flipud(lon); depth=flipud(depth);
depth=depth*-1;

