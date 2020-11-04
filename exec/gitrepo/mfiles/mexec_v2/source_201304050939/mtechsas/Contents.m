% mtechsas : Library of routines to enable access to techsas files
%
% requires the mexec setup to be configured
% This list created by BAK on JC032 9 April 2009
% 
% Programs mostly begin "mt", and mostly have their own help comments
%
% Contents.m                  this file
%
% These files generally do useful things;
%
% mtdfinfo.m                  start and end times and number of data cycles in a stream
% mtnames.m                   list all the techsas stream names and their mexec short equivalents
% mtvars.m                    list the var names and units in a stream
% mtlast.m                    list the last data cycle in a stream
% mtlookd.m                   list the erliest and latest times in a stream
% mtgetstreams.m              identify all the different streams in a techsas directory             
% mtload.m                    load techsas data into matlab          
% mtposinfo.m                 obtain position from a nav file             
% techsas_to_mstar2.m         create an mstar file from a collection of techsas stream files 
% mtgaps.m                    searches for gaps in a stream
% mtlistit.m                  list data from techsas files
%
% These files are mainly called by the programs above
%
% mtchoosefiles.m             identify full file names required to match a given time            
% mtgetdcrange.m              get range of data cycles within a file that match a given time            
% mtgetdfinfo.m               do the work for mtdfinfo  
% mtgetfiletimes.m            get the time of first and last data cycle in a techsas file      
% mtgetstreamfilenames.m      identify all the file names that match the given stream            
% mtgetvars.m                 do the work for mtvars 
% mtresolve_stream.m          if the input argument is an mexec short stream name, return the full techsas name        
%
% Files called directly by users more often
%
%   mtdfinfo             - function mtdfinfo(instream,tonly)
%   mtgaps               - function mtgaps(instream,g,dn1,dn2,mode)
%   mtgetstreams         - function tstreams = mtgetstreams
%   mtlast               - function [data units] = mtlast(instream)
%   mtlistit             - function mtlistit(instream,ints,dn1,dn2,varlist)
%   mtload               - function [vdata vunits] = mtload(instream,dn1,dn2,mode)
%   mtlookd              - function mtlookd(arg)
%   mtnames              - function matlist = mtnames
%   mtposinfo            - function [lat lon] = mtposinfo(varargin)
%   mtvars               - function mtvars(instream)
%   techsas_to_mstar2    - function techsas_to_mstar(tstream,dn1,dn2)
%
% Files called directly by users less often
%
%   mtchoosefiles        - function fnames = mtchoosefiles(instream,dn1,dn2)
%   mtgetdcrange         - function [dc1 dc2 timem] = mtgetdcrange(techsasfn,dn1,dn2)
%   mtgetdfinfo          - function [firstm lastm numdc] = mgetdfinfo(instream,tonly)
%   mtgetfiletimes       - function [firstm lastm numdc] = mgetfiletimes(fname)
%   mtgetstreamfilenames - function matnames = mtgetstreamfilenames(instream)
%   mtgetvars            - function [vars units] = mtgetvars(instream)
%   mtresolve_stream     - function tstream = mtresolve_stream(instream)

