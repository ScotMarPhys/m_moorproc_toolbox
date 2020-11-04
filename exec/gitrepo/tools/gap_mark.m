%function gapI = gap_mark(vec,gap_max,iss)

function [gapI] = gap_mark(vec,gap_max,iss)

   [a,b]    = consec_nan(vec); 
   gap      = b/iss;
   gapI     = find(gap>gap_max);
   if ~isempty(gapI)
     gapI        = igrep(sort([a(gapI) a(gapI)+b(gapI)-1]));
   else
     gapI        = [];   
   end