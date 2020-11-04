function k = m_rowcol_to_index(row,col,h,varnum)

% convert k to row and col

nrows = h.dimrows(varnum);
ncols = h.dimcols(varnum);

k = row + nrows * (col-1);
return