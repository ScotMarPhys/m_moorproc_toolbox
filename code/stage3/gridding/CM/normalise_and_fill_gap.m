function y_pred = normalise_and_fill_gap(x,y)
%          """
%          fill gap in y by normalising x scaled with y mean and y std
%          x    : 1st signal
%          y    : 2nd signal - with gabs (nan) to fill with x
% 
%          returns
%          y_pred : fill y
%          """
        % First normalise the variable
        
        
        ibad   = find(isnan(y));
        xnorm  = (x(ibad) - nanmean(x)) / (nanstd(x));
        y_pred = y;
        y_pred(ibad) = (xnorm.*nanstd(y))+nanmean(y);

end