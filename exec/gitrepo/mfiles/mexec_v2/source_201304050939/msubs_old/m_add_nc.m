function fn = m_add_nc(fn)

% Use: fn = m_add_nc(fn)
% 
% If the string argument is fn, check to see whether fn already has a
% suffic of .nc.
% If not, add .nc

len = length(fn);


if len >= 3; 
    if strmatch('.nc',fn(end-2:end),'exact') %already has .nc suffix
        return
    end
end

fn = [fn '.nc'];
return

