function new_ylim = dynamic_ylim(current_ylim, data)
    % Find the maximum absolute value in the data (handling matrices or vectors)
    max_val = max(abs(data(:)), [], 'omitnan');
    
    % Find the current boundary (assuming it is centered around zero)
    current_bound = max(abs(current_ylim));
    
    % If data exceeds current limits, expand by 10%
    if max_val > current_bound
        new_bound = max_val * 1.1;
        new_ylim = [-new_bound, new_bound];
    else
        % Keep the original limits if data fits inside
        new_ylim = current_ylim;
    end
    
    % Apply to the current axes automatically
    ylim(new_ylim);
end
