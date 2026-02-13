function [Data_corr, h_fig] = s55_hard_iron_compass_correction2(Data, flag_value, apl_cor, plot_data)
% S55_HARD_IRON_COMPASS_CORRECTION Corrects ENU velocities for magnetic interference.
%
% USAGE:
%   [Data_new, h] = s55_hard_iron_compass_correction2(Data, 0, true, true);
%
% INPUTS:
%   Data       - Structure with Magnetometer [Nx2], U, V, and Time.
%   flag_value - QC value to keep (usually 0).
%   apl_cor    - Optional, default is false
%                If true, returns full Data struct with corrected U/V.
%                If false, returns only metadata/offsets.
%   plot_data  - Optional, default is true
%                If true, generates diagnostic figure.

    % 1. Set Defaults for Optional Arguments (Most consistent method)
    if nargin < 3 || isempty(apl_cor),  apl_cor = false; end
    if nargin < 4 || isempty(plot_data), plot_data = true; end

    % 2. Initialize Outputs Properly
    h_fig = [];
    if apl_cor
        Data_corr = Data; % Initialize with full dataset
    else
        Data_corr = struct(); % Initialize as empty metadata struct
    end

    % 3. Flexible Time Check (Handles 'Average_Time' vs 'time')
    if isfield(Data, 'Average_Time'),     T = Data.Average_Time;
    elseif isfield(Data, 'time'),         T = Data.time;
    else, error('Structure "Data" must contain a "time" or "Average_Time" field.'); end

    % 4. Extract and Clean Magnetometer/Velocity Data
    mx = Data.Average_Magnetometer(:,1); 
    my = Data.Average_Magnetometer(:,2); 
    U  = Data.Average_VelEast;
    V  = Data.Average_VelNorth;
    
    mask_QC_3D = Data.mask_QC_3D;
    mask_QC_1D = Data.mask_QC_1D;
    
    % Masking
    mx(mask_QC_1D ~= flag_value) = NaN;
    my(mask_QC_1D ~= flag_value) = NaN;
    U(mask_QC_3D ~= flag_value)  = NaN;
    V(mask_QC_3D ~= flag_value)  = NaN;
    
    %% 5. Solve for circle parameters using Least Squares
    valid = ~isnan(mx) & ~isnan(my);
    
    if sum(valid) > 3 
        A = [mx(valid), my(valid), ones(sum(valid), 1)];
        B = mx(valid).^2 + my(valid).^2;
        params = A \ B;
        xc = params(1)/2;
        yc = params(2)/2;
        radius = sqrt(params(3) + xc^2 + yc^2);
    else
        warning('Insufficient valid magnetometer data for circle fitting.');
        % Ensure consistent output even on failure
        Data_corr.hard_iron_CCW_angle = NaN;
        if apl_cor
            Data_corr.U_hard_iron_corrected = U * NaN;
            Data_corr.V_hard_iron_corrected = V * NaN;
        end
        return;
    end
    
    %% 6. Calculate the CCW Correction Angle (phi)
    mx_c = mx - xc;
    my_c = my - yc;
    angle_raw  = atan2d(my, mx);
    angle_corr = atan2d(my_c, mx_c);
    angle_err  = angle_corr - angle_raw; 
    
    % Normalize to [-180, 180]
    angle_err = atan2d(sind(angle_err), cosd(angle_err));
    
    %% 7. Apply Rotation and Store Metadata
    if apl_cor
        Data_corr.U_hard_iron_corrected = U .* cosd(angle_err) - V .* sind(angle_err);
        Data_corr.V_hard_iron_corrected = U .* sind(angle_err) + V .* cosd(angle_err);
    end 
    
    % Metadata is always included in Data_corr regardless of apl_cor flag
    Data_corr.hard_iron_offset_x  = xc;
    Data_corr.hard_iron_offset_y  = yc;
    Data_corr.radius              = radius;
    Data_corr.hard_iron_CCW_angle = angle_err;

    %% 8. Diagnostic Plotting
    if plot_data
        % Command Window Summary
        err_min = min(angle_err, [], 'omitnan');
        err_max = max(angle_err, [], 'omitnan');
        deg_symbol = char(176); 
        fprintf(1, ['\n\n*** Hard-Iron Compass Correction ***\n' ...
            'Circle fit error angle varies between %3.0f%s and %3.0f%s\n' ...
            'Offsets: X = %.1f, Y = %.1f\n'], ...
            err_min, deg_symbol, err_max, deg_symbol, xc, yc);

        % Plotting
        h_fig = figure;
        hold on;
        
        theta_vec = linspace(0, 2*pi, 300);
        avg_radius = mean(sqrt(mx.^2 + my.^2), 'omitnan');

        scatter(mx, my, 5, T, 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Mag Data');
        plot(avg_radius*cos(theta_vec), avg_radius*sin(theta_vec), 'r--', 'LineWidth', 1, 'DisplayName', 'Ideal Path (at 0,0)');
        plot(xc + radius*cos(theta_vec), yc + radius*sin(theta_vec), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Fitted Circle');
        plot(xc, yc, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Fit Center');

        axis equal; grid on;
        xlabel('Mag X (mG)'); ylabel('Mag Y (mG)');
        title('Horizontal Magnetometer Fit (Hard-Iron)');

        xl = xlim; yl = ylim;
        line(xl, [0 0], 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        line([0 0], yl, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');

        hcb = colorbar;
        ylabel(hcb, 'Deployment Time');
        datetick(hcb, 'y', 'yyyy-mmm', 'keepticks');

        legend('Location', 'best');
        text(xc, yc-10, sprintf(['Offset X: %.1f mG\nOffset Y: %.1f mG\n' ...
            'Radius: %.1f mG'], xc, yc, radius), ...
            'VerticalAlignment', 'top','HorizontalAlignment','center', ...
            'FontWeight', 'bold');
        hold off;
    end
end
