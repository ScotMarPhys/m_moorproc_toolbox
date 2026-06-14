function [Data_corr, h_fig] = s55_hard_iron_compass_correction(Data, flag_value, apl_cor, plot_data, select_time_frame,outdir,filename)
    % 1. Set Defaults
    if nargin < 3 || isempty(apl_cor),  apl_cor = false; end
    if nargin < 4 || isempty(plot_data), plot_data = true; end
    if nargin < 5 || isempty(select_time_frame), select_time_frame = false; end
    if nargin < 7 || isempty(outdir) || isempty(filename)
        select_time_file = false;
    else
        select_time_file = true;
    end

    % 2. Initialize Outputs
    h_fig = [];
    if apl_cor, Data_corr = Data; else, Data_corr = struct(); end

    % 3. Time Check
    if isfield(Data, 'Average_Time'), T = Data.Average_Time;
    elseif isfield(Data, 'time'), T = Data.time;
    else, error('No time field found.'); end

    % 4. Extract and Clean Data
    mx = Data.Average_Magnetometer(:,1); 
    my = Data.Average_Magnetometer(:,2); 
    U  = Data.Average_VelEast;
    V  = Data.Average_VelNorth;
    mask_QC_time = Data.mask_QC_time;
    mask_QC_2D = Data.mask_QC_2D;
    
    mx(mask_QC_time ~= flag_value) = NaN;
    my(mask_QC_time ~= flag_value) = NaN;
    U(mask_QC_2D ~= flag_value)  = NaN;
    V(mask_QC_2D ~= flag_value)  = NaN;

    % Variables for iterative fitting
    keep_refining = true;

    if select_time_file
        csv_path = fullfile(outdir, [filename, '_hard_iron_time_selection.csv']);
        if exist(csv_path, 'file')
            pre_selected_limits = readmatrix(csv_path);
            if ~isempty(pre_selected_limits)
                time_mask = (T >= pre_selected_limits(1)) & (T <= pre_selected_limits(2));
            end
        else
            time_mask = true(size(T)); % Initially use all data
        end
    else
        time_mask = true(size(T)); % Initially use all data
    end
        

    while keep_refining
        %% 5. Solve for circle parameters using valid data in the SELECTED window
        idx_fit = ~isnan(mx) & ~isnan(my) & time_mask;
        
        if sum(idx_fit) > 3
            % Subset the data to remove NaNs
            Xi = mx(idx_fit);
            Yi = my(idx_fit);
            
            % Initial guess: [center_x, center_y, radius]
            p0 = [mean(Xi), mean(Yi), std(Xi)];
            
            % Objective function: (x-xc)^2 + (y-yc)^2 - r^2
            % Ensure this function only sees the Xi/Yi subset
            fun = @(p) (Xi - p(1)).^2 + (Yi - p(2)).^2 - p(3)^2;
            
            % Optimization Options (to suppress iterations in command window)
            opts = optimoptions('lsqnonlin','Display','off');
            
            try
                p = lsqnonlin(fun, p0, [], [], opts);
                xc = p(1); yc = p(2); radius = p(3);
            catch
                % Fallback to Linear if lsqnonlin fails
                A = [Xi, Yi, ones(size(Xi))]; B = Xi.^2 + Yi.^2;
                p_lin = A \ B;
                xc = p_lin(1)/2; yc = p_lin(2)/2;
                radius = sqrt(p_lin(3) + xc^2 + yc^2);
            end
        else
            error('Insufficient valid data points for circle fitting.');
        end
        
        %% 6. Calculate Error Angle (phi) for the WHOLE dataset
        mx_c = mx - xc;
        my_c = my - yc;
        angle_raw  = atan2d(my, mx);
        angle_corr = atan2d(my_c, mx_c);
        angle_err  = atan2d(sind(angle_corr - angle_raw), cosd(angle_corr - angle_raw));
        
        %% 7. Diagnostic Plotting
        if plot_data
            h_fig = figure; clf; hold on;
            theta_vec = linspace(0, 2*pi, 300);
            avg_radius = mean(sqrt(mx.^2 + my.^2), 'omitnan');

            scatter(mx, my, 5, T, 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'All Data');
            if select_time_frame && any(~time_mask)
                scatter(mx(idx_fit), my(idx_fit), 10, 'k', 'DisplayName', 'Data used for Fit');
            end
            plot(avg_radius*cos(theta_vec), avg_radius*sin(theta_vec), 'r--', 'DisplayName', 'Ideal');
            plot(xc + radius*cos(theta_vec), yc + radius*sin(theta_vec), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Fitted Circle');
            plot(xc, yc, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Fit Center');
            
            axis equal; grid on; 
            hcb = colorbar; 
            datetick(hcb, 'y', 'yyyy-mmm', 'keepticks');
            
            title('Hard-Iron Circle Fit'); legend('Location','best');
            hold off;
            
            % ONLY trigger refinement if select_time_frame is true
            if select_time_frame
                choice = input('\nRefine time frame for calibration? (y/n) [n]: ', 's');
                if strcmpi(choice, 'y')
                    h_sel = figure(6); clf;
                    ax(1) = subplot(2,1,1);
                    plot(T, Data.Average_Heading); 
                    title('Heading'); grid on; datetick('x','keeplimits');
                    
                    ax(2) = subplot(2,1,2);
                    plot(T, mx, T, my);
                    title('Magnetometer'); grid on; ylabel('mG'); datetick('x','keeplimits');
                    
                    linkaxes(ax, 'x');
                    fprintf('--> ZOOM to window, then press ENTER in command window.\n');
                    pause; 
                    
                    limits = xlim(ax(2));
                    time_mask = (T >= limits(1)) & (T <= limits(2));
                    fprintf('Selected: %s to %s\n', datestr(limits(1)), datestr(limits(2)));
                    close(h_sel);
                    % SAVE the new selection
                    if ~all(time_mask)
                        writematrix(limits, csv_path);
                    end
                else
                    keep_refining = false;
                end
            else
                keep_refining = false;
            end
        else
            keep_refining = false;
        end
    end

    %% 8. Final Application and Metadata
    if apl_cor
        Data_corr.U_hard_iron_corrected = U .* cosd(angle_err) - V .* sind(angle_err);
        Data_corr.V_hard_iron_corrected = U .* sind(angle_err) + V .* cosd(angle_err);
    end 
    
    Data_corr.hard_iron_offset_x  = xc;
    Data_corr.hard_iron_offset_y  = yc;
    Data_corr.hard_iron_radius    = radius;
    Data_corr.hard_iron_CCW_angle = angle_err;
    Data_corr.hard_iron_time_window  = [T(find(time_mask,1)) T(find(time_mask,1,'last'))];
end
