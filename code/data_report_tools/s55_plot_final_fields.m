% quick plot of post-processed S55 ADCP file
% you have to set path
%
% TODO
% Revise path definition to work with global MOORPROC_G
%
% KB - 20.05.2026

indir = "C:\Users\sa07kb\Projects\Moor_Data_Proc\moor_examples\osnap\data\moor\proc\rhadcp_01_2020\adcp_S55\";
infile = "rhadcp_01_2020_200044_stage2.nc";
filename = fullfile(indir,infile);
S55_final = s55_load_nc_as_struct(filename);

% 1. Get all field names inside your structure
all_fields = fieldnames(S55_final);

% 2. Initialize an empty cell array to hold the names of 2D variables
twod_variables = {};

% 3. Loop through each field and check its size
for i = 1:numel(all_fields)
    field_name = all_fields{i};
    field_size = size(S55_final.(field_name));
    
    % Check if the variable is a 2D matrix matching your Time x Depth dimensions
    if isequal(field_size, [30863, 56])
        disp(field_name)
        % 1. Choose which 2D variable you want to plot
        variable_to_plot = S55_final.(field_name); % Swap to VCUR, AMP_BEAM1, etc. as needed
        variable_name    = S55_final.([field_name,'_attrs']).long_name; 
        variable_units   = S55_final.([field_name,'_attrs']).units;

        % 2. Create the figure window
        figure('Color', 'w', 'Position', [100, 100, 1000, 500]);

        % 3. Create the 2D surface plot
        % Transposing the variable matrix ensures it aligns with TIME (X) and DEPTH (Y)
        h = pcolor(S55_final.TIME, S55_final.DEPTH, variable_to_plot');

        % 4. Format the rendering for clean visualization
        set(h, 'EdgeColor', 'none'); % Removes ugly black grid lines around 30k bins
        shading interp;             % Smooths out the color transitions
        axis tight;                 % Snaps the plot limits tightly to your data boundaries

        % 5. Handle the Orientation (Y-Axis)
        % Standard oceanographic profiles should show the surface (0m) at the top
        set(gca, 'YDir', 'reverse'); 

        % 6. Format the X-Axis for clear Date/Time labels
        datetick('x', 'yyyy-mm', 'keeplimits'); 

        % 7. Add Labels and Colorbar
        xlabel('Time (Year-Month)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Depth (m)', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s - %s - %s', field_name, S55_final.Attributes.platform_code, variable_name), 'FontSize', 12,'Interpreter', 'none');

        % Add colorbar and scale its label
        cb = colorbar;
        ylabel(cb, sprintf('[%s]', variable_units), 'FontSize', 11);

        % Optional: Set a symmetric divergent color map for current velocity vectors
        if ismember(field_name, {'UCUR', 'VCUR', 'WCUR'})
            colormap(cmocean('balance')); % Requires the Climate Data Toolbox/cmocean
            % If you do not have cmocean installed, uncomment the line below for standard MATLAB:
            % colormap(bluewhitered); 
        end
    end
end


