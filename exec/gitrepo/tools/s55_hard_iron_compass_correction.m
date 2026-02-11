function Data_corr = s55_hard_iron_compass_correction(Data,flag_value,apl_cor)
% S55_HARD_IRON_COMPASS_CORRECTION Corrects ENU velocities for magnetic interference.
%
% USAGE:
%   Data = s55_hard_iron_compass_correction(Data, QC_1D)
%
% DESCRIPTION:
%   This function calculates the hard-iron offset (center shift) of the 
%   ADCP's magnetometer data using a least-squares circle fit. It then 
%   determines the angular error (phi) introduced by this shift.
%
%   The correction is applied to the East (U) and North (V) velocities 
%   using a standard Counter-Clockwise (CCW) rotation matrix:
%       U_new = U*cos(phi) - V*sin(phi)
%       V_new = U*sin(phi) + V*cos(phi)
%
% INPUTS:
%   Data  - Structure containing Average_Magnetometer [Nx2], U, and V.
%   QC_1D - Quality control mask (0 = Good, others = NaN).
%
% OUTPUTS:
%   Data  - Structure with U_corrected, V_corrected, and hard_iron_CCW_angle.

    %% 1. Extract and Clean Magnetometer Data
    
    mx = Data.Average_Magnetometer(:,1); 
    my = Data.Average_Magnetometer(:,2); 
    U = Data.Average_VelEast;
    V = Data.Average_VelNorth;
    
    % Apply QC masking
    mask_QC_3D = Data.mask_QC_3D;
    mask_QC_1D = Data.mask_QC_1D;
    mx(mask_QC_1D ~= flag_value) = NaN;
    my(mask_QC_1D ~= flag_value) = NaN;

    U(mask_QC_3D ~= flag_value) = NaN;
    V(mask_QC_3D ~= flag_value) = NaN;
    
    %% 2. Solve for circle parameters using Least Squares
    valid = ~isnan(mx) & ~isnan(my);
    if apl_cor
        Data_corr = Data;
    end
    
    if sum(valid) > 3 
        % Linearized circle equation: (x-xc)^2 + (y-yc)^2 = R^2
        A = [mx(valid), my(valid), ones(sum(valid), 1)];
        B = mx(valid).^2 + my(valid).^2;
        
        params = A \ B;
        xc = params(1)/2;
        yc = params(2)/2;
        radius = sqrt(params(3) + xc^2 + yc^2);
    else
        warning('Insufficient valid magnetometer data for circle fitting.');
        Data_corr.U_corrected = U*NaN;
        Data_corr.V_corrected = V*NaN;
        Data_corr.hard_iron_CCW_angle = NaN;
        return;
    end
    
    %% 3. Calculate the CCW Correction Angle (phi)
    % Shifted coordinates
    mx_c = mx - xc;
    my_c = my - yc;
    
    % Calculate angles using Cartesian convention atan2(y, x)
    angle_raw  = atan2d(my, mx);
    angle_corr = atan2d(my_c, mx_c);
    
    % phi is the CCW angle needed to rotate the raw vectors to the corrected state
    angle_err = angle_corr - angle_raw; 
    
    % Normalize phi to [-180, 180] to prevent wrapping jumps (e.g., 359 to 1)
    angle_err = atan2d(sind(angle_err), cosd(angle_err));
    
    %% 4. Apply Rotation to ENU Velocities
  
    % Apply the CCW rotation matrix
    % Note: If U/V are [Time x Bins], phi [Time x 1] expands automatically
    Data_corr.U_hard_iron_corrected = U .* cosd(angle_err) - V .* sind(angle_err);
    Data_corr.V_hard_iron_corrected = U .* sind(angle_err) + V .* cosd(angle_err);

    
    %% 5. Store Metadata
    Data_corr.hard_iron_offset_x = xc;
    Data_corr.hard_iron_offset_y = yc;
    Data_corr.hard_iron_CCW_angle = angle_err;

    if apl_cor==false
        Data_corr.radius = radius;
        Data_corr.params = params;
    end
    
end