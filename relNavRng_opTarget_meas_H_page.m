function [out, H] = relNavRng_opTarget_meas_H_page(x, RCF1_N, n_elm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Collects range, range rate, relative optical, and ground feature
% measuremnts for the UKF to digest.
%
% Inputs:
%       x       state vector
%                   x(1:3)     target position relative to chaser position 
%                               in the inertial frame [km] 
%                   x(4:6)     target velocity relative to chaser velocity 
%                               in the inertial frame [km/s]
%                   x(7:9)     chaser unmodeled acceleration in inertial 
%                               frame [km/s^2]
%                   x(10:12)   chaser position in the inertial frame [km]
%                   x(13:15)   chaser velocity in the inertial frame [km/s]
%                   x(16:18)   unmodeled accelerometer bias est [km/s^2]
%                   x(19)      lidar bias
%       RCF1_N  DCM from inertial frame to camera frame directed toward
%                           target vehicle
%       n_elm   number of elements for realtive target states and chaser
%                           states 
%       
%
% Outputs:
%       out     measurement vector
%                   out(1)      range
%                   out(2)      range rate
%                   out(3:4)    A, B: azimuth and elevation of the target
%                                   vehicle in the camera frame  [rad, rad]
%
% Author: Rachel Mamich
% Created: 09/16/22
% Modified: 04/11/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n_xOriginal, n_nodes] = size(x);

if n_xOriginal < 19
    % Need to fill
    if n_elm(1) < 9
        % If we are missing the acceleration terms
        x = [x(1:6, :); zeros(3, n_nodes); x(7:end, :)];
    end
    
    if n_elm(2) < 9
        % If we are missing the bias terms
        if length(x) == 15
            % Also missing bias 
            % add filler for both accelerometer bias and lidar bias
            x = [x; zeros(4, n_nodes)];
        elseif length(x) == 16
            % not missing bias
            % Add filler for accelerometer bias while keeping lidar bias at the
            % end
            x = [x(1:15, :); zeros(3, n_nodes); x(end, :)];
        end
    else
        % accel bias terms present, but no lidar bias 
        x = [x; zeros(1, n_nodes)];
    end
end

x = reshape(x, [], 1, n_nodes);

[rng, H_rng] = relNavRangeMeas_page(x);
[azEl_target, H_target] = relNavPinHoleOpticalMeas_target_page(x, RCF1_N);

out = [rng; azEl_target];
H = [H_rng; H_target];



if n_xOriginal < 19
    if n_elm(2) < 9
        % Acceleration bias not included 
        if n_elm(1) + n_elm(2) < n_xOriginal
            % Lidar bias here, keep last element
            % Indicated we have a lidar bias present
            H = H(:, [1:15, 19], :);
        else
            % No Lidar bias, no acceleration bias
            % Indicates we have no lidar bias present
            H = H(:, 1:15, :);
        end
    else
        % Acceleration bias present
        if n_elm(1) + n_elm(2) == n_xOriginal
            % Indicates we have no lidar bias present
            H = H(:, 1:18, :);
        end
    end
    
    if n_elm(1) < 9
        % No acceleration states
        [~, n_H] = size(H);
        H = H(:, [1:6, 10:n_H], :);
    end
end

end

