function [out, H] = rng_opTarget_meas_H_page(x, RCF1_N, accelFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Collects range, range rate, relative optical, and ground feature
% measuremnts for the UKF to digest.
%
% Inputs:
%       x       state vector
%                   x(1:3)     target position in the inertial frame [km] 
%                   x(4:6)     target velocity in the inertial frame [km/s]
%                   x(7:9)     target unmodeled acceleration [km/s^2]
%                   x(10:12)   servicer position in the target lvlh frame [km]
%                   x(13:15)   servicer velocity in the target lvlh frame [km/s]
%                   x(16:18)   servicer unmodeled acceleration est [km/s^2]
%       accelFlag   if unmodeled acceleration is included, true if it is
%       frameConvertFlag    flag to convert to the inertial frame if needed
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
% Modified: 10/24/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~accelFlag
    x = [x(1:6); zeros(3, 1); x(7:12); zeros(3, 1)];
end

[rngRngRate, H_rng] = rangeMeas(x);
[azEl_target, H_target] = pinHoleOpticalMeas_target(x, RCF1_N);

out = [rngRngRate(1); azEl_target];
H = [H_rng(1, :); H_target];

if ~accelFlag
    H = H(:, [1:6, 10:15]);
end

end

