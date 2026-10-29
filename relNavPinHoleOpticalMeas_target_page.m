function [out, H] = relNavPinHoleOpticalMeas_target_page(x, RCF_N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: calculates the optical measurment generated from observing the 
% "target" or "chief" vehicle as seen by the "deputy" or "chaser" vehicle.
%
% Inputs:
%       x       state vector
%                   x(1:3)     chief position in the inertial frame [m] 
%                   x(4:6)     chief velocity in the inertial frame [m/s]
%                   x(7:9)     chief unmodeled acceleration [m/s^2]
%                   x(10:12)   deputy position in the inertial frame [m]
%                   x(13:15)   deputy velocity in the inertial frame [m/s]
%                   x(16:18)   deputy unmodeled acceleration est [m/s^2]
%                   x(19)      lidar bias [m]
%       
%
% Outputs:
%       out     A, B: azimuth and elevation in the camera frame  [rad, rad]    
%
% Author: Rachel Mamich
% Created: 09/09/22
% Modified: 12/21/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N, ~, n_nodes] = size(x);
rr = x(1:3, :, :);

rho_N = rr;
rho_CF = pagemtimes(RCF_N, rho_N);

A = atan2(rho_CF(1, :, :), rho_CF(3, :, :));
B = atan2(rho_CF(2, :, :), rho_CF(3, :, :));

out = [A; B];

H = zeros(2, N, n_nodes);

H(1,1:3, :) = pagemtimes(pagemtimes([rho_CF(3, :, :), zeros(1, 1, n_nodes), -rho_CF(1, :, :)], RCF_N), 1 ./ (rho_CF(3, :, :).^2 + rho_CF(1, :, :).^2));
H(2,1:3, :) = pagemtimes(pagemtimes([zeros(1, 1, n_nodes), rho_CF(3, :, :) -rho_CF(2, :, :)], RCF_N),  1 ./ (rho_CF(3, :, :).^2 + rho_CF(2, :, :).^2));


end

