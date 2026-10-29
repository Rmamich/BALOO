function [out, H] = relNavDirectInertial_meas_H_page(x)
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

out = x(1:3, :);
H = zeros(3, n_xOriginal, n_nodes);
H(1:3, 1:3, :) = repmat(eye(3), 1, 1, n_nodes);

end

