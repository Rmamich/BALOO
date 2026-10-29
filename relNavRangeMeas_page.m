function [out, H] = relNavRangeMeas_page(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: calculates the range and range rate of the target wrt the chaser
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
%
% Outputs:
%       rho     range [km] 
%       rhoDot  range rate [km/s]
%
% Author: Rachel Mamich
% Created: 09/09/22
% Modified: 04/11/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N, ~, n_nodes] = size(x);
rr = x(1:3, :, :);
lidarBias = x(19, :, :);

rhoVec = rr;
rho = pagenorm(rhoVec(1:3, :, :)); 


out = rho + lidarBias;

H = zeros(1, N, n_nodes);
H(1, 1:3, :) = pagemtimes(rhoVec, 'transpose', 1./rho, 'none');
H(1, end, :) = 1;




end

