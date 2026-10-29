function [ out ] = tbpJ2_J6_accelerometerGating_rpod_lidarBias( t, ic, gravParam, beta, accelerometer, n_perSat, gate, Rbody_N, RN_rot, bias_body, Ccoeff, Scoeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Sets up the Equations of Motion (EOM) to integrate with
% spherical harmonics with accelerometer inputs when the accelerometer
% reads over a gate level
%
% Inputs:
%       t             time [s]
%       ic            initial conditions:
%                       position vector [km], velocity vector [km / s]
%       gravParam     gravitational parameter for the central body [km^3 /
%                       s^2] (mu, R, J2)
%       accelerometer accelerometer measurmets in body frame [m/s^2]
%       gate          gate for including accelerometer readings or not 
%       Rbody_N       DCM that describes rotation from inertial to body
%                       frame, 3D array where each layer is the DCM for
%                       that body 
%
% Outputs:
%       out           current state:
%                       postion vector [m], velocity vector [m / s]
%
% Author: Rachel Mamich
% Created: 11/10/18
% Modified: 12/21/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(ic);
out = zeros(N, 1);
mu = gravParam(1);
R = gravParam(2);
J2 = gravParam(3);
J3 = gravParam(4);
J4 = gravParam(5);
J5 = gravParam(6);
J6 = gravParam(7);

%% Target
rx_t = ic(1);
ry_t = ic(2);
rz_t = ic(3);
vx_t = ic(4);
vy_t = ic(5);
vz_t = ic(6);

if n_perSat(1) == 9
    % accelerometer bias in body frame
    burn_t = ic(7:9);
else
    burn_t = zeros(3, 1);
end
       
r_t = [rx_t; ry_t; rz_t];
a_t = - mu / (r_t' * r_t)^(3/2) * r_t;
r_tNorm = (rx_t^2 + ry_t^2 + rz_t^2)^(1/2);
r_texp5 = r_tNorm^5;
r_texp7 = r_tNorm^7;
r_texp13 = r_tNorm^13;
r_texp15 = r_tNorm^15;

    
a_tJ2x = 3/2 * J2 * mu * R^2 / r_texp5 * (5 * (r_t(3) / (r_tNorm))^2 - 1) * r_t(1);
a_tJ2y = 3/2 * J2 * mu * R^2 / r_texp5 * (5 * (r_t(3) / (r_tNorm))^2 - 1) * r_t(2);
a_tJ2z = 3/2 * J2 * mu * R^2 / r_texp5 * (5 * (r_t(3) / (r_tNorm))^2 - 3) * r_t(3);
a_tJ2 = [a_tJ2x; a_tJ2y; a_tJ2z];


a_tJ3x = -1/2 * J3 * mu * R^3 / r_texp7 * 5*(3 - 7*(r_t(3)/r_tNorm)^2) * r_t(1)*r_t(3);
a_tJ3y = -1/2 * J3 * mu * R^3 / r_texp7 * 5*(3 - 7*(r_t(3)/r_tNorm)^2) * r_t(2)*r_t(3);
a_tJ3z = 1/2 * J3 * mu * R^3 / r_texp5 * (3 - 30*(r_t(3)/r_tNorm)^2 + 35*(r_t(3)/r_tNorm)^4);
a_tJ3 = [a_tJ3x; a_tJ3y; a_tJ3z];


a_tJ4x = 5/8 * J4 * mu * R^4 / r_texp7 *(3 - 42*(r_t(3)/r_tNorm)^2 + 63*(r_t(3)/r_tNorm)^4)*r_t(1);
a_tJ4y = 5/8 * J4 * mu * R^4 / r_texp7 *(3 - 42*(r_t(3)/r_tNorm)^2 + 63*(r_t(3)/r_tNorm)^4)*r_t(2);
a_tJ4z = 5/8 * J4 * mu * R^4 / r_texp7 *(15 - 70*(r_t(3)/r_tNorm)^2 + 63*(r_t(3)/r_tNorm)^4)*r_t(3);
a_tJ4 = [a_tJ4x; a_tJ4y; a_tJ4z];



a_tJ5x = 21/8 * J5 * mu * R^5 / r_texp13 *(5*r_t(1)^4 + 10*r_t(1)^2*r_t(2)^2 - 20*r_t(1)^2*r_t(3)^2 + 5*r_t(2)^4 - 20*r_t(2)^2*r_t(3)^2 + 8*r_t(3)^4)*r_t(1)*r_t(3);
a_tJ5y = 21/8 * J5 * mu * R^5 / r_texp13 *(5*r_t(1)^4 + 10*r_t(1)^2*r_t(2)^2 - 20*r_t(1)^2*r_t(3)^2 + 5*r_t(2)^4 - 20*r_t(2)^2*r_t(3)^2 + 8*r_t(3)^4)*r_t(2)*r_t(3);
a_tJ5z = -3/8 * J5 * mu * R^5 / r_texp13 *(5*r_t(1)^6 + 15*r_t(1)^4*r_t(2)^2 - 90*r_t(1)^4*r_t(3)^2 + 15*r_t(1)^2*r_t(2)^4 - 180*r_t(1)^2*r_t(2)^2*r_t(3)^2 + 120*r_t(1)^2*r_t(3)^4 + 5*r_t(2)^6 - 90*r_t(2)^4*r_t(3)^2 + 120*r_t(2)^2*r_t(3)^4 - 16*r_t(3)^6);
a_tJ5 = [a_tJ5x; a_tJ5y; a_tJ5z];


arg_tJ6 = 5*r_t(1)^6 + 15*r_t(1)^4*r_t(2)^2 - 120*r_t(1)^4*r_t(3)^2 + 15*r_t(1)^2*r_t(2)^4 - 240*r_t(1)^2*r_t(2)^2*r_t(3)^2 + 240*r_t(1)^2*r_t(3)^4 + 5*r_t(2)^6 - 120*r_t(2)^4*r_t(3)^2 + 240*r_t(2)^2*r_t(3)^4 - 64*r_t(3)^6;
a_tJ6x = -7/16 * J6 * mu * R^6 / r_texp15*arg_tJ6*r_t(1);
a_tJ6y = -7/16 * J6 * mu * R^6 / r_texp15*arg_tJ6*r_t(2);
a_tJ6z = -7/16 * J6 * mu * R^6 / r_texp15*(35*r_t(1)^6 + 105*r_t(1)^4*r_t(2)^2 - 210*r_t(1)^4*r_t(3)^2 + 105*r_t(1)^2*r_t(2)^4 - 420*r_t(1)^2*r_t(2)^2*r_t(3)^2 + 168*r_t(1)^2*r_t(3)^4 + 35*r_t(2)^6 - 210*r_t(2)^4*r_t(3)^2 + 168*r_t(2)^2*r_t(3)^4 - 16*r_t(3)^6)*r_t(3);
a_tJ6 = [a_tJ6x; a_tJ6y; a_tJ6z];




a_t = a_t + burn_t + a_tJ2 + a_tJ3 + a_tJ4 + a_tJ5 + a_tJ6;

if n_perSat(1) == 9
    out(1:9) = [vx_t; vy_t; vz_t; a_t(1); a_t(2); a_t(3); -beta(1)*burn_t];
else
    out(1:6) = [vx_t; vy_t; vz_t; a_t(1); a_t(2); a_t(3)];
end




%% Chaser 
rx_c = ic(n_perSat(1) + 1);
ry_c = ic(n_perSat(1) + 2);
rz_c = ic(n_perSat(1) + 3);
vx_c = ic(n_perSat(1) + 4);
vy_c = ic(n_perSat(1) + 5);
vz_c = ic(n_perSat(1) + 6);

if n_perSat(2) == 9
    % accelerometer bias in body frame
    bias_xbody = ic(n_perSat(1) + 7);
    bias_ybody = ic(n_perSat(1) + 8);
    bias_zbody = ic(n_perSat(1) + 9);
    bias_bodyi = [bias_xbody; bias_ybody; bias_zbody];
end

r_c = [rx_c; ry_c; rz_c];
a_c_TBP = - mu / (r_c' * r_c)^(3/2) * r_c;
r_cNorm = (rx_c^2 + ry_c^2 + rz_c^2)^(1/2);
r_cexp5 = r_cNorm^5;
r_cexp7 = r_cNorm^7;
r_cexp13 = r_cNorm^13;
r_cexp15 = r_cNorm^15;

    
a_cJ2x = 3/2 * J2 * mu * R^2 / r_cexp5 * (5 * (r_c(3) / (r_cNorm))^2 - 1) * r_c(1);
a_cJ2y = 3/2 * J2 * mu * R^2 / r_cexp5 * (5 * (r_c(3) / (r_cNorm))^2 - 1) * r_c(2);
a_cJ2z = 3/2 * J2 * mu * R^2 / r_cexp5 * (5 * (r_c(3) / (r_cNorm))^2 - 3) * r_c(3);
a_cJ2 = [a_cJ2x; a_cJ2y; a_cJ2z];


a_cJ3x = -1/2 * J3 * mu * R^3 / r_cexp7 * 5*(3 - 7*(r_c(3)/r_cNorm)^2) * r_c(1)*r_c(3);
a_cJ3y = -1/2 * J3 * mu * R^3 / r_cexp7 * 5*(3 - 7*(r_c(3)/r_cNorm)^2) * r_c(2)*r_c(3);
a_cJ3z = 1/2 * J3 * mu * R^3 / r_cexp5 * (3 - 30*(r_c(3)/r_cNorm)^2 + 35*(r_c(3)/r_cNorm)^4);
a_cJ3 = [a_cJ3x; a_cJ3y; a_cJ3z];


a_cJ4x = 5/8 * J4 * mu * R^4 / r_cexp7 *(3 - 42*(r_c(3)/r_cNorm)^2 + 63*(r_c(3)/r_cNorm)^4)*r_c(1);
a_cJ4y = 5/8 * J4 * mu * R^4 / r_cexp7 *(3 - 42*(r_c(3)/r_cNorm)^2 + 63*(r_c(3)/r_cNorm)^4)*r_c(2);
a_cJ4z = 5/8 * J4 * mu * R^4 / r_cexp7 *(15 - 70*(r_c(3)/r_cNorm)^2 + 63*(r_c(3)/r_cNorm)^4)*r_c(3);
a_cJ4 = [a_cJ4x; a_cJ4y; a_cJ4z];



a_cJ5x = 21/8 * J5 * mu * R^5 / r_cexp13 *(5*r_c(1)^4 + 10*r_c(1)^2*r_c(2)^2 - 20*r_c(1)^2*r_c(3)^2 + 5*r_c(2)^4 - 20*r_c(2)^2*r_c(3)^2 + 8*r_c(3)^4)*r_c(1)*r_c(3);
a_cJ5y = 21/8 * J5 * mu * R^5 / r_cexp13 *(5*r_c(1)^4 + 10*r_c(1)^2*r_c(2)^2 - 20*r_c(1)^2*r_c(3)^2 + 5*r_c(2)^4 - 20*r_c(2)^2*r_c(3)^2 + 8*r_c(3)^4)*r_c(2)*r_c(3);
a_cJ5z = -3/8 * J5 * mu * R^5 / r_cexp13 *(5*r_c(1)^6 + 15*r_c(1)^4*r_c(2)^2 - 90*r_c(1)^4*r_c(3)^2 + 15*r_c(1)^2*r_c(2)^4 - 180*r_c(1)^2*r_c(2)^2*r_c(3)^2 + 120*r_c(1)^2*r_c(3)^4 + 5*r_c(2)^6 - 90*r_c(2)^4*r_c(3)^2 + 120*r_c(2)^2*r_c(3)^4 - 16*r_c(3)^6);
a_cJ5 = [a_cJ5x; a_cJ5y; a_cJ5z];


arg_cJ6 = 5*r_c(1)^6 + 15*r_c(1)^4*r_c(2)^2 - 120*r_c(1)^4*r_c(3)^2 + 15*r_c(1)^2*r_c(2)^4 - 240*r_c(1)^2*r_c(2)^2*r_c(3)^2 + 240*r_c(1)^2*r_c(3)^4 + 5*r_c(2)^6 - 120*r_c(2)^4*r_c(3)^2 + 240*r_c(2)^2*r_c(3)^4 - 64*r_c(3)^6;
a_cJ6x = -7/16 * J6 * mu * R^6 / r_cexp15*arg_cJ6*r_c(1);
a_cJ6y = -7/16 * J6 * mu * R^6 / r_cexp15*arg_cJ6*r_c(2);
a_cJ6z = -7/16 * J6 * mu * R^6 / r_cexp15*(35*r_c(1)^6 + 105*r_c(1)^4*r_c(2)^2 - 210*r_c(1)^4*r_c(3)^2 + 105*r_c(1)^2*r_c(2)^4 - 420*r_c(1)^2*r_c(2)^2*r_c(3)^2 + 168*r_c(1)^2*r_c(3)^4 + 35*r_c(2)^6 - 210*r_c(2)^4*r_c(3)^2 + 168*r_c(2)^2*r_c(3)^4 - 16*r_c(3)^6)*r_c(3);
a_cJ6 = [a_cJ6x; a_cJ6y; a_cJ6z];






    if n_perSat(2) == 9
        if norm(accelerometer) > gate
            burn_c = Rbody_N'*(accelerometer - bias_bodyi);
        else
            burn_c = zeros(3, 1);
        end
        a_c = a_c_TBP + burn_c + a_cJ2 + a_cJ3 + a_cJ4 + a_cJ5 + a_cJ6;
        out(n_perSat(1)+1:n_perSat(1)+n_perSat(2)) = [vx_c; vy_c; vz_c; a_c(1); a_c(2); a_c(3); -beta(2)*[bias_xbody; bias_ybody; bias_zbody]];
    else
        if nargin >= 10
            bias_bodyi = bias_body;
        else
            bias_bodyi = zeros(3, 1);
        end
        if norm(accelerometer) > gate
            burn_c = Rbody_N'*(accelerometer - bias_bodyi);
        else
            burn_c = zeros(3, 1);
        end
        a_c = a_c_TBP + burn_c + a_cJ2 + a_cJ3 + a_cJ4 + a_cJ5 + a_cJ6;
        out(n_perSat(1)+1:n_perSat(1)+n_perSat(2)) = [vx_c; vy_c; vz_c; a_c(1); a_c(2); a_c(3)];
    end
    
    if N > n_perSat(1)+ n_perSat(2)
        out(end) = -beta(3)*ic(end); % Lidar bias as Markov 
    end
    
end

