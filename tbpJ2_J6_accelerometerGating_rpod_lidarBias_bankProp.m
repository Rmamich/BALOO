function [ out ] = tbpJ2_J6_accelerometerGating_rpod_lidarBias_bankProp( t, ic, gravParam, beta, accelerometer, n_perSat, gate, Rbody_N, RN_rot, bias_body, Ccoeff, Scoeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Sets up the Equations of Motion (EOM) to integrate with
% spherical harmonics with accelerometer inputs when the accelerometer
% reads over a gate level
%
% Inputs:
%       t             time [s]
%       ic            initial conditions:
%                       target position vector [m], target velocity vector
%                       [m / s], optional unmodeled target acceleration
%                       [m /s^2], chaser position vector [m], chaser 
%                       velocity vector [m / s], optional chaser 
%                       accelerometer bias [m/s^2]
%       gravParam     gravitational parameter for the central body [m^3 /
%                       s^2] (mu, R, J2, J3, J4, J5, J6)
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
% Modified: 03/25/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N, n_bank] = size(ic);
out = zeros(N, n_bank);
mu = gravParam(1);
R = gravParam(2);
J2 = gravParam(3);
J3 = gravParam(4);
J4 = gravParam(5);
J5 = gravParam(6);
J6 = gravParam(7);
ic = reshape(ic, N, 1, n_bank);


%% Target
rx_t = ic(1, :, :);
ry_t = ic(2, :, :);
rz_t = ic(3, :, :);


if n_perSat(1) == 9
    % accelerometer bias in body frame
    burn_t = ic(7:9, :, :);
else
    burn_t = zeros(3, 1, n_bank);
end

r_t = ic(1:3, :, :);
v_t = ic(4:6, :, :);

r_tNorm =  pagenorm(r_t);
r_texp5 = r_tNorm.^5;
r_texp7 = r_tNorm.^7;
r_texp13 = r_tNorm.^13;
r_texp15 = r_tNorm.^15;


rx_t2 = rx_t.^2;
rx_t4 = rx_t.^4;
rx_t6 = rx_t.^6;


ry_t2 = ry_t.^2;
ry_t4 = ry_t.^4;
ry_t6 = ry_t.^6;


rz_t2 = rz_t.^2;
rz_t4 = rz_t.^4;
rz_t6 = rz_t.^6;


a_t = - pagemtimes(mu  / (r_tNorm).^(3), r_t);

a_tJ2x = pagemtimes((3/2 * J2 * mu * R^2 ./ r_texp5), pagemtimes((5 * pagemrdivide(rz_t, r_tNorm).^2 - 1), rx_t));
a_tJ2y = pagemtimes((3/2 * J2 * mu * R^2 ./ r_texp5), pagemtimes((5 * pagemrdivide(rz_t, r_tNorm).^2 - 1), ry_t));
a_tJ2z = pagemtimes((3/2 * J2 * mu * R^2 ./ r_texp5), pagemtimes((5 * pagemrdivide(rz_t, r_tNorm).^2 - 3), rz_t));
a_tJ2 = [a_tJ2x; a_tJ2y; a_tJ2z];


a_tJ3x = pagemtimes((-1/2 * J3 * mu * R^3 ./ r_texp7),  pagemtimes(pagemtimes(5*(3 - 7*pagemrdivide(rz_t, r_tNorm).^2),  rx_t), rz_t));
a_tJ3y = pagemtimes((-1/2 * J3 * mu * R^3 ./ r_texp7),  pagemtimes(pagemtimes(5*(3 - 7*pagemrdivide(rz_t, r_tNorm).^2),  ry_t), rz_t));
a_tJ3z = pagemtimes((1/2 * J3 * mu * R^3 ./ r_texp5),  (3 - 30*pagemrdivide(rz_t, r_tNorm).^2 + 35*pagemrdivide(rz_t, r_tNorm).^4));
a_tJ3 = [a_tJ3x; a_tJ3y; a_tJ3z];


a_tJ4x = pagemtimes((5/8 * J4 * mu * R^4 ./ r_texp7), pagemtimes((3 - 42*pagemrdivide(rz_t, r_tNorm).^2 + 63*pagemrdivide(rz_t, r_tNorm).^4), rx_t));
a_tJ4y = pagemtimes((5/8 * J4 * mu * R^4 ./ r_texp7), pagemtimes((3 - 42*pagemrdivide(rz_t, r_tNorm).^2 + 63*pagemrdivide(rz_t, r_tNorm).^4), ry_t));
a_tJ4z = pagemtimes((5/8 * J4 * mu * R^4 ./ r_texp7), pagemtimes((15 - 70*pagemrdivide(rz_t, r_tNorm).^2 + 63*pagemrdivide(rz_t, r_tNorm).^4), rz_t));
a_tJ4 = [a_tJ4x; a_tJ4y; a_tJ4z];



a_tJ5x = pagemtimes((21/8 * J5 * mu * R^5 ./ r_texp13), pagemtimes(pagemtimes((5*rx_t4 + 10*pagemtimes(rx_t2, ry_t2) - 20*pagemtimes(rx_t2, rz_t2) + 5*ry_t4 - 20*pagemtimes(ry_t2,rz_t2) + 8*rz_t4), rx_t), rz_t));
a_tJ5y = pagemtimes((21/8 * J5 * mu * R^5 ./ r_texp13), pagemtimes(pagemtimes((5*rx_t4 + 10*pagemtimes(rx_t2, ry_t2) - 20*pagemtimes(rx_t2, rz_t2) + 5*ry_t4 - 20*pagemtimes(ry_t2,rz_t2) + 8*rz_t4), ry_t), rz_t));
a_tJ5z = pagemtimes((-3/8 * J5 * mu * R^5 ./ r_texp13), (5*rx_t6 + 15*pagemtimes(rx_t4, ry_t2) - 90*pagemtimes(rx_t4, rz_t2) + 15*pagemtimes(rx_t2,ry_t4) - 180*pagemtimes(pagemtimes(rx_t2, ry_t2), rz_t2) + 120*pagemtimes(rx_t2, rz_t4) + 5*ry_t6 - 90*pagemtimes(ry_t4, rz_t2) + 120*pagemtimes(ry_t2, rz_t4) - 16*rz_t6));
a_tJ5 = [a_tJ5x; a_tJ5y; a_tJ5z];


argJ6_t = 5*rx_t6 + 15*pagemtimes(rx_t4,ry_t2) - 120*pagemtimes(rx_t4,rz_t2) + 15*pagemtimes(rx_t2, ry_t4) - 240*pagemtimes(pagemtimes(rx_t2, ry_t2), rz_t2) + 240*pagemtimes(rx_t2,rz_t4) + 5*ry_t6 - 120*pagemtimes(ry_t4, rz_t2) + 240*pagemtimes(ry_t2,rz_t4) - 64*rz_t6;
a_tJ6x = pagemtimes(pagemtimes((-7/16 * J6 * mu * R^6 ./ r_texp15), argJ6_t), rx_t);
a_tJ6y = pagemtimes(pagemtimes((-7/16 * J6 * mu * R^6 ./ r_texp15), argJ6_t), ry_t);
a_tJ6z = pagemtimes((-7/16 * J6 * mu * R^6 ./ r_texp15),...
    pagemtimes(35*rx_t6 + 105*pagemtimes(rx_t4, ry_t2)...
    - 210*pagemtimes(rx_t4, rz_t2) + 105*pagemtimes(rx_t2,ry_t4)...
    - 420*pagemtimes(pagemtimes(rx_t2, ry_t2), rz_t2)...
    + 168*pagemtimes(rx_t2, rz_t4) + 35*ry_t6 ...
    - 210*pagemtimes(ry_t4, rz_t2) + 168*pagemtimes(ry_t2, rz_t4) - 16*rz_t6, rz_t));
a_tJ6 = [a_tJ6x; a_tJ6y; a_tJ6z];


a_t = a_t + burn_t + a_tJ2 + a_tJ3 + a_tJ4 + a_tJ5 + a_tJ6;

if n_perSat(1) == 9
    out(1:3, :) = reshape(v_t, 3, n_bank);
    out(4:6, :) = reshape(a_t, 3, n_bank);
    a_tUnmodeled = -beta(1).*burn_t;
    out(7:9, :) = reshape(a_tUnmodeled, 3, n_bank);
else
    out(1:3, :) = reshape(v_t, 3, n_bank);
    out(4:6, :) = reshape(a_t, 3, n_bank);
end




%% Chaser 
rx_c = ic(n_perSat(1) + 1, :, :);
ry_c = ic(n_perSat(1) + 2, :, :);
rz_c = ic(n_perSat(1) + 3, :, :);


if n_perSat(2) == 9
    % accelerometer bias in body frame
    % bias_xbody = ic(n_perSat(1) + 7, :, :);
    % bias_ybody = ic(n_perSat(1) + 8, :, :);
    % bias_zbody = ic(n_perSat(1) + 9, :, :);
    bias_bodyi = ic(n_perSat(1) + (7:9), :, :);
end


r_c = ic(n_perSat(1) + (1:3), :, :);
v_c = ic(n_perSat(1) + (4:6), :, :);

r_cNorm =  pagenorm(r_c);
r_cexp5 = r_cNorm.^5;
r_cexp7 = r_cNorm.^7;
r_cexp13 = r_cNorm.^13;
r_cexp15 = r_cNorm.^15;


rx_c2 = rx_c.^2;
rx_c4 = rx_c.^4;
rx_c6 = rx_c.^6;


ry_c2 = ry_c.^2;
ry_c4 = ry_c.^4;
ry_c6 = ry_c.^6;


rz_c2 = rz_c.^2;
rz_c4 = rz_c.^4;
rz_c6 = rz_c.^6;


a_c_TBP  = - pagemtimes(mu  / (r_cNorm).^(3), r_c);

a_cJ2x = pagemtimes((3/2 * J2 * mu * R^2 ./ r_cexp5), pagemtimes((5 * pagemrdivide(rz_c, r_cNorm).^2 - 1), rx_c));
a_cJ2y = pagemtimes((3/2 * J2 * mu * R^2 ./ r_cexp5), pagemtimes((5 * pagemrdivide(rz_c, r_cNorm).^2 - 1), ry_c));
a_cJ2z = pagemtimes((3/2 * J2 * mu * R^2 ./ r_cexp5), pagemtimes((5 * pagemrdivide(rz_c, r_cNorm).^2 - 3), rz_c));
a_cJ2 = [a_cJ2x; a_cJ2y; a_cJ2z];


a_cJ3x = pagemtimes((-1/2 * J3 * mu * R^3 ./ r_cexp7),  pagemtimes(pagemtimes(5*(3 - 7*pagemrdivide(rz_c, r_cNorm).^2),  rx_c), rz_c));
a_cJ3y = pagemtimes((-1/2 * J3 * mu * R^3 ./ r_cexp7),  pagemtimes(pagemtimes(5*(3 - 7*pagemrdivide(rz_c, r_cNorm).^2),  ry_c), rz_c));
a_cJ3z = pagemtimes((1/2 * J3 * mu * R^3 ./ r_cexp5),  (3 - 30*pagemrdivide(rz_c, r_cNorm).^2 + 35*pagemrdivide(rz_c, r_cNorm).^4));
a_cJ3 = [a_cJ3x; a_cJ3y; a_cJ3z];


a_cJ4x = pagemtimes((5/8 * J4 * mu * R^4 ./ r_cexp7), pagemtimes((3 - 42*pagemrdivide(rz_c, r_cNorm).^2 + 63*pagemrdivide(rz_c, r_cNorm).^4), rx_c));
a_cJ4y = pagemtimes((5/8 * J4 * mu * R^4 ./ r_cexp7), pagemtimes((3 - 42*pagemrdivide(rz_c, r_cNorm).^2 + 63*pagemrdivide(rz_c, r_cNorm).^4), ry_c));
a_cJ4z = pagemtimes((5/8 * J4 * mu * R^4 ./ r_cexp7), pagemtimes((15 - 70*pagemrdivide(rz_c, r_cNorm).^2 + 63*pagemrdivide(rz_c, r_cNorm).^4), rz_c));
a_cJ4 = [a_cJ4x; a_cJ4y; a_cJ4z];



a_cJ5x = pagemtimes((21/8 * J5 * mu * R^5 ./ r_cexp13), pagemtimes(pagemtimes((5*rx_c4 + 10*pagemtimes(rx_c2, ry_c2) - 20*pagemtimes(rx_c2, rz_c2) + 5*ry_c4 - 20*pagemtimes(ry_c2,rz_c2) + 8*rz_c4), rx_c), rz_c));
a_cJ5y = pagemtimes((21/8 * J5 * mu * R^5 ./ r_cexp13), pagemtimes(pagemtimes((5*rx_c4 + 10*pagemtimes(rx_c2, ry_c2) - 20*pagemtimes(rx_c2, rz_c2) + 5*ry_c4 - 20*pagemtimes(ry_c2,rz_c2) + 8*rz_c4), ry_c), rz_c));
a_cJ5z = pagemtimes((-3/8 * J5 * mu * R^5 ./ r_cexp13), (5*rx_c6 + 15*pagemtimes(rx_c4, ry_c2) - 90*pagemtimes(rx_c4, rz_c2) + 15*pagemtimes(rx_c2,ry_c4) - 180*pagemtimes(pagemtimes(rx_c2, ry_c2), rz_c2) + 120*pagemtimes(rx_c2, rz_c4) + 5*ry_c6 - 90*pagemtimes(ry_c4, rz_c2) + 120*pagemtimes(ry_c2, rz_c4) - 16*rz_c6));
a_cJ5 = [a_cJ5x; a_cJ5y; a_cJ5z];


argJ6_c = 5*rx_c6 + 15*pagemtimes(rx_c4,ry_c2) - 120*pagemtimes(rx_c4,rz_c2) + 15*pagemtimes(rx_c2, ry_c4) - 240*pagemtimes(pagemtimes(rx_c2, ry_c2), rz_c2) + 240*pagemtimes(rx_c2,rz_c4) + 5*ry_c6 - 120*pagemtimes(ry_c4, rz_c2) + 240*pagemtimes(ry_c2,rz_c4) - 64*rz_c6;
a_cJ6x = pagemtimes(pagemtimes((-7/16 * J6 * mu * R^6 ./ r_cexp15), argJ6_c), rx_c);
a_cJ6y = pagemtimes(pagemtimes((-7/16 * J6 * mu * R^6 ./ r_cexp15), argJ6_c), ry_c);
a_cJ6z = pagemtimes((-7/16 * J6 * mu * R^6 ./ r_cexp15),...
    pagemtimes(35*rx_c6 + 105*pagemtimes(rx_c4, ry_c2)...
    - 210*pagemtimes(rx_c4, rz_c2) + 105*pagemtimes(rx_c2,ry_c4)...
    - 420*pagemtimes(pagemtimes(rx_c2, ry_c2), rz_c2)...
    + 168*pagemtimes(rx_c2, rz_c4) + 35*ry_c6 ...
    - 210*pagemtimes(ry_c4, rz_c2) + 168*pagemtimes(ry_c2, rz_c4) - 16*rz_c6, rz_c));
a_cJ6 = [a_cJ6x; a_cJ6y; a_cJ6z];



    if n_perSat(2) == 9
        if norm(accelerometer) > gate
            % Acceleration detected -- apply and subtract out bias
            burn_c = pagemtimes(Rbody_N', (accelerometer - bias_bodyi));
        else % no acceleration detected, do not use in prop
            burn_c = zeros(3, 1, n_bank);
        end

        a_c = a_c_TBP  + burn_c + a_cJ2 + a_cJ3 + a_cJ4 + a_cJ5 + a_cJ6;
        % out(n_perSat(1)+1:n_perSat(1)+n_perSat(2)) = [vx_c; vy_c; vz_c; a_c(1); a_c(2); a_c(3); -beta(2)*[bias_xbody; bias_ybody; bias_zbody]];
        
        out(n_perSat(1)+ (1:3), :) = reshape(v_c, 3, n_bank);
        out(n_perSat(1)+ (4:6), :) = reshape(a_c, 3, n_bank);
        b_cDot =  -beta(2).*bias_bodyi;
        out(n_perSat(1)+ (7:9), :) = reshape(b_cDot, 3, n_bank);
    else
        if nargin >= 10
            bias_bodyi = bias_body;
        else
            bias_bodyi = zeros(3, 1, n_bank);
        end
        if norm(accelerometer) > gate
            burn_c = pagemtimes(Rbody_N', (accelerometer - bias_bodyi));
        else
            burn_c = zeros(3, 1);
        end
        a_c = a_c_TBP + burn_c + a_cJ2 + a_cJ3 + a_cJ4 + a_cJ5 + a_cJ6;
        % out(n_perSat(1)+1:n_perSat(1)+n_perSat(2)) = [vx_c; vy_c; vz_c; a_c(1); a_c(2); a_c(3)];
        out(n_perSat(1)+ (1:3), :) = reshape(v_c, 3, n_bank);
        out(n_perSat(1)+ (4:6), :) = reshape(a_c, 3, n_bank);
    end
    
    if N > n_perSat(1)+ n_perSat(2)
        out(end, :, :) = -beta(3).*ic(end, :, :); % Lidar bias as Markov 
    end
    out = reshape(out, N, n_bank);
end

