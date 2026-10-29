function [out] = CRTBP_targetChaser_bankProp(t, ic, gravParam, beta, accelerometer, n_perSat, gate, Rbody_R, RN_rot, bias_body, Ccoeff, Scoeff )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Sets up the Equations of Motion (EOM) to integrate the circular
% restricted three body problem for a given three body system -- taken wrt
% rotatimg frame 
%
% Inputs:
%       t             nondimensional time [-]
%       ic            initial conditions:
%                       nondimenasional target position vector [-], 
%                       nondimenaisonal target velocity vector [-],
%                       nondimensional target unmodeled acceleration [-],
%                       nondimenasional chaser position vector [-], 
%                       nondimenaisonal chaser velocity vector [-],
%                       nondimensional chaser accelerometer bias [-],
%       gravParam       gravitational parameter mass ratio [-],
%                           characteristic length [km], mean motion [1/s],
%                           characteristic time [s]
%       beta           nondimensional 
%       accelerometer  nondimesaional chaser acceleration in the body frame
%                       [-]
%       n_perSat       number of states for each body 2x1 vec
%       gate           accelerometer gate thershold 
%       Rbody_N        DCM from inertial to body frame of each s/c
%       RN_rot         n/a for 3BP
%       bias_body      accelerometer bias -- if not included in ic 
%       Ccoeff         n/a for 3BP
%       Scoeff         n/a for 3BP
%
% Outputs:
%       out           current state:
%                       postion vector [km], velocity vector [km / s]
%
% Author: Rachel Mamich
% Created: 03/12/19
% Modified: 08/01/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N, n_bank] = size(ic);
ic = reshape(ic, N, 1, n_bank);
out = zeros(N, 1, n_bank);
mu = gravParam(1); % t, L, n, tCharacteristic
L = gravParam(2);
n = gravParam(3);
tCharacteristic = gravParam(4);


%% Target
rx_t = ic(1, :, :);
ry_t = ic(2, :, :);
rz_t = ic(3, :, :);
vx_t = ic(4, :, :);
vy_t = ic(5, :, :);
vz_t = ic(6, :, :);



r_t = ic(1:3, :, :);
v_t = ic(4:6, :, :);

if n_perSat(1) == 9
    % accelerometer bias in body frame
    burn_t = ic(7:9, :, :);
else
    burn_t = zeros(3, 1, n_bank);
end


r1Vec_t = [rx_t + mu; ry_t; rz_t];
r2Vec_t = [rx_t-1+mu; ry_t; rz_t];


r1_t = pagenorm(r1Vec_t);
r1_texpneg3 = r1_t.^(-3);

r2_t = pagenorm(r2Vec_t);
r2_texpneg3 = r2_t.^(-3);


rx_tmu = rx_t+mu;
rx_t1mu = rx_t-1+mu;

onemu = 1-mu;

% r1_t = sqrt((rx_t + mu)^2 + ry_t^2 + rz_t^2);
% r2_t = sqrt((rx_t - 1 + mu)^2 + ry_t^2 + rz_t^2);
ax_t = 2*vy_t  + rx_t -onemu*pagemtimes(rx_tmu,r1_texpneg3) - mu*pagemtimes(rx_t1mu,r2_texpneg3);
ay_t = -2*vx_t + ry_t -onemu*pagemtimes(ry_t,r1_texpneg3) -mu*pagemtimes(ry_t,r2_texpneg3);
az_t = -onemu*pagemtimes(rz_t,r1_texpneg3) -mu*pagemtimes(rz_t,r2_texpneg3);

a_t = [ax_t; ay_t; az_t] + burn_t;

if n_perSat(1) == 9
    out(1:9, :, :) = [vx_t; vy_t; vz_t; a_t; -beta(1)*burn_t];
else
    out(1:6, :, :) = [vx_t; vy_t; vz_t; a_t];
end


% r1=sqrt((rx_t+mu)^2+ry_t^2+rz_t^2);
% r2=sqrt((rx_t-1+mu)^2+ry_t^2+rz_t^2);
% c1=(1-mu)/r1^3;
% c2=mu/r2^3;
% f=[vx_t; vy_t; vz_t; 2*vy_t+rx_t-c1*(rx_t+mu)-c2*(rx_t-1+mu)+ burn_t(1); -2*vx_t+ry_t-c1*ry_t-c2*ry_t+ burn_t(2); -c1*rz_t-c2*rz_t+ burn_t(3)];


%% Chaser
rx_c = ic(n_perSat(1) + 1, :, :);
ry_c = ic(n_perSat(1) + 2, :, :);
rz_c = ic(n_perSat(1) + 3, :, :);
vx_c = ic(n_perSat(1) + 4, :, :);
vy_c = ic(n_perSat(1) + 5, :, :);
vz_c = ic(n_perSat(1) + 6, :, :);

if n_perSat(2) == 9
    % accelerometer bias in body frame
    bias_xbody = ic(n_perSat(1) + 7, :, :);
    bias_ybody = ic(n_perSat(1) + 8, :, :);
    bias_zbody = ic(n_perSat(1) + 9, :, :);
    bias_bodyi = [bias_xbody; bias_ybody; bias_zbody];
end


r1Vec_c = [rx_c + mu; ry_c; rz_c];
r2Vec_c = [rx_c-1+mu; ry_c; rz_c];


r1_c = pagenorm(r1Vec_c);
r1_cexpneg3 = r1_c.^(-3);

r2_c = pagenorm(r2Vec_c);
r2_cexpneg3 = r2_c.^(-3);


rx_cmu = rx_c+mu;
rx_c1mu = rx_c-1+mu;

onemu = 1-mu;


ax_c = 2*vy_c  + rx_c -onemu*pagemtimes(rx_cmu,r1_cexpneg3) - mu*pagemtimes(rx_c1mu,r2_cexpneg3);
ay_c = -2*vx_c + ry_c -onemu*pagemtimes(ry_c,r1_cexpneg3) -mu*pagemtimes(ry_c,r2_cexpneg3);
az_c = -onemu*pagemtimes(rz_c,r1_cexpneg3) -mu*pagemtimes(rz_c,r2_cexpneg3);


a_c = [ax_c; ay_c; az_c];



if n_perSat(2) == 9
    if norm(accelerometer) > gate
        burn_c = pagemtimes(Rbody_R', (accelerometer - bias_bodyi));
        % [out_acc] = nonDimRotational_dimInertial_conversion(burn_c, t, L, n, tCharacteristic, 'accelN2R_posVelR');
        % burn_c = out_acc(7:9); 
    else
        burn_c = zeros(3, 1, n_bank);
    end
    a_c = a_c + burn_c;
    out(n_perSat(1)+1:n_perSat(1)+n_perSat(2), :, :) = [vx_c; vy_c; vz_c; a_c; -beta(2)*[bias_xbody; bias_ybody; bias_zbody]];
else % TODO: Edit this else statement to work properly -- might never use it, so not doing it now -- like I haven't touched antthing in this else 
    % if nargin >= 10
    %     bias_bodyi = bias_body;
    % else
    %     bias_bodyi = zeros(3, 1);
    % end
    if norm(accelerometer) > gate
        burn_c = pagemtimes(Rbody_R',(accelerometer - bias_bodyi));
    else
        burn_c = zeros(3, 1, n_bank);
    end
    a_c = a_c + burn_c;
    out(n_perSat(1)+1:n_perSat(1)+n_perSat(2), :, :) = [vx_c; vy_c; vz_c; a_c];
end

if N > n_perSat(1)+ n_perSat(2) % haven't touched this either
    out(end, :, :) = -beta(3)*ic(end, :, :); % Lidar bias as Markov
end


out = reshape(out, N, n_bank);
end

