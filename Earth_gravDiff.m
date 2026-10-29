%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple Simulation Truth and Measurment Gen About Earth
% No burns, but noisy relative target state 
% 
%
% Created: 04/03/23
% Last Modified: 04/03/23
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%house keeping
if ~exist('flagSetFlags', 'var')
    clear
    close all
    flagScenario = 'Earth';
    fileBase = fullfile('.', 'MAT', flagScenario);
    measNum = 3;
    flagGravModelTruth = 'j2On';                      % Sets label for saving data with the truth data model
    plotCheck = true;
    beta = 0*[0, 1/3600, 0];                            % time constants for Markov Processes
    bias = [(300*9.81e-6); (300*9.81e-6); (300*9.81e-6)]; % initial bias (300*9.81e-6) LN200S
    biasLidar = 0;
    gate = 1e-7;
    eom = @tbpJ2_accelerometerGating_rpod_lidarBias;
    flagTruthNoisyBurn = false;
    iiTruth = 1;
    expNum = [8, 0];  
    eomTruth = @tbpJ2_accelerometerGating_rpod_lidarBias;
    flagNoiseonJerk = false;
end
rng(1263 + iiTruth+1)
fprintf('Running Truth %d...\n', iiTruth)
%% Initialization

% Initial Chief Spacecraft Keplarian Elements
ac = (6378+500)*1e3; % m
ec = 0.0;
incc = 95*pi/180;  % rad
Omegac = 0; % rad
omegac = 0; % rad
Mc = 0; % rad
uc = Mc + omegac;


% Initial Deputy Spacecraft Keplarian Elements

ad = ac;
ed = ec;
incd = incc;
Omegad = Omegac;
omegad = omegac;
Md = 0.000005; % rad


% Convert to inertial elements
mu = 3.986e5*(1e3)^3; % m^3/s^2 
J2 = .0010826269;
J3 = -0.0000025323;
J4 = -0.0000016204;
J5 = -2.27753590730836e-7;
J6 = 5.40666576283813e-7;
J7 = -3.50551795713742e-7;
J8 = -2.03993125929884e-7;
RBody = 6378*1e3; % m

P = 2*pi*sqrt(ac^3/mu);
dt = 1; % 2.5*3600/30;
time = 0:dt:5500; % P*1.5;
N = length(time);

truthState_N = zeros(19, N);

Phi_true = zeros(19, 19, N);


inputFlag = "kep";
outputFlag = "cartesian";

[ xc ] = orbitalElmConverter( inputFlag, [ac; ec; incc*180/pi; Omegac*180/pi; omegac*180/pi; Mc*180/pi], 0, outputFlag, mu );
[ xd ] = orbitalElmConverter( inputFlag, [ad; ed; incd*180/pi; Omegad*180/pi; omegad*180/pi; Md*180/pi], 0, outputFlag, mu );



burnStart = {round(2*P/dt); round(3*P/dt)}; % time step 803
burnEnd = {burnStart{1} + 60/dt; burnStart{2} + 60/dt}; % time step 806
burn = {0;0};
burnMagSigma = 0.15;
burnDirSigma = 0.5; % deg

n_burn = [length(burnStart{1}), length(burnStart{2})];


if strcmp(flagGravModelTruth, '8x8On')
    gravParamTruth = [mu, RBody, J2, J3, J4, J5, J6, J7, J8, 1];
elseif strcmp(flagGravModelTruth, 'j2On')
    gravParamTruth = [mu, RBody, J2, 0, 0, 0, 0, 0, 0, 0];
elseif strcmp(flagGravModelTruth, 'j2j3On')
    gravParamTruth = [mu, RBody, J2, J3, 0, 0, 0, 0, 0, 0];
elseif strcmp(flagGravModelTruth, 'j2j3j4On')
    gravParamTruth = [mu, RBody, J2, J3, J4, 0, 0, 0, 0, 0];
elseif strcmp(flagGravModelTruth, 'j2j3j4j5On')
    gravParamTruth = [mu, RBody, J2, J3, J4, J5, 0, 0, 0, 0];
elseif strcmp(flagGravModelTruth, 'j2j3j4j5j6On')
    gravParamTruth = [mu, RBody, J2, J3, J4, J5, J6, 0, 0, 0];
elseif strcmp(flagGravModelTruth, 'j2j3j4j5j6j7On')
    gravParamTruth = [mu, RBody, J2, J3, J4, J5, J6, J7, 0, 0];
elseif strcmp(flagGravModelTruth, 'j2j3j4j5j6j7j8On')
    gravParamTruth = [mu, RBody, J2, J3, J4, J5, J6, J7, J8, 0];
end

burnHistCommand = zeros(6, N);
burnHistActual = zeros(6, N);
accelerometer = zeros(3, N);




calcBurn = false;

initialState = [xc; zeros(3, 1); xd; bias; biasLidar];
truthState_N(:, 1) = initialState';
iminus1 = 1;
if dt > 1
    dtProp = 1;
else
    dtProp = dt;
end

% DCM Info for Body Frame
omegaEarth = 360 /(24*3600); % deg/s
lam0 = 0;
RBody1_N = NaN(3, 3, N);
RBody2_N = NaN(3, 3, N);
RCF1_Body1_N = NaN(3, 3, N);
RCF2_Body1_N = NaN(3, 3, N);
RCF1_Body2_N = NaN(3, 3, N);
RCF2_Body2_N = NaN(3, 3, N);

RN_rot = NaN(3, 3, N);
RN_rot(:, :, 1) = eye(3);

pointNoise = eye(3)*(1e-8)*0;
if norm(pointNoise) ~= 0
    pointNoise = chol(pointNoise);
end
RBody1_N(:, :, 1) = cameraPoint(truthState_N(1:6, 1), truthState_N(10:15, 1), pointNoise);
RBody2_N(:, :, 1) = cameraPoint(truthState_N(10:15, 1), truthState_N(1:6, 1), pointNoise);

burnCount = [1, 0];

n_exp = length(expNum);
Q0 = zeros(3*n_exp);
for i = 1:n_exp
    if expNum(i) == 0
        Q0(((i*3-2):i*3), ((i*3-2):i*3)) = zeros(3);
    else
        Q0(((i*3-2):i*3), ((i*3-2):i*3)) = eye(3)*10^(-expNum(i)); % 1e-6 is good for no accel -- 1e-12 is good for accel on
    end
end

if norm(Q0(1:3, 1:3)) > 0
    Qchol_t = chol(Q0(1:3, 1:3), 'lower');
else 
    Qchol_t = zeros(3);
end
if norm(Q0(4:6, 4:6))
    Qchol_c = chol(Q0(4:6, 4:6), 'lower');
else
    Qchol_c = zeros(3);
end

for i = 2:N
    % DCMs
    RN_rot(:, :, i) = rotz(lam0 + omegaEarth*time(i));
    RBody1_N(:, :, i) = cameraPoint(truthState_N(1:6, iminus1), truthState_N(10:15, iminus1), pointNoise);
    RBody2_N(:, :, i) = cameraPoint(truthState_N(10:15, iminus1), truthState_N(1:6, iminus1), pointNoise);
    
    timeInt = time(iminus1:i);
    burnVecCommand = zeros(6, 1);
    burnVecActual = zeros(6, 1);
    burnHistActual(:, i) = burnVecActual;
    burnHistCommand(:, i) = burnVecCommand;
    timeProp = timeInt(1):dtProp:timeInt(2);
    Nprop = length(timeProp)-1;
    initialState = truthState_N(:, iminus1);
    % if ~flagNoiseonJerk
    %     initialState(7:9) = Qchol_t*randn(3, 1);
    % % else
    % %     initialState(7:9) = Qchol_t*randn(3, 1);
    % end
    accelerometer(:, i) = burnVecActual(4:6) + initialState(16:18);% + Qchol_c*randn(3, 1);
    Rbody_N = RBody2_N(:, :, iminus1);
    Phi_true(:, :, i) = eye(19);
    for j = 1:Nprop
        if flagNoiseonJerk
            state = rk4(initialState, timeProp(j:j+1), eomTruth, gravParamTruth, beta, accelerometer(:, i), gate, Rbody_N, 0*randn(3, 1), [9, 9], {});
        else
            state = rk4(initialState, timeProp(j:j+1), eomTruth, gravParamTruth, beta, accelerometer(:, i), gate, Rbody_N, RN_rot, [9, 9], {});
        end
        [F_i] = buildPhi(initialState, gravParamTruth, dt, beta, [9, 9], true);
        Phi_true(:, :, i) = F_i*Phi_true(:, :, i);
        initialState = state;
    end
    truthState_N(:, i) = state;
    RBody1_N(:, :, i) = cameraPoint(truthState_N(1:6, i), truthState_N(10:15, i), pointNoise);
    RBody2_N(:, :, i) = cameraPoint(truthState_N(10:15, i), truthState_N(1:6, i), pointNoise);
    
    iminus1 = i;
    if mod(i, 1000) == 0
       fprintf('%d\n', i) 
    end
end

truthState = truthState_N;
truthState(1:6, :) = truthState_N(1:6, :) - truthState_N(10:15, :);



%% Verify Traj
if plotCheck
    figure(201)
    plot3(truthState_N(1, :), truthState_N(2, :), truthState_N(3, :))
    hold on
    grid minor
    plot3(truthState_N(10, :), truthState_N(11, :), truthState_N(12, :))
end

%% Generate Measuremnts
% From D'Amico's paper 
focalLength = 7.6e-3; % m

% From ASC Flash LiDAR -- IEEE paper --amzajerdian2013lidar
fov_rad = [5; 5]*pi/180; % rad
numPixel = [128; 128]; % -

R_opNav = fov_rad./numPixel;

sigma_r = 0.08; % [-] 
sigma_rDot = 0.001; % [-]

% From Jena Optronik RVS 3000-X Lidar -- GNSki 2024
% R_opNav = ones(2, 1)*(0.1*pi/180); % rad^2 -- 

if measNum == 6
    R = diag([sigma_r, sigma_rDot, R_opNav', R_opNav'].^2);
    R_rangeRangeRate_chol = chol(R(1:2, 1:2), 'lower');
    R_rangeRangeRate = R(1:2, 1:2);
    R_opNav_chol = chol(R(3:4, 3:4), 'lower');
elseif measNum == 4
    R = diag([sigma_r, sigma_rDot, R_opNav'].^2);
    R_rangeRangeRate_chol = chol(R(1:2, 1:2), 'lower');
    R_rangeRangeRate = R(1:2, 1:2);
    R_opNav_chol = chol(R(3:4, 3:4), 'lower');
elseif measNum == 3
    if contains(flagScenario, 'direct')
        R = diag([sigma_r, R_opNav'].^2);

        sigmarho = sqrt(R(1, 1));
        sigmaA = sqrt(R(2, 2));
        sigmaB = sqrt(R(3, 3));
        sigmaz = sqrt(sigmarho^2*cos(sigmaA)^2 / (1 + cos(sigmaA)^2*tan(sigmaB)^2));

        denomsigmaZ = (1 + cos(sigmaA)^2*tan(sigmaB)^2);
        dzdrho = sigmarho*cos(sigmaA)^2/denomsigmaZ/sigmaz;
        dzdA = sigmarho^2/sigmaz*cos(sigmaA)*sin(sigmaA)*(cos(sigmaA)^2/denomsigmaZ^2*tan(sigmaB)^2 - 1/denomsigmaZ);
        dzdB = -sigmarho^2/denomsigmaZ^2*tan(sigmaB)*cos(sigmaA)^4*sec(sigmaB)^2/sigmaz;


        dxdrho = tan(sigmaA)*dzdrho;
        dxdA = dzdA*tan(sigmaA) + sigmaz*sec(sigmaA)^2;
        dxdB = dzdB*tan(sigmaA);

        dydrho = tan(sigmaB)*dzdrho;
        dydA = tan(sigmaB)*dzdA;
        dydB = tan(sigmaB)*dzdB + sec(sigmaB)^2*sigmaz;
        
        T = [dxdrho, dxdA, dxdB; dydrho, dydA, dydB; dzdrho, dzdA, dzdB];

        R = T*R*T'; 
        R_direct_chol = chol(R, 'lower');
    else
        % range, relative bearing
        R = diag([sigma_r, R_opNav'].^2);
        R_rangeRangeRate_chol = sqrt(R(1, 1));
        R_rangeRangeRate = R(1, 1);
        R_opNav_chol = chol(R(2:3, 2:3), 'lower');
    end
elseif measNum == 5
    % range, relative bearing 
    R = diag([sigma_r, R_opNav', R_opNav'].^2);
    R_rangeRangeRate_chol = sqrt(R(1, 1));
    R_rangeRangeRate = R(1, 1);
    R_opNav_chol = chol(R(2:3, 2:3), 'lower');
end

meas = NaN(measNum, N);
meas_truth = NaN(measNum, N);


% Assume that ACAF and ACI are coaligned at the start time
% Assume that ACAF rotates as a constant rate of omegaErros about the z
% axis
% Assume both s/c are keeping boresight approximately towards the other s/c
% groundFeature(x, RCF_N, RN_rot, L)
N_feature = 1;
% [L_set] = genFeatures(REarth, N_feature);
L_set = [0;0;0];


% Assumes boresight is in the z direction of the camera frame
RCF1_Body1 = rotx(-90); % looking at body 2
RCF2_Body1 = roty(90); % looking toward the central body (ish)
RCF1_Body2 = rotx(-90); % looking at body 1
RCF2_Body2 = roty(90); % looking toward the central body (ish)

H = NaN(measNum, 19, N);

featureIndex = NaN(1, N);

for i = 1:N
    RCF1_Body1_N(:, :, i) = RCF1_Body1*RBody1_N(:, :, i);
    RCF2_Body1_N(:, :, i) = RCF2_Body1*RBody1_N(:, :, i);
    RCF1_Body2_N(:, :, i) = RCF1_Body2*RBody2_N(:, :, i);
    RCF2_Body2_N(:, :, i) = RCF2_Body2*RBody2_N(:, :, i); 
    
    if measNum == 4
        [meas_i, H_i] = rngRngRate_opTarget_meas_H(truthState(:, i), RCF1_Body2_N(:, :, i), [9, 9]);
        meas_truth(:, i) = meas_i;
        meas(1:2, i) = meas_i(1:2) + R_rangeRangeRate_chol*randn(2, 1);
        meas(3:4, i) = meas_i(3:4) + R_opNav_chol*randn(2, 1);
        H(:, :, i) = H_i;
    elseif measNum == 5
        featureIndex(i) = 1;
        [meas_i, H_i] = relNavRngRngRate_opTarget_feature_meas_H(truthState(:, i), RCF1_Body2_N(:, :, i), RCF2_Body2_N(:, :, i), RN_rot(:, :, i), L_set(:, 1), [9, 9]);
        meas_truth(:, i) = meas_i;
        meas(1, i) = meas_i(1) + R_rangeRangeRate_chol*randn(1, 1);
        meas(2:3, i) = meas_i(2:3) + R_opNav_chol*randn(2, 1);
        meas(4:5, i) = meas_i(4:5) + R_opNav_chol*randn(2, 1);
        H(:, :, i) = H_i;
    elseif measNum == 3
        if contains(flagScenario, 'Rel')
            if contains(flagScenario, 'direct')
                [meas_i, H_i] = relNavDirectInertial_meas_H(truthState(:, i));
                meas_truth(:, i) = meas_i;
                meas(:, i) = meas_i + R_direct_chol*randn(3, 1);
                H(:, :, i) = H_i;
            else
                [meas_i, H_i] = relNavRng_opTarget_meas_H(truthState(:, i), RCF1_Body2_N(:, :, i), [9, 9]);
                meas_truth(:, i) = meas_i;
                meas(1, i) = meas_i(1) + R_rangeRangeRate_chol*randn(1, 1);
                meas(2:3, i) = meas_i(2:3) + R_opNav_chol*randn(2, 1);
                H(:, :, i) = H_i;
            end
        else 
            [meas_i, H_i] = rng_opTarget_meas_H(truthState_N(:, i), RCF1_Body2_N(:, :, i), [9, 9]);
            meas_truth(:, i) = meas_i;
            meas(1, i) = meas_i(1) + R_rangeRangeRate_chol*randn(1, 1);
            meas(2:3, i) = meas_i(2:3) + R_opNav_chol*randn(2, 1);
            H(:, :, i) = H_i;
        end
    end
    

end

% meas = meas_truth;

%% Save Off Info
if ~exist(fileBase, 'dir')
    mkdir(fileBase)
end

if ~contains(flagScenario, 'Rel')
    truthState = truthState_N;
end

save(fullfile(fileBase, strcat('truthState.mat')), 'truthState')
save(fullfile(fileBase, 'meas_.mat'), 'meas')
save(fullfile(fileBase, 'truemeas.mat'), 'meas_truth')
save(fullfile(fileBase, 'groundFeatIndex.mat'), 'featureIndex')
save(fullfile(fileBase, 'acaf2aci.mat'), 'RN_rot')
save(fullfile(fileBase, 'aci2body1CF1.mat'), 'RCF1_Body1_N')
save(fullfile(fileBase, 'aci2body1CF2.mat'), 'RCF2_Body1_N')
save(fullfile(fileBase, 'aci2body2CF1.mat'), 'RCF1_Body2_N')
save(fullfile(fileBase, 'aci2body2CF2.mat'), 'RCF2_Body2_N')
save(fullfile(fileBase, 'body12N.mat'), 'RBody1_N')
save(fullfile(fileBase, 'body22N.mat'), 'RBody2_N')
save(fullfile(fileBase, 'timeMat.mat'), 'time')
save(fullfile(fileBase, 'burnHistActual.mat'), 'burnHistActual')
save(fullfile(fileBase, 'burnHistCommand.mat'), 'burnHistCommand')
save(fullfile(fileBase, 'accelerometerTruth.mat'), 'accelerometer')
save(fullfile(fileBase, 'Rmat.mat'), 'R')
save(fullfile(fileBase, 'Hmat.mat'), 'H')
save(fullfile(fileBase, 'featureSet.mat'), 'L_set');
save(fullfile(fileBase, 'trueSTM.mat'), 'Phi_true');
save(fullfile(fileBase, strcat('truthStateMC', num2str(iiTruth), '.mat')), 'truthState')
save(fullfile(fileBase, strcat('truemeasMC', num2str(iiTruth),'.mat')), 'meas_truth')

%% Functoins Called Up Top
function [RBody_N] = cameraPoint(xc, xd, pointNoise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Points body y at the other spacecraft with the body x
% approximately towards the central body
%
% Assumptions: Assumes that boresight is in the y direction on the body
% frame
%
% Inputs:
%   xc          state of the chief satellite
%   xd          state of the deputy satellite (on board)
%   pointNoise  noise associated with the vehicle pointing
%
% Output:
%   RCF_N       DCM that defines the rotation from the inertial frame to
%                   the camera frame
%
% Author: Rachel Mamich
% Creater: 09/12/2022
% Modified: 09/12/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rc = xc(1:3);
rd = xd(1:3);
yaxis = rc - rd;
yaxis = yaxis/norm(yaxis);

rdUnit_N = rd/norm(rd);
zaxis = cross(rdUnit_N, yaxis);
zaxis = zaxis/norm(zaxis);
xaxis = cross(yaxis, zaxis);

if norm(pointNoise) == 0
    noise = eye(3);
else
    noise = pointNoise*randn(3, 1);
    noise = PRV2C(noise);
end
RBody_N = [xaxis, yaxis, zaxis]'*noise;

end


function [RCF_N] = cameraPoint_atFeature(x, pointNoise, L, RN_rot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Points body y at the other spacecraft with the body x
% approximately towards the central body
%
% Assumptions: Assumes that boresight is in the y direction on the body
% frame
%
% Inputs:
%   x           state of the satellite
%   pointNoise  noise associated with the vehicle pointing
%
% Output:
%   RCF_N       DCM that defines the rotation from the inertial frame to
%                   the camera frame
%
% Author: Rachel Mamich
% Creater: 09/20/2022
% Modified: 09/20/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = (RN_rot*L - x(1:3));
v = x(4:6);
zaxis = r;
zaxis = zaxis/norm(zaxis);

xaxis = cross(r', v')';
xaxis = xaxis/norm(xaxis);

yaxis = cross(zaxis', xaxis')';
yaxis = yaxis/norm(yaxis);

if norm(pointNoise) == 0
    noise = eye(3);
else
    noise = pointNoise*randn(3, 1);
    noise = PRV2C(noise);
end

RCF_N = [xaxis, yaxis, zaxis]'*noise;

end

function [C] = PRV2C(q)
%     """
%     PRV2C
%
%     	C = PRV2C(Q) returns the direction cosine
%     	matrix in terms of the 3x1 principal rotation vector
%     	Q.
%     """
% Modified from Basilisk

q0 = norm(q);
if q0 == 0.0
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    fprintf('Quaternion not Valid\n')
    C = eye(3);
else
    q1 = q(1) / q0;
    q2 = q(2) / q0;
    q3 = q(3) / q0;
    cp = cos(q0);
    sp = sin(q0);
    d1 = 1 - cp;
    C = zeros(3);
    C(1, 1) = q1 * q1 * d1 + cp;
    C(1, 2) = q1 * q2 * d1 + q3 * sp;
    C(1, 3) = q1 * q3 * d1 - q2 * sp;
    C(2, 1) = q2 * q1 * d1 - q3 * sp;
    C(2, 2) = q2 * q2 * d1 + cp;
    C(2, 3) = q2 * q3 * d1 + q1 * sp;
    C(3, 1) = q3 * q1 * d1 + q2 * sp;
    C(3, 2) = q3 * q2 * d1 - q1 * sp;
    C(3, 3) = q3 * q3 * d1 + cp;
end
end


function [L_set] = genFeatures(R, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Generates N random points on a sphere of radius R
%
%
% Inputs:
%   R       radius of sphere [m]
%   N       number of points
%
% Output:
%   L_set   points on the surface of the sphere [m]
%
% Author: Rachel Mamich
% Creater: 09/13/2022
% Modified: 09/13/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angles = unifrnd(0, 2*pi, [2, N]);
[~, index] = sort(angles(2, :));
angles(1, :) = angles(1, index);
angles(2, :) = angles(2, index);
L_set = NaN(3, N);
for i = 1:N
    L_set(:, i) = R*[cos(angles(1, i))*cos(angles(2, i)); cos(angles(1, i))*sin(angles(2, i)); sin(angles(1, i))];
end

end

function [meas] = inFOV(RN_rot, L, x, RCF_N, fov)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Generates N random points on a sphere of radius R
%
%
% Inputs:
%   RN_rot  DCM from cb centered cb fixed frame to inertial
%   L       feature location in cd centered cb fixed frame
%   x       state
%   RCF_N   DCM from inertial to camera frame
%   fov     field of view of the camera
%
% Output:
%   flag    if true, feature is in fov
%
% Author: Rachel Mamich
% Creater: 09/13/2022
% Modified: 09/13/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rd_N = x(10:12); % m
L_N = RN_rot*L; % m

meas = NaN(2, 1);
if wrapToPi(acos(dot(rd_N, L_N) / (norm(rd_N)*norm(L_N)))) < pi/2
    [angles, ~] = pinHoleOpticalMeas_groundFeature(x, RCF_N, RN_rot, L);
    if angles(1) < fov(1)
        if angles(2) < fov(2)
            meas = angles;
        end
    end
end

end


function [Pi]= normalizationCoeff_Pi(L, m)
if m == 0
    k = 1.0;
else
    k = 2.0;
end
Pi = sqrt((factorial(L + m)) / (factorial((L - m))*k*(2*L + 1)));
end