%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dual Burn Rendezvous Scenario IC Script:
% Sets flags for the Dual Burn Rendezvous Scenario including tuning
% parameters
%
%
% Created: 03/07/24
% Last Modified: 03/07/24
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: Set this up to run

nMC = 1;
nMCTruth = 100;% 
firstPlot = 1;

 % 'EarthRel_gravDiff_burnRendezvous'

flagFilterType = 'EKF';                            % 'GSKF' : Gaussian Second-order Kalman Filter, 'EKF': Extended Kalman Filter '1stProp2ndMeas': first order prop, 2nd order meas update

iiTruth = 1;

% Initial covariance study
flagInitCovStudy = true;                                    % Flag to run an initial covariance study -- vary Chaser init, hold al else the same
P0Scales = [0, 1e-10, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
P0str = {'zero', '1eneg10', '1eneg5', '1eneg4', '1eneg3', '1eneg2', '1eneg1', 'none'};
iiP0 = 8;
if flagInitCovStudy
    fileBase = fullfile('.', 'MAT', flagScenario, 'InitCovFactor', P0str{iiP0}, flagFilterType);                           % Sets base file path for Figure .mat file read in
else
    fileBase = fullfile('.', 'MAT', flagScenario, flagFilterType);                           % Sets base file path for Figure .mat file read in
end

if flagInitCovStudy
    n_P0 = length(P0Scales);
else
    n_P0 = 0;
end




measNum = 3;
plotCheck = false;



% Flags for Truth model dynamics
flagGravModelTruth = 'j2j3j4On';                        % Sets grav model for truth
if strcmp(flagGravModelTruth, 'j2j3j4j5j6On')
    eomTruth = @tbpJ2_J6_accelerometerGating_rpod_lidarBias;
else
    eomTruth = @tbpJ2_J4_accelerometerGating_rpod_lidarBias;
end


mu = 3.986e5*(1e3)^3; % m^3/s^2 
J2 = .0010826269;
J3 = -0.0000025323;
J4 = -0.0000016204;
J5 = -2.27753590730836e-7;
J6 = 5.40666576283813e-7;
J7 = -3.50551795713742e-7;
J8 = -2.03993125929884e-7;
RBody = 6378*1e3; % m


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




flagTruthNoisyBurn = false;

%% Flags for Filter/MC
expNum = [6, 12];                                         % Process noise (1e-expNum)
underweight = 0;                                    % Amount of underweighting -- 0 for no underweighting
flagUnderweight = 'fixed'; %                           % Falg for underweighting, if fixed, will remain what it is set to, if 'variable' will update due to the convention outlined inUnderweighting Nonlinear Measurements by Zanetti, DeMars, and Bishop, or 'switchOff' to have it switch off after a certain uncertainty
plotExclude = []; % 7, 28

flagEstTargetAccel = true;                         % adds acceleration to state when true
flagEstChaserAccelBias = true;                      % Estimates bias on chaser accelerometer when true
flagSampleInit = true;                              % Samples from initial distribution when true
flagGravModel = 'j2j3j4On';                             % Sets gravity model for filter
flagTruthMeas = false;                              % Passes filter truth measurments when true
flagBurnOn = true;                                  % Zeros our burn hist when set to false
flagAccelerometerPerfectMatch = false;               % When false, adds noise to accelermoneter measurement based on VRW
flagRunMeasUpdate = true;                           % Runs the measurment update when true
flagQ = 'fixed';% 'CorrMethodInf'  %'CovMatchTrad' %'CorrMethodInf';    % Sets the process noise handling
adaptBuffer = 15;                                   % If running adaptive process noise, sets buffer size -- helps control how quickly we can adjust -- 20 is good for linear meas 
fullEnoughSplit = 10; % ,aybe switch to nees 
fullEnoughAdapt = 15;

alphaNIS = 0.01; % significance level
Nnx = fullEnoughSplit*3; % 3 meas
% Compute Intervals:
r1x = chi2inv(alphaNIS/2, Nnx )./ Nnx;
r2x = chi2inv(1-alphaNIS/2, Nnx )./ Nnx;
% epsilon_split = r2x*1.2; % 0.04;
 epsilon_split = r2x*35; % 0.04;
epsilon_prune = 1e-12;
flagConsider = [false(1, 3), true(1, 3), false(1, 1)]; % [false(1, 2), false, false(1, 3)];                         % Sets parameters to be consider parameters when true (lumps all 3 components together (pos, vel, accel) false(6, 1); %
beta = [0, 1/3600, 0];                            % time constants for Markov Processes (target acceleration, chaser accelerometer, lidar bias)
bias = [(300*9.81e-6); (300*9.81e-6); (300*9.81e-6)]/1000; % initial bias (300*9.81e-6) LN200S
biasLidar = 0;                                      % Lidar Bias
accelerometerVRW = eye(3)*(35e-6*9.81)^2/1000;           % (m/s/srt(s))^2 - VRW^2 -- LN200S for numbers
flagLidarBias = false;                              % Includes Lidar bias in state space when true
if contains(flagScenario, 'trueJacobians')
    flagUseTrueJacobians = true;                       % Debugging flag to give EKF and adaptation true jacobians
else
    flagUseTrueJacobians = false;                       % Debugging flag to give EKF and adaptation true jacobians
end

if flagEstTargetAccel
    flagNoiseonJerk = true;
    % eom = @tbpJ2_accelerometerGating_rpod_lidarBias_jerkNoise;
    if strcmp(flagGravModel, 'j2j3j4j5j6On')
        eom = @tbpJ2_J6_accelerometerGating_rpod_lidarBias_bankProp;
    elseif strcmp(flagGravModel, 'j2j3j4On')
        eom = @tbpJ2_J4_accelerometerGating_rpod_lidarBias_bankProp;
    else
        eom = @tbpJ2_J8_accelerometerGating_rpod_lidarBias_bankProp;
    end
else
    flagNoiseonJerk = false;
    eom = @tbpJ2_accelerometerGating_rpod_lidarBias;
end

gate = 1e-7;
Earth_NoQ_DualBurnRendezvous_DifferentOrbits_j6_scratch;

gate = 1e-4;

P0_1sat = diag([ones(3, 1)*(25/3); ones(3, 1)*(.1); ones(3, 1)*1e-6].^2);
P0_1satArray = zeros(9, n_P0);
for ii_P0 = 1:n_P0
    P0_1satArray(:, ii_P0) = diag(P0_1sat)*P0Scales(ii_P0);
end
% P0 = diag([ones(3, 1)*1e3; ones(3, 1)*1; ones(3, 1)*1e-4; ones(3, 1)*1e3; ones(3, 1)*1; ones(3, 1)*1e-4].^2);

% Buffer of 10
% Pqkminus0 = eye(3)*1e-8; % May want to come back to this at some point
% qkminus0 = ones(3, 1)*1e-8; % 10^(-expNum(1)/2);% ones(3, 1)*Q0(1, 1) + sqrt((Pqkminus*randn(3, 1)).^2);
% QW = diag([1e3, 1e3, 1e3]); %diag([1e8, 1e8, 1e8]) eye(3)*1e-17*0; good for direct meas wiht direct est of q % diag([1e10, 0, 1e10]) good for nl -- play with this
% RW = 1e-2; %1e-9 1e-9 good for direct meas with direct est of sqrt q at 1e-4 % 1e-5 good for nl meas % 1e-6 good for sqrt q 1e-3 % 1e-5 good for nl
% % RW_array = [1e4 1e5, 1e4];
% RW_array = ones(1, 3)*(RW^-1);
% 
% % Buffer of 15
% Pqkminus0 = eye(3)*1e-8; % May want to come back to this at some point
% qkminus0 = ones(3, 1)*1e-8; % 10^(-expNum(1)/2);% ones(3, 1)*Q0(1, 1) + sqrt((Pqkminus*randn(3, 1)).^2);
% QW = diag([1e4, 1e4, 1e4]); %diag([1e8, 1e8, 1e8]) eye(3)*1e-17*0; good for direct meas wiht direct est of q % diag([1e10, 0, 1e10]) good for nl -- play with this
% RW = 1e-6; %1e-9 1e-9 good for direct meas with direct est of sqrt q at 1e-4 % 1e-5 good for nl meas % 1e-6 good for sqrt q 1e-3 % 1e-5 good for nl
% RW_array = [1e1, 8e6, 1e4];

% % Buffer of 30 -- zero beta
Pqkminus0 = eye(3)*1e-6; % May want to come back to this at some point
qkminus0 = ones(3, 1)*1e-10; % 10^(-expNum(1)/2);% ones(3, 1)*Q0(1, 1) + sqrt((Pqkminus*randn(3, 1)).^2);
QW = diag([1e8, 1e8, 1e8])/100; %
RW = 1e-7; %-10 is ok 
% RW_array = [1e2, 1e6, 1e4];
RW_array = ones(3, 1)*RW^-1;

% RW_array(1) = 1e6;

% 04/17/25 -- calling it a night, but for adaptBuff=15=fullEnoughAdapt,
% fullEnoughSplit=10, epsAplit=0.025, QW 1e9, RW1e-7 performs well on run 1
% -- lol good luck getting more than two runs to work 


if contains(flagQ, 'Inf')
    Lambda0 = Pqkminus0\eye(3);
    dk0 = Lambda0*qkminus0;
    % QW = QW\eye(3);
    RW = 1/RW;
end



n_elm = zeros(2, 1);
if flagEstTargetAccel
    n_elm(1) = 9;
else
    n_elm(1) = 6;
end

if flagEstChaserAccelBias
    n_elm(2) = 9;
else
    n_elm(2) = 6;
end


% RW_array = [1e5, 1e8, 1e2]; % good for nl w/ buff of 75 noise on accel
% RW_array = [1e2, 1e8, 1e-1]; % good for nl w/ buff of 30 noise on accel
% RW_array = [1e-1, 1e8, 1e-2]; % noise on jerk 75 buffer
% RW_array = [1e2, 1e8, 1e1];
% RW_array = ones(1, 3)*RW;
% RW_array = [1e1, 1e8, 1e-3]; 

% RW_array = [1e1, 5e10, 1e-3]
RW = eye(adaptBuffer*n_elm(1)^2 + n_elm(1)*(n_elm(1) + 1)/2)*RW;

if n_elm(1) == 6
    firstIndex = [1, 7, 12; 16, 19, 21];
    for i = 1:3
        for k = 1:2
            RW(firstIndex(k, i), firstIndex(k, i)) = RW_array(i);
        end
    end

    xidex = 21 + [1, 22]';
    yidex = 21 + [8, 29]';
    zidex = 21 + [15, 36]';
    xyzIndex = [xidex, yidex, zidex];
    for i = 1:adaptBuffer
        for j = 1:3 % column
            for k = 1:2 % row
                RW(xyzIndex(k, j), xyzIndex(k, j)) = RW_array(j);
            end
        end
        xyzIndex = xyzIndex + 81;
    end
elseif n_elm(1) == 9
    firstIndex = [1, 10, 18;... % x, y, z
        25, 31, 36;... % x, y, z
        40, 43, 45];  % x, y, z
    for j = 1:3
        for k = 1:3
            RW(firstIndex(k, j), firstIndex(k, j)) = RW_array(j);
        end
    end

    xidex = 45 + [1, 31, 61]';
    yidex = 45 + [11, 41, 71]';
    zidex = 45 + [21, 51, 81]';
    xyzIndex = [xidex, yidex, zidex];
    for i = 1:adaptBuffer
        for j = 1:3 % column
            for k = 1:3 % row
                RW(xyzIndex(k, j), xyzIndex(k, j)) = RW_array(j);
            end
        end
        xyzIndex = xyzIndex + 81;
    end
end

if flagInitCovStudy
    P0_1sat = P0_1sat*P0Scales(iiP0);
end