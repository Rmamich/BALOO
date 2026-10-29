%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BALOO 
% Venting Scenario IC Script:
% Sets flags for the Venting Scenario including tuning parameters
%
%
% Created: 08/06/24
% Last Modified: 08/06/24
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Flags to Set
% 'NRHO_Venting_RelNav_nonDim'
flagFilterType = 'EKF';                            % 'GSKF' : Gaussian Second-order Kalman Filter, 'EKF': Extended Kalman Filter '1stProp2ndMeas': first order prop, 2nd order meas update
nMC = 1;
nMCTruth = 100;
firstPlot = 1; % 70-61 is a bit weird -- warning flag went off for plotting  81-90 is fully horrible 
plotExclude = []; % 7, 28

% NRHOVenting_Scratch_keepitNonDim;
% From Lorin's AAS GNC Paper 2024  
mu = 0.012150584270; % mass ratio from paper
L = 384400; % km -- characteristic length 
L_m = L*1000;
tCharacteristic = 3.751902619517e5; % s -- characteristic time from paper  

PNondim = 1.511143593137; % period non dimensional from paper

tDim = PNondim*tCharacteristic; % dimensional period -- convert as nondimensional period * characteristic time 
P = tDim; %
n = 2*pi/P; % 1/s 


gravParamTruth = [mu, L, n, tCharacteristic]; 

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
flagGravModelTruth = 'CRTBP';                        % Sets grav model for truth

eom = @CRTBP_targetChaser_bankProp;        % Can swap out for different dynamics, but will enforce that EKF prop functions go with EKF and GSKF goes with GSKF
eomTruth = @CRTBP_targetChaser;
flagTruthNoisyBurn = false;

%% Flags for Filter/MC
expNum = [12, 0];  % Fixed EKF runs well with [2, 0]   % Process noise (1e-expNum)
% expNum = [6, 0];
underweight = 0.0;                                    % Amount of underweighting -- 0 for no underweighting
flagUnderweight = 'fixed'; %                           % Falg for underweighting, if fixed, will remain what it is set to, if 'variable' will update due to the convention outlined inUnderweighting Nonlinear Measurements by Zanetti, DeMars, and Bishop, or 'switchOff' to have it switch off after a certain uncertainty


flagEstTargetAccel = true;                         % adds acceleration to state when true
flagEstChaserAccelBias = true;                      % Estimates bias on chaser accelerometer when true
flagSampleInit = true;    %                           % Samples from initial distribution when true
flagGravModel = 'CRTBP';                             % Sets gravity model for filter
flagTruthMeas = false;                              % Passes filter truth measurments when true
flagBurnOn = true;                                  % Zeros our burn hist when set to false
flagAccelerometerPerfectMatch = false;               % When false, adds noise to accelermoneter measurement based on VRW
flagRunMeasUpdate = true;                           % Runs the measurment update when true
flagQ = 'CorrMethodInf'; %'fixed';%% 'CorrMethodInf'  %'CovMatchTrad' %'CorrMethodInf';    % Sets the process noise handling
adaptBuffer = 15;                                   % If running adaptive process noise, sets buffer size -- helps control how quickly we can adjust -- 20 is good for linear meas 
fullEnoughSplit = 10; % ,aybe switch to nees 
fullEnoughAdapt = 15;

alphaNIS = 0.01; % significance level
Nnx = fullEnoughSplit*3; % 3 meas
% Compute Intervals:
r1x = chi2inv(alphaNIS/2, Nnx )./ Nnx;
r2x = chi2inv(1-alphaNIS/2, Nnx )./ Nnx;
% epsilon_split = r2x*1.2; % 0.04;
 epsilon_split = 50;%r2x*60; % 0.04; 1.68;%
epsilon_prune = 1e-12;
flagConsider = [false(1, 3), true(1, 2), false(1, 1), true(1, 1)]; % false(1, 7);%  [false(1, 2), false, false(1, 3)];                         % Sets parameters to be consider parameters when true (lumps all 3 components together (pos, vel, accel) false(6, 1); %
beta = [0, 1e-5, 0]*0;                            % time constants for Markov Processes (target acceleration, chaser accelerometer, lidar bias)
bias = [(300*9.81e-6); (300*9.81e-6); (300*9.81e-6)]/100/(L_m)*tCharacteristic^2; % initial bias (300*9.81e-6) LN200S
% bias = nonDimRotational_dimInertial_conversion(bias, 0, L, n, tCharacteristic, 'accelerometer');
biasLidar = 0;                                      % Lidar Bias
gate = 1e-4/(L_m)*tCharacteristic^2;
% gate = nonDimRotational_dimInertial_conversion(gate, 0, L, n, tCharacteristic, 'accelerometer');
accelerometerVRW = eye(3)*(35e-6*9.81/100)^2 / (L_m/(tCharacteristic^(3/2)))^2*1e-6;           % (m/s/srt(s))^2 - VRW^2 -- LN200S for numbers [-]
flagLidarBias = false;                              % Includes Lidar bias in state space when true
if contains(flagScenario, 'trueJacobians')
    flagUseTrueJacobians = true;                       % Debugging flag to give EKF and adaptation true jacobians
else
    flagUseTrueJacobians = false;                       % Debugging flag to give EKF and adaptation true jacobians
end
flagNoiseonJerk = false;

% look into partials closer 

% P0_1sat = diag([ones(3, 1)*(25/3)/L_m; ones(3, 1)*(.1)/L_m*tCharacteristic; ones(3, 1)*1e-6/(L_m)*tCharacteristic^2].^2)/100;
P0_1sat = diag([ones(3, 1)*2.5e-9; ones(3, 1)*2.5e-11; ones(3, 1)*3.5e-7].^2); % from Dalton's paper on Cislunar OD 
P0_1satArray = zeros(9, n_P0);
for ii_P0 = 1:n_P0
    P0_1satArray(:, ii_P0) = diag(P0_1sat)*P0Scales(ii_P0);
end
% P0 = diag([ones(3, 1)*1e3; ones(3, 1)*1; ones(3, 1)*1e-4; ones(3, 1)*1e3; ones(3, 1)*1; ones(3, 1)*1e-4].^2);

% TODO: Come back here for process noise analysis -- abandoned when working
% on fixed EKF 
% Pqkminus0 = eye(3)*1e-4; % acceleration units [-]
% qkminus0 = ones(3, 1)*10^(-expNum(1)/2);  %  sqrt(acceleration) units [-]
% QW = diag([1e8, 1e8, 1e8]);  % acceleration units [-]
% RW = 1e-16;

Pqkminus0 = eye(3)*1e-10; % acceleration units [-]
qkminus0 = ones(3, 1)*10^(-expNum(1)/2);  %  sqrt(acceleration) units [-]
QW = diag([1e17, 1e17, 1e17]);  % acceleration units [-]
RW = 1e-1; % 1e-4 is good

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
% RW_array = [1e11, 1e11, 1e11];
RW_array =ones(1, 3)*RW^-1;
% RW_array = [1e18, 1e18, RW^-1];

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

