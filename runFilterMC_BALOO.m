%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Filter MC
%
%
% Created: 10/12/22
% Last Modified: 05/31/23
%       Modification to add in GSKF as option
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%house keeping
if ~exist('flagSetFlags', 'var')
    clear
    close all
    expNum = 14;                                       % Process noise (1e-expNum)
    underweight = 0;
    nMC = 100;
    flagUnderweight = 'fixed';
    flagEstTargetAccel = false;                       % adds target acceleration to state when true
    if flagEstTargetAccel
        flagNoiseonJerk = true;
    else
        flagNoiseonJerk =false;
    end
    flagEstChaserAccelBias = false;                   % adds chaser accelerometer bias to state space when true
    flagSampleInit = false;                            % Samples from initial distribution when true
    flagGravModel = 'j2On';                           % Sets up grav param to be on or off
    flagGravModelTruth = 'j2On';                      % Sets label for saving data with the truth data model
    flagScenario = 'KayhanRelWSE2SE'; % 'Erros';         % Sets central body for sim you want to run
    flagFilterType = 'EKF'; % EKF by default 
    flagTruthMeas = false;                            % Passes filter truth measurments when true
    flagBurnOn = true;                                % Zeros our burn hist when set to false
    flagAccelerometerPerfectMatch = true;             % When true, uses imperfect burn that matches the truth when false, matches commanded burn
    flagRunMeasUpdate = true;                         % Runs the measurment update when true
    flagInitCovStudy = false;
    fileBase = fullfile('.', 'MAT', flagScenario, flagFilterType); % Sets base file path for Figure .mat file read in
    flagQ = 'fixed'; %'fixed'; %'CovMatchTrad' % 'CorrMethod' % Sets the process noise handling
    flagConsider = true(1, 6); % [true, true,  false(1, 4)]; % Sets parameters to be consider parameters when true (lumps all 3 components together (pos, vel, accel) false(6, 1); %
    plotCheck = true;
    beta = 0*[0, 1/3600, 0]; % time constan -- 1 hour
    accelBias = zeros(6, 1); %[0; 0; 0; (300*9.81e-6); (300*9.81e-6); (300*9.81e-6)]; % initial bias (300*9.81e-6) LN200S
    gate = 1e-7;
    accelerometerVRW = zeros(3); % eye(3)*(35e-6*9.81*3)^2;       % (m/s/srt(s))^2 - VRW^2 -- LN200S for numbers 
    eom = @tbpJ2_accelerometerGating_rpod_lidarBias;
    adaptBuffer = 50;
    flagUseTrueJacobians = false;                       % Debugging flag to give EKF and adaptation true jacobians 
    P0_1sat = diag([ones(3, 1)*5; ones(3, 1)*.5; ones(3, 1)*1e-6].^2); 

    % P0 = diag([ones(3, 1)*1e3; ones(3, 1)*1; ones(3, 1)*1e-4; ones(3, 1)*1e3; ones(3, 1)*1; ones(3, 1)*1e-4].^2);
    if contains(flagQ, 'CorrMethod')
        Pqkminus0 = eye(3)*1e-7; % May want to come back to this at some point
        qkminus0 = ones(3, 1)*1e-9;% ones(3, 1)*Q0(1, 1) + sqrt((Pqkminus*randn(3, 1)).^2);
        QW = eye(3)*1e-12;
        RW= 1*1e-4;
        if contains(flagQ, 'Inf')
            Lambda0 = Pqkminus0\eye(3);
            dk0 = Lambda0*qkminus0;
            QW = QW\eye(3);
            RW = 1/RW;
        end
    end
    iiTruth = 1;
    flagCompensateMean = false;
    CRTBPFlag = false;
end

rng(18)


%% Load in Truth and Measurment Data
load(fullfile(fileBase, 'truthState.mat'))
load(fullfile(fileBase, 'meas_.mat'))
load(fullfile(fileBase, 'truemeas.mat'))
load(fullfile(fileBase, 'acaf2aci.mat'))
load(fullfile(fileBase, 'aci2body2CF1.mat'))
load(fullfile(fileBase, 'aci2body1CF2.mat'))
load(fullfile(fileBase, 'aci2body2CF1.mat'))
load(fullfile(fileBase, 'aci2body2CF2.mat'))
load(fullfile(fileBase, 'body12N.mat')) % RBody1_N'
load(fullfile(fileBase, 'body22N.mat')) % 'RBody2_N'
load(fullfile(fileBase, 'timeMat.mat'))
load(fullfile(fileBase, 'burnHistActual.mat'));
load(fullfile(fileBase, 'burnHistCommand.mat'));
load(fullfile(fileBase, 'accelerometerTruth.mat')) % 'accelerometer'
load(fullfile(fileBase, 'Rmat.mat'))
load(fullfile(fileBase, 'Hmat.mat'))
load(fullfile(fileBase, 'groundFeatIndex.mat'))
load(fullfile(fileBase, 'featureSet.mat'));
% load(fullfile(fileBase, 'impulsiveBurns.mat')); 
load(fullfile(fileBase, 'trueSTM.mat')); % 'Phi_true'
H_true = H;

%% Implement Flags

if flagSampleInit
    initName = '_noisyInit';
else
    initName = '_trueInit';
end
noUpdateCount = 0;

% CM Bounds
procNoiseName = '';
[n_nodesQ, n_exp] = size(expNum);
Q0 = zeros(3*n_exp, 1, n_nodesQ);
for ii_node = 1:n_nodesQ
    for i = 1:n_exp
        if expNum(ii_node, i) == 0
            Q0(((i*3-2):i*3), :, ii_node) = zeros(3, 1);
            procNoiseName = strcat(procNoiseName, '_zeroQ');
        else
            Q0(((i*3-2):i*3), :, ii_node) = ones(3, 1)*10^(-expNum(ii_node, i)); % 1e-6 is good for no accel -- 1e-12 is good for accel on
            procNoiseName = strcat('_nonzeroQ_eneg', num2str(expNum(ii_node, i)));
        end
    end
end



%% TODO: Implement other measurment types
% Only the relative navigation version of range, az, el has been coded up
% to perfrom in a page setting


[measNum, ~] = size(meas);
if measNum == 5 % have NOT implemented this 
    measFunc = @relNavRngRngRate_opTarget_feature_meas_H;
    R_rangeRangeRate_chol = sqrt(R(1, 1));
    R_opNav_chol = chol(R(2:3, 2:3), 'lower');
    measFlag = 'relRngRelAzElFeatAzEl';
elseif measNum == 3 
    if contains(flagScenario, 'direct') % have implemented this
        measFunc = @relNavDirectInertial_meas_H_page;
        R_direct_chol = chol(R, 'lower');
        measFlag = 'direct';
    else % have implemented this
        measFunc = @relNavRng_opTarget_meas_H_page;
        R_rangeRangeRate_chol = sqrt(R(1, 1, 1));
        R_opNav_chol = chol(R(2:3, 2:3, 1), 'lower');
        measFlag = 'relRngRelAzEl';
        if ~contains(flagScenario, 'Rel') % have NOT implemented this
            measFunc = @rng_opTarget_meas_H_page;
        end
    end
else % have NOT implemented this
    measFunc = @rngRngRate_opTarget_Ground_meas_H; 
    R_rangeRangeRate_chol = chol(R(1:2, 1:2), 'lower');
    R_opNav_chol = chol(R(3:4, 3:4), 'lower');
    measFlag = 'relRngRelAzElFeatAzEl';
end

if contains(flagScenario, 'Rel')
    relNavFlag = true;
    if flagInitCovStudy || norm(P0_1sat) == 0
        P0_1satTarget = P0_1sat*2;
        P0 = [P0_1satTarget, P0_1sat; P0_1sat, P0_1sat];
        % P0(16:18, 16:18) = P0(16:18, 16:18)*1e-8;
        P0(7:9, 16:18) = 0;
        P0(16:18, 7:9) = 0;
    else
        P0 = [P0_1sat*2, P0_1sat; P0_1sat, P0_1sat];
    end
else
    relNavFlag = false;
    P0 = blkdiag(P0_1sat, P0_1sat);
end




if flagLidarBias
    P0 = blkdiag(P0, 0.006);
end


if contains(flagScenario, 'Erros')
    % Erros
    mu = 4.4627e5; % m^3/s^2 -- From: A Global Solution for the Gravity Field, Rotation, Landmarks, and Ephemeris of Eros
    J2 = 0.052478; % From: Determination of Shape, Gravity, and Rotational State of Asteroid 433 Eros
    J3 = -0.0014; 
    J4 = -0.012134;
    J5 = 0;
    J6 = 0;
    J7 = 0; 
    J8 = 0;
    RBody = 16e3; % m % From: Determination of Shape, Gravity, and Rotational State of Asteroid 433 Eros
else
    mu = 3.986e5*(1e3)^3; % m^3/s^2 
    J2 = .0010826269; 
    J3 = -0.0000025323;
    J4 = -0.0000016204;
    J5 = -2.27753590730836e-7;
    J6 = 5.40666576283813e-7;
    J7 = -3.50551795713742e-7;
    J8 = -2.03993125929884e-7;
    RBody = 6378*1e3; % m 
end

CRTBPFlag = false;
if strcmp(flagGravModel, '8x8On')
    gravParam = [mu, RBody, J2, J3, J4, J5, J6, J7, J8, 1];
elseif strcmp(flagGravModel, 'j2On')
    gravParam = [mu, RBody, J2, 0, 0, 0, 0, 0, 0, 0];
elseif strcmp(flagGravModel, 'j2j3On')
    gravParam = [mu, RBody, J2, J3, 0, 0, 0, 0, 0, 0];
elseif strcmp(flagGravModel, 'j2j3j4On')
    gravParam = [mu, RBody, J2, J3, J4, 0, 0, 0, 0, 0];
elseif strcmp(flagGravModel, 'j2j3j4j5On')
    gravParam = [mu, RBody, J2, J3, J4, J5, 0, 0, 0, 0];
elseif strcmp(flagGravModel, 'j2j3j4j5j6On')
    gravParam = [mu, RBody, J2, J3, J4, J5, J6, 0, 0, 0];
elseif strcmp(flagGravModel, 'j2j3j4j5j6j7On')
    gravParam = [mu, RBody, J2, J3, J4, J5, J6, J7, 0, 0];
elseif strcmp(flagGravModel, 'j2j3j4j5j6j7j8On')
    gravParam = [mu, RBody, J2, J3, J4, J5, J6, J7, J8, 0];
elseif contains(flagScenario, 'NRHO')
    mu = 0.012150584270; % mass ratio from paper
    L = 384400; % km -- characteristic length 
    tCharacteristic = 3.751902619517e5; % s -- characteristic time from paper
    n_motion = 2*pi/P; % 1/s 
    gravParam = [mu, L, n_motion, tCharacteristic];
    CRTBPFlag = true;
end


considerParam = false(19, 1);
for i = 1:6
    if flagConsider(i)
        considerParam((i*3-2):(i*3)) = true;
    end
end
if flagConsider(7)
    considerParam(19) = true;
end

if flagAccelerometerPerfectMatch
    burnName = '_noNoiseAccelerometer';
else
    burnName = '_noisyAccelerometer';
end

burnHist = accelerometer;

if ~flagBurnOn
    burnHist = burnHistCommand*0;
    burnName = '_noBurn';
end


if ~flagEstTargetAccel && ~flagEstChaserAccelBias
    accelName = 'accelOffBiasOff';
    if flagLidarBias
        P0 = P0([1:6, 10:15, 19], [1:6, 10:15, 19]);
        truthState = truthState([1:6, 10:15, 19], :);
        [n, ~] = size(P0);
        H = H(:, [1:6, 10:15, 19], :);
        considerParam = considerParam([1:6, 10:15, 19]);
        Phi_true = Phi_true([1:6, 10:15, 19], [1:6, 10:15, 19], :);
        H_true = H_true(:, [1:6, 10:15, 19], :);
    else
        P0 = P0([1:6, 10:15], [1:6, 10:15]);
        truthState = truthState([1:6, 10:15], :);
        [n, ~] = size(P0);
        H = H(:, [1:6, 10:15], :);
        considerParam = considerParam([1:6, 10:15]);
        Phi_true = Phi_true([1:6, 10:15], [1:6, 10:15], :);
        H_true = H_true(:, [1:6, 10:15], :);
    end
    targetStates = 1:6;
    chaserStates = 7:12;
    accelBias = mean(accelerometer(:, 1:3600), 2);
    n_elm = [6, 6];
elseif ~flagEstTargetAccel && flagEstChaserAccelBias
    if flagLidarBias
        P0 = P0([1:6, 10:19], [1:6, 10:19]);
        truthState = truthState([1:6, 10:19], :);
        [n, ~] = size(P0);
        H = H(:, [1:6, 10:19], :);
        considerParam = considerParam([1:6, 10:19]);
        Phi_true = Phi_true([1:6, 10:19], [1:6, 10:19], :);
        H_true = H_true(:, [1:6, 10:19], :);
    else
        P0 = P0([1:6, 10:18], [1:6, 10:18]);
        truthState = truthState([1:6, 10:18], :);
        [n, ~] = size(P0);
        H = H(:, [1:6, 10:18], :);
        considerParam = considerParam([1:6, 10:18]);
        Phi_true = Phi_true([1:6, 10:18], [1:6, 10:18], :);
        H_true = H_true(:, [1:6, 10:18], :);
    end
    targetStates = 1:6;
    chaserStates = 7:15;
    accelName = 'accelOffBiasOn';
    if N < 3600
        accelBias = mean(accelerometer(:, 1:N), 2);
    else
        accelBias = mean(accelerometer(:, 1:3600), 2);
    end
    n_elm = [6, 9];
elseif flagEstTargetAccel && ~flagEstChaserAccelBias
    if flagLidarBias
        P0 = P0([1:15, 19], [1:15, 19]);
        truthState = truthState([1:15, 19], :);
        [n, ~] = size(P0);
        H = H(:, [1:15, 19], :);
        considerParam = considerParam([1:15, 19]);
        Phi_true = Phi_true([1:15, 19], [1:15, 19], :);
        H_true = H_true(:, [1:15, 19], :);
    else
        P0 = P0(1:15, 1:15);
        truthState = truthState(1:15, :);
        [n, ~] = size(P0);
        H = H(:, 1:15, :);
        considerParam = considerParam(1:15);
        Phi_true = Phi_true(1:15, 1:15, :);
        H_true = H_true(:, 1:15, :);
    end
    targetStates = 1:9;
    chaserStates = 10:15;
    accelName = 'accelOnBiasOff';
    accelBias = mean(accelerometer(:, 1:3600), 2);
    n_elm = [9, 6];
else
    if ~flagLidarBias
        P0 = P0(1:18, 1:18);
        truthState = truthState(1:18, :);
        [n, ~] = size(P0);
        H = H(:, 1:18, :);
        considerParam = considerParam(1:18);
        Phi_true = Phi_true(1:18, 1:18, :);
        H_true = H_true(:, 1:18, :);
    else
        n = length(P0);
    end
    targetStates = 1:9;
    chaserStates = 10:18;
    accelName = 'accelOnBiasOn';
    n_elm = [9, 9];
end

if flagLidarBias
    accelName = strcat(accelName, 'LidarBiasOn');
else
    accelName = strcat(accelName, 'LidarBiasOff');
%     eom = @tbpJ2_accelerometerGating_rpod;
end

if any(diag(P0)==0)
    P0chol = P0*0;
    P0chol(targetStates, targetStates) = chol(P0(targetStates, targetStates), 'lower');
else
    P0chol = chol(P0, 'lower');
end
if norm(accelerometerVRW) == 0
    accelChol = accelerometerVRW;
else
    accelChol = chol(accelerometerVRW, 'lower'); % m/s^2
end


if contains(flagScenario, 'interacting')
    flagInteracting = true;
else
    flagInteracting = false;
end

if contains(flagScenario, 'BALOO')
    flagBALOOsplit = true;
else
    flagBALOOsplit = false;
end

if contains(flagScenario, '2modelsNoiseDifference')
    flagIMMnoiseDiff = true;
else
    flagIMMnoiseDiff = false;
end

if flagInteracting && ~flagBALOOsplit
    if ~flagIMMnoiseDiff
        P0(7:9, 7:9) = 1e-36;
    end
end

N = length(truthState);
NEES = NaN(nMC, N);
NEESdiag = NaN(nMC, N);
NEESfirst6 = NaN(nMC, N);
NEESfirst6diag = NaN(nMC, N);
[m, ~, N_r] = size(R);

if N_r > 1
    R_all = R;
    R = R_all(:, :, 1);
end
Nwindow = 30;

%% Run EKF
for i = 1:nMC
    %% Initialize Variables
    % Initialize output structures that contain standard single state
    % output and covariance for combined nodes
    % Initialize new set of variables with structure for layers, always 
    % initialize to just one node, that will change if there are new nodes
    % added. This is the set of variables that things will funciton on.
    rng(211105+i+iiTruth+1)
    % Sample initial condition
    if flagSampleInit
        dx = randn(n, 1);
        % fprintf('initial dx %f, %f, %f\n', dx(7:9))
        % fprintf('initial dx %f, %f, %f\n', dx(10:12))
        x0 = truthState(:, 1) + P0chol*dx;
        if flagEstTargetAccel
            x0(7:9, 1) = 0;
        end
%         if flagLidarBias
%             x0(end, 1) = biasLidar;
%         end
    else
        x0 = truthState(:, 1); 
    end
    
    % Output variables -- combine the nodes together to export solution.
    x = NaN(n, N);
    xminus = NaN(n, N);
    P = NaN(n, n, N);
    Pminus = NaN(n, n, N);
    riMinus = NaN(m, N);
    riPlus = NaN(m, N);
    Pri = NaN(m, m, N);
    PriMinus = NaN(m, m, N);
    qCountDiag = NaN(3, N);
    qCountDiagCov = NaN(3, 3, N);
    
    %% TODO: Revert this section to initialize with just 1 node
    if flagInteracting && ~flagBALOOsplit
        nodeIDs = 1:n_nodes;
    else
        n_nodes = 1;
        nodeIDs = 1;
    end
    % Working variables -- contain an additional dimension to keep track of
    % all nodes -- initialized as having one node
    x_nodes = NaN(n, n_nodes, N);
    P_nodes = NaN(n, n, n_nodes, N);
    Pminus_nodes = NaN(n, n, n_nodes, N);
    riMinus_nodes = NaN(m, n_nodes, N);
    riPlus_nodes = NaN(m, n_nodes, N);
    Pri_nodes = NaN(m, m, n_nodes, N);
    PriMinus_nodes = NaN(m, m, n_nodes, N);
    qCountDiag_nodes = NaN(3, n_nodes, N);
    qCountDiagCov_nodes = NaN(3, 3, n_nodes, N);
    vBarHist = zeros(3, n_nodes,  N); % 3x1 vector for 1 node over time
    Phi = NaN(n, n, n_nodes, N); % nxn matrix for 1 node over time
    Phi_window = NaN(n, n, n_nodes, Nwindow); % nxn matrix for 1 node over window 
    H_hist = NaN(m, n, n_nodes, N); % mxn matrix for 1 node over time 
    w_nodes = ones(n_nodes, N)/n_nodes; % Filter bank weights, weight of 1 for one node over time
    splitWait = 75;

    % All nodes over time -- size up for each new node, never sizes down 
    % -- contain an additional dimension to keep track of
    % all nodes -- initialized as having one node
    x_nodes_all = NaN(n, n_nodes, N);
    P_nodes_all = NaN(n, n, n_nodes, N);
    Pminus_nodes_all = NaN(n, n, n_nodes, N);
    riMinus_nodes_all = NaN(m, n_nodes, N);
    riPlus_nodes_all = NaN(m, n_nodes, N);
    Pri_nodes_all = NaN(m, m, n_nodes, N);
    PriMinus_nodes_all = NaN(m, m, n_nodes, N);
    qCountDiag_nodes_all = NaN(3, n_nodes, N);
    qCountDiagCov_nodes_all = NaN(3, 3, n_nodes, N);
    vBarHist_all = zeros(3, n_nodes,  N); % 3x1 vector for 1 node over time
    Phi_all = NaN(n, n, n_nodes, N); % nxn matrix for 1 node over time
    Phi_window_all = NaN(n, n, n_nodes, Nwindow); % nxn matrix for 1 node over window 
    H_hist_all = NaN(m, n, n_nodes, N); % mxn matrix for 1 node over time 
    w_nodes_all = ones(n_nodes, N)/n_nodes; % Filter bank weights, weight of 1 for one node over time


    % initialize output variables 
    x(:, 1) = x0;
    P(:, :, 1) = P0;
    Pminus(:, :, 1) = P0;
    xminus(:, 1) = x0;


    % Split variabes
    if flagBALOOsplit || ~contains(flagScenario, 'fixed')
        rWindow_nodes = NaN(m*adaptBuffer+1, n_nodes);
        yWindow = NaN(m*adaptBuffer+1, 1);
        % uWindow_nodes = zeros(3*(adaptBuffer +1), n_nodes);
        nA = length(targetStates);
        Fwindow_nodes = zeros(nA, nA, n_nodes, adaptBuffer);
        Hwindow_nodes = zeros(m*(adaptBuffer+1), nA, n_nodes);
    end
    

    % initialize working variables 
    for ii_node = 1:n_nodes
        x_nodes(:, ii_node, 1) = x0;
        P_nodes(:, :, ii_node, 1) = P0;
        Pminus_nodes(:, :, ii_node,  1) = P0;
        Phi(:, :, ii_node, 1) = eye(n);
        x_nodes_all(:, ii_node, 1) = x0;
        P_nodes_all(:, :, ii_node, 1) = P0;
        Pminus_nodes_all(:, :, ii_node,  1) = P0;
        Phi_all(:, :, ii_node, 1) = eye(n);
    end

    if ~flagBALOOsplit && flagInteracting
        if ~flagIMMnoiseDiff
            for ii_node = 1:n_nodes
                x_nodes(7:9, ii_node, 1) = accelModels(:, ii_node);
            end
        end
    end

    %TODO: Come back here to set this with some more thought

    if ~exist('epsilon_split', 'var')
        epsilon_split = 10; 
    end
    if ~exist('epsilon_prune', 'var')
        epsilon_prune = -1; 
    end

    num_nodesPlot = ones(1, N);

    %TODO: Bring Buffer building out of CorrMethods block so I can use
    %buffer to split plain EKF later? 

    
    iminus1 = 1;
    accelerometerMeas = burnHist;

    dt = time(2) - time(1);
    flagIncMean = false;
    if dt > 1
        dtProp = 1;
    else
        dtProp = dt;
    end

    Qhat = zeros(n, n, n_nodes);
    if n_elm(1) > 6
        Xk = Xk_relInertOnly_jerkNoise(time(2) - time(1));
        if contains(flagScenario, 'trying')
            Xk = Xk_relInertOnly_accelNoise(time(2) - time(1));
        end
    else
        [Xk] = Xk_relInertOnly_ZOH(time(2) - time(1));
    end
    Xk_nodes = repmat(Xk, 1, 1, n_nodes);
    QhatUT = pagemtimes(Xk_nodes, Q0(1:3, :, :));
    QhatTemp = reverseVech_page(QhatUT);
    Qhat(targetStates, targetStates, :) = QhatTemp;
    QhatUT = pagemtimes(Xk_nodes, Q0(4:6, :, :));
    QhatTemp = reverseVech_page(QhatUT);
    Qhat(chaserStates, chaserStates, :) = QhatTemp;
    
    if flagInteracting && ~flagBALOOsplit
        if ~flagIMMnoiseDiff
            Qhat(7:9, 7:9) = 0;
        end
    end
    

    if flagEstChaserAccelBias
        Qhat(chaserStates(7:9), chaserStates(7:9), :) = repmat(accelerometerVRW*dt/100, 1, 1, n_nodes);
    end

    if length(Qhat) < n
        Qhat = blkdiag(Qhat, 0);
    end
    %% Process Noise Stuff 
    if strcmp(flagQ, 'CovMatchTrad') % This case isn't currently oppertaitonal 4/3/25
        Q_hist = zeros(n, n, Nwindow);
    elseif strcmp(flagQ, 'CovMatch') % This case isn't currently oppertaitonal 4/3/25
        % D'Amico paper 
        [Xk] = Xk_ADMC(beta(1), time(2) - time(1));
        [n_w, ~] = size(diag(vech(eye(n))));
        W_hist = zeros(n_w, n_w, Nwindow);
    elseif contains(flagQ, 'CorrMethod') % in progress
        if n_elm(1) > 6
            Xk = Xk_relInertOnly_jerkNoise(time(2) - time(1));
            if contains(flagScenario, 'trying')
                Xk = Xk_relInertOnly_accelNoise(time(2) - time(1));
            end
        else
            [Xk] = Xk_relInertOnly_ZOH(time(2) - time(1));
        end
        Pqkminus = Pqkminus0; % Gets updated each time adaptive method employed -- easy to change to be layered
   
        
        deviationFromqkminus = ones(3, 1)*10^(-randi([4, 11], 1, 1));%chol(Pqkminus)*randn(3, 1);
        qkminus = qkminus0 + deviationFromqkminus;
        QhatUT = Xk*(qkminus.^2);
        QhatTemp = reverseVech(QhatUT);
        Qhat(targetStates, targetStates) = QhatTemp;
        if contains(flagQ, 'Inf')
            % Generate initial conditions for information form 
            Lambdakminus = Lambda0;
            dkminus = Lambdakminus*qkminus;
            % Overwrites cov form for later use 
            qkminus = dkminus;
            Pqkminus = Lambdakminus;
        end

        %TODO: Must remove this chunk when done testing
        qkminus = repmat(qkminus, 1, n_nodes);
        Qhat = repmat(Qhat, 1, 1, n_nodes); % each node has its own Qhat
        Lambdakminus = repmat(Lambdakminus, 1, 1, n_nodes); % each node has its own Informaiton matrix 
        dkminus = repmat(dkminus, 1, n_nodes); % each node has its own information state 
        Pqkminus = repmat(Pqkminus, 1, 1, n_nodes); % each node has its own covariance
        Xk_nodes = repmat(Xk, 1, 1, n_nodes); % repeat for ease later -- will need this line when resizing number of nodes later 
        
        % rWindow = NaN(m, adaptBuffer+1);
        % uWindow = zeros(3, adaptBuffer +1);
        % nA = length(targetStates);
        % Fwindow = zeros(nA, nA, adaptBuffer);
        % Hwindow = zeros(m, nA, adaptBuffer+1);
        % Awindow = zeros(nA, m*(adaptBuffer+1), adaptBuffer+1);
        % Bwindow = zeros(nA, nA, adaptBuffer, adaptBuffer);
        % Zwindow = zeros(nA, adaptBuffer+1);
        % GWindow = [];

        % rWindow_nodes = NaN(m*adaptBuffer+1, n_nodes);
        % yWindow = NaN(m*adaptBuffer+1, 1);
        uWindow_nodes = zeros(3*(adaptBuffer +1), n_nodes);
        nA = length(targetStates);
        % Fwindow_nodes = zeros(nA, nA, n_nodes, adaptBuffer);
        % Hwindow_nodes = zeros(m*(adaptBuffer+1), nA, n_nodes);
        Awindow_nodes = zeros(nA*(adaptBuffer+1), m*(adaptBuffer+1), n_nodes);
        Bwindow_nodes = zeros(nA*(adaptBuffer), nA*(adaptBuffer), n_nodes);
        Zwindow_nodes = zeros(nA*(adaptBuffer+1), n_nodes);
        GWindow_nodes = [];

        fullEnoughAdaptCheck = false(1, n_nodes);
        fullEnoughSplitCheck = false(1, n_nodes);
        if ~exist('fullEnoughSplit', 'var')
            fullEnoughSplit = 15;
        end
        if ~exist('fullEnoughAdapt', 'var')
            fullEnoughAdapt = 15;
        end
    end
    %%
    meas_run = NaN(m, N);
    num_nodesPlot = ones(1, N)*n_nodes;
    jminus1 = 1;
    accelBias = zeros(3, 1);
    fprintf('Running EKF MC run %d...\n', iiTruth)
    for j = 2:N
        % fprintf('%d\n', j)
        if j == 2916
            hfjewhfew = 0;
        end
        %% Interacting?
        if flagInteracting
            % Making the choice (for now??) to not mix the process noise
            % estimates in an interacting multiple MILTON scenario -- I
            % don't think it makes sesne to have the process noise
            % estimates interact -- at least not at first? Probably makes
            % sense anyway because one difference in IMM that is possible
            % is to just have different process noise values 
            if n_nodes > 1
                if flagBALOOsplit % if it's not BALOO, it won't need this on the fly b/c the numner of nodes is fixed 
                    if exist('staysTheSame', 'var')
                        [transitionProbMat] = transitionalProbabilityMatIMM(n_nodes, staysTheSame);
                    else
                        [transitionProbMat] = transitionalProbabilityMatIMM(n_nodes);
                    end
                end
                mixingProb = transitionProbMat*w_nodes(:, jminus1);
                mixingProb = mixingProb/sum(mixingProb);
                [xMixjminus1, PMixjminus1] = outputGMM(x_nodes(:, :, jminus1), P_nodes(:, :, :, jminus1), mixingProb, n_nodes);
                if flagIMMnoiseDiff
                    % Want to mix acceleration estimates when the "model
                    % Difference" is the npoise level in the system
                    xMixjminus1 = repmat(xMixjminus1, 1, n_nodes);
                    PMixjminus1 = repmat(PMixjminus1, 1, 1, n_nodes);
                    x_nodes(:, :, jminus1) = xMixjminus1;
                    P_nodes(:, :, :, jminus1) = PMixjminus1;
                else
                    xMixjminus1 = repmat(xMixjminus1, 1, n_nodes);
                    PMixjminus1 = repmat(PMixjminus1, 1, 1, n_nodes);
                    x_nodes(1:6, :, jminus1) = xMixjminus1(1:6, :);
                    % x_nodes(chaserStates, :, jminus1) = xMixjminus1(chaserStates, :);
                    P_nodes(1:6, 1:6, :, jminus1) = PMixjminus1(1:6, 1:6, :);
                    % w_nodes(:, jminus1) = mixingProb;
                end
            end
        end

        %% Propagation 
        tj = time(j);
        times = [time(jminus1), tj];
        if ~flagAccelerometerPerfectMatch
            if jminus1 == 1
                accelerometerMeas(:, jminus1) = accelerometerMeas(:, jminus1) +  accelChol*randn(3, 1)*sqrt(dtProp);
                accelerometerMeas(:, j) = accelerometerMeas(:, j) +  accelChol*randn(3, 1)*sqrt(dtProp);
            else
                accelerometerMeas(:, j) = accelerometerMeas(:, j) +  accelChol*randn(3, 1)*sqrt(dtProp);
            end
        end
        if norm(accelerometerMeas(:, j)) > gate
            % fprintf('Burn Detected, time step %d:%d\n', jminus1:j)
           if ~flagEstChaserAccelBias
               accelBias = mean(accelerometerMeas(:, (j-121):(j-1)), 2); % body -- averaging 2000 non maneuverg points to get bias
           end
        end
        Phi(:, :, :, j) = repmat(eye(n), 1, 1, n_nodes);
        if j ~=1
            if j == 119
                fjdkl  =0;
            end
            timeProp = times(1):dtProp:times(2);
            Nprop = length(timeProp)-1;
            state = x_nodes(:, :, jminus1);
            % P(:, :, j) = P(:, :, jminus1);
            P_nodes(:, :, :, j) = P_nodes(:, :, :, jminus1);
            Rbody_N = RBody2_N(:, :, j);
            for k = 1:Nprop
                initialState = state;
                if contains(flagScenario, 'NRHO')
                    [state, P_nodes(:, :, :, j), Phi_k] = GMM_EKFprop(timeProp(k:k+1), initialState, P_nodes(:, :, :, j), eom, Qhat, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, vBarHist(:, :, j-1), relNavFlag, zeros(3, 1), n_elm, accelerometerVRW*dt, flagUseTrueJacobians, squeeze(Phi_true(:, :, j)), CRTBPFlag);
                    %                            EKFprop(times,           x,            P,          eom, Qhat, gravParam, beta_markov, accelerometer,    gate, Rbody_N, RN_rot,           relNavFlag, accelBias, n_elm, accelVRW,            flagUseTrueJacobians, Phi_true,                   CRTBPFlag)
                else
                    % weight should remain unchanged through prop if using
                    % GMM 
                    [state, P_nodes(:, :, :, j), Phi_k] = GMM_EKFprop(timeProp(k:k+1), initialState, P_nodes(:, :, :, j), eom, Qhat, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, vBarHist(:, :, j-1), relNavFlag, accelBias, n_elm, accelerometerVRW*dt, flagUseTrueJacobians, squeeze(Phi_true(:, :, j)));
                    % eom, Qhat, gravParam, beta_markov, accelerometer, gate, Rbody_N, RN_rot,relNavFlag, accelBias, n_elm
                end
                Phi(:, :, :, j) = pagemtimes(Phi_k, Phi(:, :, :, j));
            end
            Pminus_nodes(:, :, :, j) = P_nodes(:, :, :, j);
            x_nodes(:, :, j) = state;

            x_nodes_all(:, nodeIDs, j) = x_nodes(:, :, j);
            Pminus_nodes_all(:, :, nodeIDs, j) = Pminus_nodes(:, :, :, j);
            Phi_all(:, :, nodeIDs, j) = Phi(:, :, :, j);
            
            % weight should remain unchanged through prop
            % cobmining nodes for output
            [xminus(:, j), Pminus(:, :, j)] = outputGMM(x_nodes(:, :, j), Pminus_nodes(:, :, :, j), w_nodes(:, jminus1), n_nodes);
            P(:, :, j) = Pminus(:, :, j);
            x(:, j) = xminus(:, j);

        end
        if flagRunMeasUpdate
            %% Measurment Step
            % Add noise to measurments
            if N_r > 1
                if (R(1, 1) - R_all(1, 1, j) ~=0) 
                    R = R_all(:, :, j);
                    R_rangeRangeRate_chol = sqrt(R(1, 1));
                    R_opNav_chol = chol(R(2:3, 2:3), 'lower');
                end
            end
            meas_j = meas_truth(:, j);
            if ~flagTruthMeas
                if strcmp(measFlag, 'relRngRelAzEl')
                    meas_j(1) = meas_j(1) + R_rangeRangeRate_chol*randn(1, 1);
                    meas_j(2:3) = meas_j(2:3) + R_opNav_chol*randn(2, 1);
                elseif strcmp(measFlag, 'relRngRelAzElFeatAzEl')
                    meas_j(1) = meas_j(1) + R_rangeRangeRate_chol*randn(1, 1);
                    meas_j(2:3) = meas_j(2:3) + R_opNav_chol*randn(2, 1);
                    meas_j(4:5) = meas_j(4:5) + R_opNav_chol*randn(2, 1);
                elseif strcmp(measFlag, 'direct')
                    meas_j = meas_j + R_direct_chol*randn(3, 1);
                else
                    meas_j(1:2) = meas_j(1:2) + R_rangeRangeRate_chol*randn(2, 1);
                    meas_j(3:4) = meas_j(3:4) + R_opNav_chol*randn(2, 1);
                    if m > 4
                        meas_j(5:6) = meas_j(5:6) + R_opNav_chol*randn(2, 1);
                    end
                end
            end
            meas_run(:, j) = meas_j;
            % if j < 10 
            %     fprintf('%f, %f, %f\n', (meas_run(:, j) - meas_truth(:, j))./ [diag(R_rangeRangeRate_chol); diag(R_opNav_chol)])
            % end
            if isnan(featureIndex(j))
                L_j = NaN(3, 1);
            else
                L_j = L_set(:, featureIndex(j));
            end


            if flagBALOOsplit || contains(flagQ, 'CorrMethod')
                % Generate prefit measurment residuals
                [riMinus_nodes(:, :, j), PriMinus_nodes(:, :, :, j), H] = GMM_EKFmeas_prefitOutput(x_nodes(:, :, j), P_nodes(:, :, :, j), meas_j, R, measFunc, L_j, RCF1_Body2_N(:, :, j), RCF2_Body2_N(:, :, j), vBarHist(:, :, j-1), n_elm, measFlag, flagUseTrueJacobians, squeeze(H_true(:, :, j)));

                % Add prefit information to buffers
                H_hist(:, :, :, j) = H(1:m, :, :);

                rWindow_nodes((m+1):(m*(adaptBuffer+1)), :) =  rWindow_nodes(1:(m*(adaptBuffer)), :); % [riMinus_nodes(:, :, j), rWindow_nodes(:, :, 1:end-1)];
                rWindow_nodes(1:m, :) = riMinus_nodes(:, :, j);

                yWindow((m+1):(m*(adaptBuffer+1))) =  yWindow(1:(m*(adaptBuffer)));
                yWindow(1:m, :) = meas_j;

                Fwindow_nodes(:, :, :, 2:end) = Fwindow_nodes(:, :, :, 1:end-1);
                Fwindow_nodes(:, :, :, 1) = Phi(targetStates, targetStates, :, j); %Phi_true(targetStates, targetStates, j); %

                Hwindow_nodes((m*(2-1)+1):end, :, :) = Hwindow_nodes(1:(m*(adaptBuffer)), :, :); % move other H matricies further down
                Hwindow_nodes(1:m, :, :) = H_hist(1:m, targetStates, :, j); %  overwrite the first H in buffer

                % update full enough criteria
                NaNcheck = sum(isnan(rWindow_nodes));
                fullEnoughAdaptCheck = NaNcheck <= (m*adaptBuffer - (fullEnoughAdapt*m));
                fullEnoughSplitCheck = false(1, n_nodes);
                % fullEnoughSplitCheck(splitWait > 0) = false;
                fullEnoughAdaptCheck(w_nodes(:, jminus1) < 1e-10) = false;
                fullEnoughSplitCheck = splitWait <= 0;
                if flagBALOOsplit
                    if contains(flagScenario, 'fixed') && ~contains(flagScenario, 'nonDim')
                        [x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, n_nodes, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait] = BALOOnoMILTON_splitCheck(fullEnoughSplit, epsilon_split, adaptBuffer, time, x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, yWindow, H_hist, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait, RCF1_Body2_N, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, fullEnoughSplitCheck, eom, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, vBarHist(:, :, j-1), relNavFlag, accelBias, n_elm, R);
                    elseif contains(flagScenario, 'fixed') && contains(flagScenario, 'nonDim')
                        [x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, n_nodes, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait] = BALOOnoMILTON_splitCheck_NRHO(fullEnoughSplit, epsilon_split, adaptBuffer, time, x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, yWindow, H_hist, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait, RCF1_Body2_N, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, fullEnoughSplitCheck, eom, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, vBarHist(:, :, j-1), relNavFlag, accelBias, n_elm, R);
                    elseif contains(flagScenario, 'nonDim')
                        [x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, n_nodes, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait,   qkminus, Pqkminus, Awindow_nodes, Bwindow_nodes, Zwindow_nodes, vBarHist, uWindow_nodes, Xk_nodes, Qhat] = BALOO_splitCheck_NRHO(fullEnoughSplit, epsilon_split, adaptBuffer, time, x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, yWindow, H_hist, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait, RCF1_Body2_N, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, fullEnoughSplitCheck, eom, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, vBarHist(:, :, j-1), relNavFlag, accelBias, n_elm,  qkminus, Pqkminus, qkminus0, Pqkminus0, Awindow_nodes, Bwindow_nodes, Zwindow_nodes, vBarHist, uWindow_nodes, Xk, Qhat, R);
                    else
                        [x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, n_nodes, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait,   qkminus, Pqkminus, Awindow_nodes, Bwindow_nodes, Zwindow_nodes, vBarHist, uWindow_nodes, Xk_nodes, Qhat] = BALOO_splitCheck(fullEnoughSplit, epsilon_split, adaptBuffer, time, x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, yWindow, H_hist, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait, RCF1_Body2_N, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, fullEnoughSplitCheck, eom, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, vBarHist(:, :, j-1), relNavFlag, accelBias, n_elm,  qkminus, Pqkminus, qkminus0, Pqkminus0, Awindow_nodes, Bwindow_nodes, Zwindow_nodes, vBarHist, uWindow_nodes, Xk, Qhat, R);
                        % [x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, n_nodes, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait,   qkminus, Pqkminus, Awindow_nodes, Bwindow_nodes, Zwindow_nodes, vBarHist, uWindow_nodes, Xk_nodes, Qhat] = BALOO_splitCheck_threeSplit(fullEnoughSplit, epsilon_split, adaptBuffer, time, x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, yWindow, H_hist, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait, RCF1_Body2_N, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, fullEnoughSplitCheck, eom, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, vBarHist(:, :, j-1), relNavFlag, accelBias, n_elm,  qkminus, Pqkminus, qkminus0, Pqkminus0, Awindow_nodes, Bwindow_nodes, Zwindow_nodes, vBarHist, uWindow_nodes, Xk, Qhat, R);
                    end
                    % x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all
                    % num_nodesPlot(j) = n_nodes;
                    NaNcheck = sum(isnan(rWindow_nodes));
                    fullEnoughAdaptCheck = NaNcheck <= (m*adaptBuffer - (fullEnoughAdapt*m));
                    % fullEnoughSplitCheck = false(1, n_nodes);
                    % fullEnoughSplitCheck(splitWait > 0) = false;
                    fullEnoughAdaptCheck(w_nodes(:, jminus1) < 1e-10) = false;
                    % fullEnoughSplitCheck(splitWait <= 0) = true;
                    fullEnoughSplitCheck = splitWait <= 0;
                end

            end
            num_nodesPlot(j) = n_nodes;
            % Run measurment update
            [x_nodes(:, :, j), P_nodes(:, :, :, j), w_nodesPlus, riPlus_nodes(:, :, j), riMinus_nodes(:, :, j), Pri_nodes(:, :, :, j), PriMinus_nodes(:, :, :, j), H] = GMM_EKFmeas(x_nodes(:, :, j), P_nodes(:, :, :, j), w_nodes(:, jminus1), meas_j, R, measFunc, L_j, RCF1_Body2_N(:, :, j), RCF2_Body2_N(:, :, j), vBarHist(:, :, j-1), considerParam, n_elm, measFlag, underweight, accelerometerMeas(:, j), gate, accelerometerVRW*dt, flagUnderweight, flagUseTrueJacobians, squeeze(H_true(:, :, j)), flagInteracting);
            
            % Calculating output values
            [riMinus(:, j), PriMinus(:, :, j)] = outputGMM(riMinus_nodes(:, :, j), PriMinus_nodes(:, :, :, j), w_nodes(:, jminus1), n_nodes);
            [riPlus(:, j), Pri(:, :, j)] = outputGMM(riPlus_nodes(:, :, j), Pri_nodes(:, :,:,  j), w_nodesPlus, n_nodes);
            [x(:, j), P(:, :, j)] = outputGMM(x_nodes(:, :, j), P_nodes(:, :, :, j), w_nodesPlus, n_nodes);
            
            w_nodes(:, j) = w_nodesPlus;

            

            x_nodes_all(:, nodeIDs, j) = x_nodes(:, :, j);
            P_nodes_all(:, :, nodeIDs, j) = P_nodes(:, :, :, j);
            w_nodes_all(nodeIDs, j) = w_nodesPlus;


            if flagBALOOsplit
                % Overwrite H -- should be the same if no new nodes
                H_hist(:, :, :, j) = H(1:m, :, :);
                % Prune Check
                if n_nodes > 1
                    if contains(flagScenario, 'fixed')
                        [x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, n_nodes, nodeIDs, splitWait] = pruneZeroWeightNodes_BALOOnoMILTON(epsilon_prune, x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, nodeIDs, splitWait, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes);
                    else
                        [x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, n_nodes, nodeIDs, splitWait, qkminus, Pqkminus,  Awindow_nodes, Bwindow_nodes, Zwindow_nodes, vBarHist, uWindow_nodes, Xk_nodes, Qhat] = pruneZeroWeightNodes(epsilon_prune, x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, nodeIDs, splitWait, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, qkminus, Pqkminus, Awindow_nodes, Bwindow_nodes, Zwindow_nodes, vBarHist, uWindow_nodes, Xk, Qhat);
                    end
                end
                
                NaNcheck = sum(isnan(rWindow_nodes));
                fullEnoughAdaptCheck = NaNcheck <= (m*adaptBuffer - (fullEnoughAdapt*m));
                fullEnoughSplitCheck = false(1, n_nodes);
                % fullEnoughSplitCheck(splitWait > 0) = false;
                fullEnoughAdaptCheck(w_nodes(:, jminus1) < 1e-10) = false;
                % fullEnoughSplitCheck(splitWait <= 0) = true;
                fullEnoughSplitCheck = splitWait <= 0;
            end
            num_nodesPlot(j) = n_nodes;
            % Calculating output values
            [riPlus(:, j), Pri(:, :, j)] = outputGMM(riPlus_nodes(:, :, j), Pri_nodes(:, :,:,  j), w_nodes(:, j), n_nodes);
            [x(:, j), P(:, :, j)] = outputGMM(x_nodes(:, :, j), P_nodes(:, :, :, j), w_nodes(:, j), n_nodes);

            if contains(flagQ, 'CorrMethod')
                % Should not need to shift buffer back by 1 because it was
                % already shifted back for this time step prior to the
                % measurment update 

                % Should only need to overwirte the most recent window
                % entry  
                rWindow_nodes(1:m, :) = riMinus_nodes(:, :, j);
                Fwindow_nodes(:, :, :, 1) = Phi(targetStates, targetStates, :, j); 
                Hwindow_nodes(1:m, :, :) = H_hist(1:m, targetStates, :, j);


                % update full enough criteria
                NaNcheck = sum(isnan(rWindow_nodes));
                fullEnoughAdaptCheck = NaNcheck <= (m*adaptBuffer - (fullEnoughAdapt*m));
                % fullEnoughSplitCheck = false(1, n_nodes);
                % fullEnoughSplitCheck(splitWait > 0) = false;
                fullEnoughAdaptCheck(w_nodes(:, jminus1) <1e-10) = false;
                % fullEnoughSplitCheck(splitWait <= 0) = true;
                fullEnoughSplitCheck = splitWait <= 0;


                if j > 100
                    flagIncMean = true;
                end

                if j > adaptBuffer+1
                    if j < adaptBuffer+2
                        flagBufferReqMet = false(1, n_nodes);
                    elseif mod(j, adaptBuffer+1) == 0
                        flagBufferReqMet = true(1, n_nodes);
                        % Check to see if each buffer is full
                        % NaNcheck = sum(isnan(rWindow_nodes));
                        % if buffer is not full, don't update Q 
                        % -- will happen when splitting on times that 
                        % don't line up with refill
                        flagBufferReqMet(~fullEnoughAdaptCheck) = false; 
                        if j < 300
                            % fprintf('update on time step: %d\n', j)
                            % mean_buff = mean(reshape(rWindow_nodes, m, adaptBuffer+1, []), 2, 'omitnan');
                            % norm_mean_buff = reshape(pagenorm(mean_buff), 1, []);
                            % fprintf('norm rWindow %0.2f, w/ split value %0.1f\n', norm_mean_buff, epsilon_split)
                        end
                    else  
                        flagBufferReqMet = false(1, n_nodes);
                        % This will save off converted info to cov form
                        qkplus = qkminus;
                        Pqkplus = Pqkminus;
                        [qCountDiag(:, j), qCountDiagCov(:, :, j)] = outputGMM(qkplus, Pqkplus, w_nodes(:, j), n_nodes); % is thig legal?
                    end


                    [qkplus, Pqkplus, Awindow_nodes, Bwindow_nodes, Zwindow_nodes, vBarHist(:, :, j)] = adaptiveProcNoise_estPSD_page(qkminus, Pqkminus, Xk_nodes, R, rWindow_nodes, uWindow_nodes, Fwindow_nodes, Hwindow_nodes, Awindow_nodes, Bwindow_nodes, Zwindow_nodes, adaptBuffer, flagBufferReqMet, QW, RW, beta(1), dt, flagQ, flagIncMean);
                    % [~, ~, AwindowFull, BwindowFull, ZwindowFull] = adaptiveProcNoise_estPSD(qkminus, Pqkminus, Xk, R, rWindow, uWindow, FwindowFull, HwindowFull, AwindowFull, BwindowFull, ZwindowFull, adaptBuffer, flagBufferReqMet, QW, RW, beta(1), dt, flagQ);


                    % Updates "minus step" for either formulation to
                    % its consistent represtnation
                    qkminus = qkplus;
                    Pqkminus = Pqkplus;

                    if contains(flagQ, 'Inf')
                        % Converts information out of information form
                        % to cov form to apply to larger EKF sim
                        dkplus = qkplus;
                        Lambdakplus = Pqkplus;
                        Pqkplus = pagemldivide(Lambdakplus, eye(3));
                        qkplus = reshape(pagemtimes(Pqkplus, reshape(dkplus, [], 1, n_nodes)), [], n_nodes);
                    end
                    % This will save off converted info to cov form
                    [qCountDiag(:, j), qCountDiagCov(:, :, j)] = outputGMM(qkplus, Pqkplus, w_nodes(:, j), n_nodes); % is thig legal?
                    % qCountDiagCov(:, :, j) = Pqkplus;
                    % qCountDiag(:, j) = qkplus;

                    if mod(j, 1000) == 0
                        % Prints updates every 1000 steps
                        % fprintf('qest = %e, %e, %e\n', qkplus)
                    end

                    % Updates discrete time proces noise used in larger
                    % EKF
                    QhatUT = pagemtimes(Xk_nodes, (reshape(qkplus.^2, [], 1, n_nodes)));
                    QhatTemp = reverseVech_page(QhatUT);
                    Qhat(targetStates, targetStates, :) = QhatTemp;
                    % if norm(QhatTemp) > 10
                    %     fprintf('Something is terribly wrong, norm(QhatTemp)>10\n')
                    % end
                end
            elseif contains(flagQ, 'CovMatchTrad') % possibly implement covariance matching here... 
                % 1/N_window * sum(K*r*r'*K')
            end
            if mod(j, 1000) == 0
                fprintf('Time step %d\n', j)
            end
        end 
        jminus1 = j;
    end
    %% Save NEES
    warning('off','all')
    for j = 1:N
        NEES(i, j) = (truthState(:, j) - x(:, j))'*(P(:, :, j)^-1)*(truthState(:, j) - x(:, j));
        NEESdiag(i, j) = (truthState(:, j) - x(:, j))'*(diag(diag(P(:, :, j)))^-1)*(truthState(:, j) - x(:, j));
        NEESfirst6(i, j) = (truthState(targetStates, j) - x(targetStates, j))'*(P(targetStates, targetStates, j)^-1)*(truthState(targetStates, j) - x(targetStates, j));
        NEESfirst6diag(i, j) = (truthState(targetStates, j) - x(targetStates, j))'*(diag(diag(P(targetStates, targetStates, j)))^-1)*(truthState(targetStates, j) - x(targetStates, j));
    end
    warning('on','all')
    %% Save Data for Plotting
    if nMC > 1
        MCcounter = i;
    % elseif exist('iiP0', 'var')
    %     % Will set to sstudy counter over truth 
    %     MCcounter = iiP0;
    elseif exist('iiTruth', 'var')
        MCcounter = iiTruth;
    else
        MCcounter = i;
    end
    fileNameSave = fullfile(fileBase, strcat(accelName, flagGravModelTruth, flagGravModel, burnName), strcat('MC_', num2str(MCcounter)));
    if ~exist(fileNameSave, 'dir')
        mkdir(fileNameSave)
    end
    save(fullfile(fileNameSave, 'xhist.mat'), 'x')
    save(fullfile(fileNameSave, 'PPlushist.mat'), 'P')
    save(fullfile(fileNameSave, 'Pminushist.mat'), 'Pminus')
    save(fullfile(fileNameSave, 'riPlushist.mat'), 'riPlus')
    save(fullfile(fileNameSave, 'riMinushist.mat'), 'riMinus')
    save(fullfile(fileNameSave, 'PriPlushist.mat'), 'Pri')
    save(fullfile(fileNameSave, 'PriMinushist.mat'), 'PriMinus')
    save(fullfile(fileNameSave, 'meas_runhist.mat'), 'meas_run')
    save(fullfile(fileNameSave, 'qCountDiag_runhist.mat'), 'qCountDiag')
    save(fullfile(fileNameSave, 'qCov.mat'), 'qCountDiagCov')
    save(fullfile(fileNameSave, 'numNodes.mat'), 'num_nodesPlot')
end

save(fullfile(fileBase, strcat(accelName, flagGravModelTruth, flagGravModel, burnName), 'NEEShist.mat'), 'NEES')
save(fullfile(fileBase, strcat(accelName, flagGravModelTruth, flagGravModel, burnName), 'NEESdiaghist.mat'), 'NEESdiag')

save(fullfile(fileBase, strcat(accelName, flagGravModelTruth, flagGravModel, burnName), 'NEESTarghist.mat'), 'NEESfirst6')
save(fullfile(fileBase, strcat(accelName, flagGravModelTruth, flagGravModel, burnName), 'NEESTargdiaghist.mat'), 'NEESfirst6diag')

% fprintf("Couldn't update %d times out of %d\n", noUpdateCount, N)
% 
% if nMC == 1
%     figure(112)
%     plot3(x(1, :), x(2, :), x(3, :))
%     hold on
%     grid minor
%     plot3(truthState(1, :), truthState(2, :), truthState(3, :))
%     xlabel('X [m]')
%     ylabel('Y [m]')
%     zlabel('Z [m]')
%     if contains(flagScenario, 'Rel')
%         title('Target Relative to Chaser Est and Truth in Inertial Space')
%     else
%         title('Target Est and Truth in Inertial Space')
%     end
% 
%     figure(113)
%     plot3(x(chaserStates(1), :), x(chaserStates(2), :), x(chaserStates(3), :))
%     hold on
%     grid minor
%     plot3(truthState(chaserStates(1), :), truthState(chaserStates(2), :), truthState(chaserStates(3), :))
%     xlabel('X [m]')
%     ylabel('Y [m]')
%     zlabel('Z [m]')
%     title('Chaser Est and Truth in Inertial Space')
% end

% figure(375)
% hold on
% grid minor
% plot(time, vBarHist(1, :), '.')
% plot(time, vBarHist(2, :), '.')
% plot(time, vBarHist(3, :), '.')


% figure(375)
% hold on
% plot(time, num_nodesPlot, 'o')
% xlabel('Time')
% ylabel('Nodes []')
% 
% figure(376)
% hold on
% plot(time(100:300), num_nodesPlot(100:300), 'o')
% xlabel('Time')
% ylabel('Nodes []')
% % 
% cndw = 0;
