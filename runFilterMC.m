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
n_exp = length(expNum);
Q0 = zeros(3*n_exp);
for i = 1:n_exp
    if expNum(i) == 0
        Q0(((i*3-2):i*3), ((i*3-2):i*3)) = zeros(3);
        procNoiseName = strcat(procNoiseName, '_zeroQ');
    else
        Q0(((i*3-2):i*3), ((i*3-2):i*3)) = eye(3)*10^(-expNum(i)); % 1e-6 is good for no accel -- 1e-12 is good for accel on
        procNoiseName = strcat('_nonzeroQ_eneg', num2str(expNum(i)));
    end
end




[measNum, ~] = size(meas);
if measNum == 5
    measFunc = @relNavRngRngRate_opTarget_feature_meas_H;
    R_rangeRangeRate_chol = sqrt(R(1, 1));
    R_opNav_chol = chol(R(2:3, 2:3), 'lower');
    measFlag = 'relRngRelAzElFeatAzEl';
elseif measNum == 3
    if contains(flagScenario, 'direct')
        measFunc = @relNavDirectInertial_meas_H;
        R_direct_chol = chol(R, 'lower');
        measFlag = 'direct';
    else
        measFunc = @relNavRng_opTarget_meas_H;
        R_rangeRangeRate_chol = sqrt(R(1, 1));
        R_opNav_chol = chol(R(2:3, 2:3), 'lower');
        measFlag = 'relRngRelAzEl';
        if ~contains(flagScenario, 'Rel')
            measFunc = @rng_opTarget_meas_H;
        end
    end
else
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
    deputyStates = 7:12;
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
    deputyStates = 7:15;
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
    deputyStates = 10:15;
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
    deputyStates = 10:18;
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
% TODO: Remove this next line
 % P0 = P0*0;
N = length(truthState);
NEES = NaN(nMC, N);
NEESdiag = NaN(nMC, N);
NEESfirst6 = NaN(nMC, N);
NEESfirst6diag = NaN(nMC, N);
[m, ~] = size(R);
Nwindow = 30;

%% Run EKF
for i = 1:nMC
    %% Initialize Variables
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
    vBarHist = zeros(3, N);
    x = NaN(n, N);
    P = NaN(n, n, N);
    Phi = NaN(n, n, N);
    Phi_window = NaN(n, n, Nwindow);
    Pminus = NaN(n, n, N);
    Pri = NaN(m, m, N);
    PriMinus = NaN(m, m, N);
    H_hist = NaN(m, n, N);
    x(:, 1) = x0;
    P(:, :, 1) = P0;
    Pminus(:, :, 1) = P0;
    Phi(:, :, 1) = eye(n);
    qCountDiag = NaN(3, N);
    qCountDiagCov = NaN(3, 3, N);
    iminus1 = 1;
    % TODO: Take the *0 out 
    accelerometerMeas = burnHist;
    riMinus = NaN(m, N);
    riPlus = NaN(m, N);
    Qiminus1tilde = Q0;
    dt = time(2) - time(1);
    flagIncMean = false;
    if dt > 1
        dtProp = 1;
    else
        dtProp = dt;
    end
    
    Qhat_old = continuious2discreteProcNoise(dtProp, beta, Qiminus1tilde, flagEstTargetAccel, flagEstChaserAccelBias, relNavFlag, flagNoiseonJerk );
    Qhat = zeros(n);
    if n_elm(1) > 6
        Xk = Xk_relInertOnly_jerkNoise(time(2) - time(1));
        if contains(flagScenario, 'trying')
            Xk = Xk_relInertOnly_accelNoise(time(2) - time(1));
        end
    else
        [Xk] = Xk_relInertOnly_ZOH(time(2) - time(1));
    end
    QhatUT = Xk*(diag(Qiminus1tilde(1:3, 1:3)));
    QhatTemp = reverseVech(QhatUT);
    Qhat(targetStates, targetStates) = QhatTemp;

    Xk_chaser = Xk_relInertOnly_ZOH(time(2) - time(1));
    QhatUT = Xk_chaser*(diag(Qiminus1tilde(4:6, 4:6)));
    QhatTemp = reverseVech(QhatUT);
    QhatTemp = blkdiag(QhatTemp, zeros(3));
    Qhat(deputyStates, deputyStates) = QhatTemp;

    % TODO: Take this out
    % [Xk] = Xk_relInertOnly_ZOH(time(2) - time(1));
    % QhatUT = 4*Xk*(qkminus0.^2);
    % QhatTemp = reverseVech(QhatUT);
    % Qhat(targetStates, targetStates) = QhatTemp;
    

    if length(Qhat) < n
        Qhat = blkdiag(Qhat, 0);
    end
    %% Process Noise Stuff 
    if strcmp(flagQ, 'CovMatchTrad')
        Q_hist = zeros(n, n, Nwindow);
    elseif strcmp(flagQ, 'CovMatch')
        [Xk] = Xk_ADMC(beta(1), time(2) - time(1));
        [n_w, ~] = size(diag(vech(eye(n))));
        W_hist = zeros(n_w, n_w, Nwindow);
    elseif contains(flagQ, 'CorrMethod')
        if n_elm(1) > 6
            Xk = Xk_relInertOnly_jerkNoise(time(2) - time(1));
            if contains(flagScenario, 'trying')
                Xk = Xk_relInertOnly_accelNoise(time(2) - time(1));
            end
        else
            [Xk] = Xk_relInertOnly_ZOH(time(2) - time(1));
        end
        Pqkminus = Pqkminus0;
        % TODO: Change this back later, making a change to keep the proc
        % noise overestimated to start -- change rand to randn to change
        % back
        deviationFromqkminus = chol(Pqkminus)*randn(3, 1);
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
        
        rWindow = NaN(m, adaptBuffer+1);
        uWindow = zeros(3, adaptBuffer +1);
        nA = length(targetStates);
        Fwindow = zeros(nA, nA, adaptBuffer);
        Hwindow = zeros(m, nA, adaptBuffer+1);
        % % TODO: Take this out
        % FwindowFull = zeros(12, 12, adaptBuffer);
        % HwindowFull = zeros(m, 12, adaptBuffer+1);
        % AwindowFull = zeros(12, m, adaptBuffer+1, adaptBuffer+1);
        % BwindowFull = zeros(12, 12, adaptBuffer, adaptBuffer);
        % ZwindowFull = zeros(12, adaptBuffer+1);
        % GWindowFull = [];
        Awindow = zeros(nA, m, adaptBuffer+1, adaptBuffer+1);
        Bwindow = zeros(nA, nA, adaptBuffer, adaptBuffer);
        Zwindow = zeros(nA, adaptBuffer+1);
        GWindow = [];
        noUpdateCount = 0;
        xnom = x(:, 1);
    end
    %%
    meas_run = NaN(m, N);
    
    jminus1 = 1;
    accelBias = zeros(3, 1);
    fprintf('Running EKF MC run %d...\n', iiTruth)
    for j = 2:N
        %% Propagation 
        tj = time(j);
        times = [time(jminus1), tj];
        if ~flagAccelerometerPerfectMatch
            accelerometerMeas(:, j) = accelerometerMeas(:, j) +  accelChol*randn(3, 1)*sqrt(dtProp);
        end
        if norm(accelerometerMeas(:, j)) > gate
            fprintf('Burn Detected, time step %d:%d\n', jminus1:j)
           if ~flagEstChaserAccelBias
               accelBias = mean(accelerometerMeas(:, (j-121):(j-1)), 2); % body -- averaging 2000 non maneuverg points to get bias
           end
        end
        Phi(:, :, j) = eye(n);
        if j ~=1
            timeProp = times(1):dtProp:times(2);
            Nprop = length(timeProp)-1;
            state = x(:, jminus1);
            P(:, :, j) = P(:, :, jminus1);
            Rbody_N = RBody2_N(:, :, j);
            for k = 1:Nprop
                initialState = state;
                if strcmp(flagFilterType, 'GSKF')
                    [state, P(:, :, j), Phi_k] = GSKFprop(timeProp(k:k+1), initialState, P(:, :, j), eom, Qhat, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, RN_rot, relNavFlag, accelBias, n_elm, accelerometerVRW*dt);
                elseif contains(flagScenario, 'NRHO')
                    [state, P(:, :, j), Phi_k] = EKFprop(timeProp(k:k+1), initialState, P(:, :, j), eom, Qhat, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, vBarHist(:, j-1), relNavFlag, zeros(3, 1), n_elm, accelerometerVRW*dt, flagUseTrueJacobians, squeeze(Phi_true(:, :, j)), CRTBPFlag);
                    %                            EKFprop(times,           x,            P,          eom, Qhat, gravParam, beta_markov, accelerometer,    gate, Rbody_N, RN_rot,           relNavFlag, accelBias, n_elm, accelVRW,            flagUseTrueJacobians, Phi_true,                   CRTBPFlag)
                else                               
                    [state, P(:, :, j), Phi_k] = EKFprop(timeProp(k:k+1), initialState, P(:, :, j), eom, Qhat, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, vBarHist(:, j-1), relNavFlag, accelBias, n_elm, accelerometerVRW*dt, flagUseTrueJacobians, squeeze(Phi_true(:, :, j)));
                end
                Phi(:, :, j) = Phi_k*Phi(:, :, j);
            end
            x(:, j) = state;
            xminus = state;
            Pminus(:, :, j) = P(:, :, j);
        end
        if flagRunMeasUpdate
            %% Measurment Step
            % Add noise to measurments
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
            if strcmp(flagFilterType, 'EKF')
                if contains(flagScenario, 'NRHO') % CRTBPFlag
                    [x(:, j), P(:, :, j), riPlus(:, j), riMinus(:, j), Pri(:, :, j), PriMinus(:, :, j), H] = EKFmeas(x(:, j), P(:, :, j), meas_j, R, measFunc, L_j, RCF1_Body2_N(:, :, j), RCF2_Body2_N(:, :, j), vBarHist(:, j-1), considerParam, n_elm, measFlag, underweight, accelerometerMeas(:, j), gate, accelerometerVRW*dt, flagUnderweight, flagUseTrueJacobians, squeeze(H_true(:, :, j)), CRTBPFlag, L);
                else
                    [x(:, j), P(:, :, j), riPlus(:, j), riMinus(:, j), Pri(:, :, j), PriMinus(:, :, j), H] = EKFmeas(x(:, j), P(:, :, j), meas_j, R, measFunc, L_j, RCF1_Body2_N(:, :, j), RCF2_Body2_N(:, :, j), vBarHist(:, j-1), considerParam, n_elm, measFlag, underweight, accelerometerMeas(:, j), gate, accelerometerVRW*dt, flagUnderweight, flagUseTrueJacobians, squeeze(H_true(:, :, j)));
                end
              % [x,       P,          riPlus,       ri, Pri, PriMinus, H]
            else % default to 2nd order meas update 
                [x(:, j), P(:, :, j), riPlus(:, j), riMinus(:, j), Pri(:, :, j), PriMinus(:, :, j), H] = GSKFmeas(x(:, j), P(:, :, j), meas_j, R, measFunc, L_j, RCF1_Body2_N(:, :, j), RCF2_Body2_N(:, :, j), vBarHist(:, j-1), considerParam, n_elm, measFlag, underweight, accelerometerMeas(:, j), gate, accelerometerVRW*dt);
            end
            H_hist(:, :, j) = H(1:m, :);

            ObsMat = [H(1:m, targetStates); H(1:m, targetStates)*Phi(targetStates, targetStates, j); H(1:m, targetStates)*(Phi(targetStates, targetStates, j)*Phi(targetStates, targetStates, j))];
            % testObs = rank(ObsMat);
            if contains(flagQ, 'CorrMethod')
                % if relNavFlag
                %     % If Conversion needed
                %     xnomr = xnom(1:n_elm(1));
                %     xnomc = xnom(n_elm(1) + 1:end);
                %     xnomt = xnomr;
                %     xnomt(1:6) = xnomc(1:6) + xnomr(1:6);
                %     initialStatenom = [xnomt; xnomc];
                % else
                %     initialStatenom = xnom;
                % end
                % % statenom = rk4(initialStatenom, times, eom, gravParam, beta, accelerometerMeas(:, j), gate, Rbody_N, RN_rot, n_elm, accelBias, {});

                % if relNavFlag
                %     xnomt = statenom(1:n_elm(1));
                %     xnomc = statenom(n_elm(1)+1:end);
                %     xnomr = xnomt;
                %     xnomr(1:6) = xnomt(1:6) - xnomc(1:6);
                % 
                %     xnom = [xnomr; xnomc];
                % else
                %     xnom = statenom;
                % end
                % [Fnom] = buildPhi(initialStatenom, gravParam, dt, beta, n_elm,relNavFlag);
                % [ynom, Hnom] = measFunc(xminus, RCF1_Body2_N(:, :, j), n_elm);
%                 [ynom, Hnom] = measFunc(xnom, RCF1_Body2_N(:, :, j), n_elm);
                % + H(1:m, targetStates)*(x(targetStates, j) -xminus(targetStates))
                rWindow = [riMinus(:, j), rWindow(:, 1:end-1)]; %[riMinus(1:m, j), rWindow(:, 1:end-1)]; %[meas_j - meas_truth(:, j), rWindow(:, 1:end-1)]; %[riMinus(1:m, j), rWindow(:, 1:end-1)]; % [riMinus(1:m, j) - H(1:m, targetStates)*(x(targetStates, j) -xminus(targetStates)), rWindow(:, 1:end-1)]; % + H(1:m, targetStates)*(x(targetStates, j) -xminus(targetStates)) [riMinus(1:m, j), rWindow(:, 1:end-1)]; % % + H(1:m, targetStates)*(x(targetStates, j) -xminus(targetStates)); % [meas_j - meas_truth(:, j), rWindow(:, 1:end-1)]; %[riMinus(1:m, j), rWindow(:, 1:end-1)]; % [riMinus(1:m, j) - H(1:m, targetStates)*(x(targetStates, j) -xminus(targetStates)), rWindow(:, 1:end-1)]; % + H(1:m, targetStates)*(x(targetStates, j) -xminus(targetStates)) [riMinus(1:m, j), rWindow(:, 1:end-1)]; %
                Fwindow(:, :, 2:end) = Fwindow(:, :, 1:end-1);
                Fwindow(:, :, 1) = Phi(targetStates, targetStates, j); %Phi_true(targetStates, targetStates, j); %  
                Hwindow(:, :, 2:end) = Hwindow(:, :, 1:end-1);
                Hwindow(:, :, 1) = H_hist(1:m, targetStates, j); %  H_true(1:m, targetStates, j); %
                % FwindowFull(:, :, 2:end) = FwindowFull(:, :, 1:end-1);
                % FwindowFull(:, :, 1) = Phi(1:end-3, 1:end-3, j); %Phi_true(targetStates, targetStates, j); %  
                % HwindowFull(:, :, 2:end) = HwindowFull(:, :, 1:end-1);
                % HwindowFull(:, :, 1) = H_hist(1:m, 1:end-3, j); %  H_true(1:m, targetStates, j); %
                if j == 118
                    hvios = 0;
                elseif j == 130
                    vdhjos = 0;
                end
                if j > 100
                    flagIncMean = true;
                end

                if j > adaptBuffer+1
                    if j < adaptBuffer+2
                        flagBufferReqMet = false;
                    elseif mod(j, adaptBuffer+1) == 0
                        flagBufferReqMet = true;
                        if j < 300
                            fprintf('update on time step: %d\n', j)
                        end
                    else
                        flagBufferReqMet = false;
                        % if contains(flagQ, 'Inf')
                        %     % Converts information out of information form
                        %     % to cov form to apply to larger EKF sim 
                        %     dkplus = qkplus;
                        %     Lambdakplus = Pqkplus;
                        %     Pqkplus = Lambdakplus\eye(3);
                        %     qkplus = Pqkplus*dkplus;
                        % end
                        % This will save off converted info to cov form
                        qkplus = qkminus;
                        Pqkplus = Pqkminus;
                        qCountDiagCov(:, :, j) = Pqkplus;
                        qCountDiag(:, j) = qkplus;
                        
                    end
                    if isempty(GWindow)
                        [qkplus, Pqkplus, Awindow, Bwindow, Zwindow, vBarHist(:, j)] = adaptiveProcNoise_estPSD(qkminus, Pqkminus, Xk, R, rWindow, uWindow, Fwindow, Hwindow, Awindow, Bwindow, Zwindow, adaptBuffer, flagBufferReqMet, QW, RW, beta(1), dt, flagQ, flagIncMean);
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
                            Pqkplus = Lambdakplus\eye(3);
                            qkplus = Pqkplus*dkplus;
                        end
                        % This will save off converted info to cov form
                        qCountDiagCov(:, :, j) = Pqkplus;
                        qCountDiag(:, j) = qkplus;

                        if mod(j, 1000) == 0
                            % Prints updates every 1000 steps
                            % fprintf('qest = %e, %e, %e\n', qkplus)
                        end
                        
                        % Updates discrete time proces noise used in larger
                        % EKF
                        QhatUT = Xk*(qkplus.^2);
                        QhatTemp = reverseVech(QhatUT);
                        Qhat(targetStates, targetStates) = QhatTemp;
                        % if norm(QhatTemp) > 10
                        %     fprintf('Something is terribly wrong, norm(QhatTemp)>10\n')
                        % end
                    else
                        [qkplus, Pqkplus, Awindow, Bwindow, Zwindow, vBarHist(:, j)] = adaptiveProcNoise_estPSD(qkminus, Pqkminus, Xk, R, rWindow, uWindow, Fwindow, Hwindow, Awindow, Bwindow, Zwindow, adaptBuffer, flagBufferReqMet, QW, RW, beta(1), dt, flagQ, GWindow);
                        % Updates "minus step" for either formulation to
                        % its consistent represtnation
                        qkminus = qkplus;
                        Pqkminus = Pqkplus;
                        
                        if contains(flagQ, 'Inf')
                            % Converts information out of information form
                            % to cov form to apply to larger EKF sim 
                            dkplus = qkplus;
                            Lambdakplus = Pqkplus;
                            Pqkplus = Lambdakplus\eye(3);
                            qkplus = Pqkplus*dkplus;
                        end
                        % This will save off converted info to cov form
                        qCountDiagCov(:, :, j) = Pqkplus;
                        qCountDiag(:, j) = qkplus;

                        if mod(j, 1000) == 0
                            % Prints updates every 1000 steps
                            % fprintf('qest = %e, %e, %e\n', qkplus)
                        end
                        
                        % Updates discrete time proces noise used in larger
                        % EKF
                        QhatUT = Xk*(qkplus.^2);
                        QhatTemp = reverseVech(QhatUT);
                        if norm(QhatTemp) < 10
                            Qhat(targetStates, targetStates) = QhatTemp;
                        end
                    end
                end
            elseif contains(flagQ, 'CovMatchTrad')
                % 1/N_window * sum(K*r*r'*K')
            end
            if mod(j, 5000) == 0
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
%     plot3(x(deputyStates(1), :), x(deputyStates(2), :), x(deputyStates(3), :))
%     hold on
%     grid minor
%     plot3(truthState(deputyStates(1), :), truthState(deputyStates(2), :), truthState(deputyStates(3), :))
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

cndw = 0;
