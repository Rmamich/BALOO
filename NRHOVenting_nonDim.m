%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Venting Scenario 
% 
%
% Created: 07/29/24
% Last Modified: 07/29/24
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%house keeping
if ~exist('flagSetFlags', 'var')
    clear
    close all
    flagScenario = 'NRHO_VentingRel_nonDim';
    fileBase = fullfile('.', 'MAT', flagScenario);
    measNum = 3;
    flagGravModelTruth = 'j2On';                      % Sets label for saving data with the truth data model
    plotCheck = true;
    beta = 0*[0, 1/3600, 0];                            % time constants for Markov Processes
    bias = 0*[(300*9.81e-6); (300*9.81e-6); (300*9.81e-6)]; % initial bias (300*9.81e-6) LN200S
    biasLidar = 0;
    gate = 1e-7*0;
    eom = @tbpJ2_accelerometerGating_rpod_lidarBias;
    flagTruthNoisyBurn = false;
    iiTruth = 1;
    expNum = [0, 0];  
    eomTruth = @CRTBP_targetChaser;
    flagNoiseonJerk = false;
end
rng(123 + iiTruth+1)
fprintf('Running Truth %d...\n', iiTruth)
%% Initialization
fileBaseBurns = fullfile('.', 'MAT', 'NRHOVentingScenario_nonDim');

load(fullfile(fileBaseBurns, 'truth_forwardBackward.mat'));
load(fullfile(fileBaseBurns, 'burnMagRendezvous.mat'))
load(fullfile(fileBaseBurns, 'burnDirRendezvous.mat'))
load(fullfile(fileBaseBurns, 'burnStartRendezvous.mat'))
load(fullfile(fileBaseBurns, 'burnEndRendezvous.mat'))
load(fullfile(fileBaseBurns, 'timeRendezvous.mat'))
load(fullfile(fileBaseBurns, 'rbody1_n.mat'))
load(fullfile(fileBaseBurns, 'rbody2_n.mat'))




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


xc = truthState_N_new(1:6, 1);
xd = truthState_N_new(10:15, 1);

dt = time(2)-time(1); 
N = length(time);

truthState_N = zeros(19, N);

Phi_true = zeros(19, 19, N);


burnMagSigma = 0.15;
burnDirSigma = 0.5; % deg

n_burn = [length(burnStart{1}), length(burnStart{2})];



burnHistCommand = zeros(6, N);
burnHistActual = zeros(6, N);
accelerometer = zeros(3, N);
posStates = [1:3; 10:12];
velStates = [4:6; 13:15];



calcBurn = false(2, 1);

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
% RBody1_N = NaN(3, 3, N);
% RBody2_N = NaN(3, 3, N);
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
% RBody1_N(:, :, 1) = cameraPoint(truthState_N(1:6, 1), truthState_N(10:15, 1), pointNoise);
% RBody2_N(:, :, 1) = cameraPoint(truthState_N(10:15, 1), truthState_N(1:6, 1), pointNoise);

burnCount = [1, 1];

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

Qaccelerometer = accelerometerVRW*dt;
Qaccelerometer_chol = chol(Qaccelerometer/100, 'lower');

burnVecCommand = zeros(6, 1);
burnVecActual = zeros(6, 1);


for i = 2:N
    % DCMs
    RN_rot(:, :, i) = rotz(lam0 + omegaEarth*time(i));
    % RBody1_N(:, :, i) = cameraPoint(truthState_N(1:6, iminus1), truthState_N(10:15, iminus1), pointNoise);
    % RBody2_N(:, :, i) = cameraPoint(truthState_N(10:15, iminus1), truthState_N(1:6, iminus1), pointNoise);
    % Rbody_N = RBody2_N(:, :, i);
    timeInt = time(iminus1:i);
    for k = 1:2
        if burnCount(k) <= n_burn(k)
            burnStart_temp = burnStart{k};
            burnStart_temp = burnStart_temp(burnCount(k));
            burnEnd_temp = burnEnd{k};
            burnEnd_temp = burnEnd_temp(burnCount(k));
            burn_mag = burn{k};
            burn_mag = burn_mag(burnCount(k));
            burn_dir = burnDirection_body{k};
            burn_dir = burn_dir(:, burnCount(k));
            if timeInt(2) >= burnStart_temp && timeInt(1) < burnEnd_temp
                % timeInt(1) < burnTimeEndc && timeInt(2)>= burnTimeInitc
                burnUnitNoisy = -rotz(randn(1)*burnDirSigma)*roty(randn(1)*burnDirSigma)*rotx(randn(1)*burnDirSigma)*burn_dir;
                if k == 1
                    % burnVecCommand(((k-1)*3+1):k*3) = burn_mag*(RBody1_N(:, :, i)'*burn_dir); % target burns represented and saved in inertial 
                    burnVecCommand(((k-1)*3+1):k*3) = burn_mag*(burn_dir);
                else
                    burnVecCommand(((k-1)*3+1):k*3) = burn_mag*RBody2_N(:, :, i)*burn_dir; % chaser burns represented and saved in body frame
                    % burnVecCommand(((k-1)*3+1):k*3) = burn_mag*burn_dir; % chaser burns represented and saved in body frame
                    % fprintf('%d burning %d, %d\n', k, iminus1:i)
                end
                burnVecActual(((k-1)*3+1):k*3) = burnVecCommand(((k-1)*3+1):k*3);
                calcBurn(1) = true;
                % fprintf('%d burning %d, %d\n', k, timeInt/dt)
            elseif timeInt(1) >= burnEnd_temp
                burnCount(k) = burnCount(k) + 1;
                calcBurn(k) = false;
                burnVecCommand(((k-1)*3+1):k*3) = zeros(3, 1);
                burnVecActual(((k-1)*3+1):k*3) = zeros(3, 1);
            end
        end
    end
    burnHistActual(:, i) = burnVecActual;
    burnHistCommand(:, i) = burnVecCommand;
    timeProp = timeInt(1):dtProp:timeInt(2);
    Nprop = length(timeProp)-1;
    initialState = truthState_N(:, iminus1);
    initialState(7:9) = burnVecActual(1:3);
    %TODO: uncomment noise
    accelerometer(:, i) = burnVecActual(4:6) + initialState(16:18);% + Qchol_c*randn(3, 1);
    Rbody_N = RBody2_N(:, :, i);
    Phi_true(:, :, i) = eye(19);
    for j = 1:Nprop
        state = rk4(initialState, timeProp(j:j+1), eomTruth, gravParamTruth, beta, accelerometer(:, i), gate, Rbody_N, RN_rot, [9, 9], {});
        [F_i] = buildPhi(initialState, gravParamTruth, dt, beta, [9, 9], true, true, timeProp(j));
        Phi_true(:, :, i) = F_i*Phi_true(:, :, i);
        initialState = state;
    end
    state(16:18) = state(16:18) + Qaccelerometer_chol*randn(3, 1);
    truthState_N(:, i) = state; % rot and nondim 
    % % truthState_R_new(:, i) = initialState_nonDimR;
    % truthState_N(:, i) = nonDimRotational_dimInertial_conversion(state, timeProp(end), L_m, n, tCharacteristic, 'fullState_nondimR2dimN');
    % RBody1_N(:, :, i) = cameraPoint(truthState_N(1:6, i), truthState_N(10:15, i), pointNoise);
    % RBody2_N(:, :, i) = cameraPoint(truthState_N(10:15, i), truthState_N(1:6, i), pointNoise);
    iminus1 = i;
    if mod(i, 10000) == 0
       % fprintf('%d\n', i) 
    end
end

truthState = truthState_N; % rot and non dim dual 
truthState(1:6, :) = truthState_N(1:6, :) - truthState_N(10:15, :); % rot and non dim relative 
truthState_dual = truthState_N;



%% Verify Traj
if plotCheck
    figure(201)
    plot3(truthState_N(1, :), truthState_N(2, :), truthState_N(3, :))
    hold on
    grid minor
    plot3(truthState_N(10, :), truthState_N(11, :), truthState_N(12, :))

    figure(203)
    for k = 1:3
        subplot(3, 1, k)
        hold on
        grid minor
        plot(time, truthState_N(k, :))
        plot(time, truthState_N(k+9, :))
        title('both plotted')
    end

    figure(204)
    for k = 1:3
        subplot(3, 1, k)
        hold on
        grid minor
        plot(time, truthState_N(k, :)-truthState_N(k+9, :))
        title('target - chaser')
    end
    outputFlag = "kep";
    inputFlag = "cartesian";

    [ oec ] = orbitalElmConverter( inputFlag, truthState_N(1:6, end), 0, outputFlag, mu );
    [ oed ] = orbitalElmConverter( inputFlag, truthState_N(10:15, end), 0, outputFlag, mu );


    figure(205)
    plot3(truthState_N(1, :), truthState_N(2, :), truthState_N(3, :))
    hold on
    plot3(truthState_N(10, :), truthState_N(11, :), truthState_N(12, :))
    % plot3(dimStateTimeHist(1, :), dimStateTimeHist(2, :), dimStateTimeHist(3, :))
    % % hold on
    % plot3(dimStateTimeHist(10, :), dimStateTimeHist(11, :), dimStateTimeHist(12, :))
    title(strcat('Plot of NRHO in the Inertial Frame'))
    xlabel('Inertial X [km]')
    ylabel('Inertial Y [km]')
    zlabel('Inertial Z [km]')
    plot3(0, 0, 0, '.', 'MarkerSize', 10)
    axis equal
    grid minor
    % legend('complexProp Target', 'complexProp Chaser', 'simpleProp Target', 'simpleProp Chaser', 'Earth')รท
    legend('complexProp Target', 'complexProp Chaser', 'Earth')
    
    
    figure(206)
    for k = 1:3
        subplot(3, 1, k)
        hold on
        grid minor
        plot(time, truthState_N(k+6, :))
    end
    title('Acceleration')

end

%% Generate Measuremnts
% focalLength = 7.6e-3; % m
% fov_rad = [44; 32]*pi/180; % rad
% numPixel = [2592; 1944]; % -
% 
% R_opNav = fov_rad./numPixel/2;
% From Jena Optronik RVS 3000-X Lidar -- GNSki 2024
% 128 x 128 pixles, 10 deg fov
% fov_rad = [10; 10]*pi/180; % rad
% numPixel = [128; 128]; % -

% R_opNav = fov_rad./numPixel/2;

% R_opNav = ones(2, 1)*(0.01*pi/180); % rad^2 -- 

% % From ASC Flash LiDAR -- IEEE paper --amzajerdian2013lidar
% fov_rad = [5; 5]*pi/180; % rad
% numPixel = [128; 128]; % -
% % R_opNav = fov_rad./numPixel;
% 
% R_opNav = ones(2, 1)*(0.1*pi/180);
% 
% sigma_r = 0.03/L_m; % [-] 
% sigma_rDot = 0.001/L_m*tCharacteristic; % [-]

% From ASC Flash LiDAR -- IEEE paper --amzajerdian2013lidar

% This set of values works and I know it-- but I want to try what is
% uncommented to see how things work when it matches the LEO scenario 
% fov_rad = [5; 5]*pi/180; % rad
% numPixel = [128; 128]; % -
% R_opNav = [0.04; 0.04]*pi/180;%fov_rad./numPixel;
% 
% sigma_r = 0.08/L_m; % [-] 
% sigma_rDot = 0.001/L_m*tCharacteristic; % [-]



fov_rad = [5; 5]*pi/180; % rad
numPixel = [128; 128]; % -

R_opNav = fov_rad./numPixel;

sigma_r = 0.01/L_m; % [-] 
sigma_rDot = 0.001/L_m*tCharacteristic; % [-]



% % From ASC Flash LiDAR -- IEEE paper --amzajerdian2013lidar
% fov_rad = [5; 5]*pi/180; % rad
% numPixel = [128; 128]; % -
% 
% R_opNav = fov_rad./numPixel;
% 
% sigma_r = 0.01/L_m; % [-] 
% sigma_rDot = 0.001/L_m*tCharacteristic; % [-]

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

        R = T*R*T'/10000; 
        R_direct_chol = chol(R, 'lower');
        
    else
        % range, relative bearing
        R = diag([sigma_r, R_opNav'].^2);
        R_rangeRangeRate_chol = sqrt(R(1, 1));
        R_rangeRangeRate = R(1, 1);
        R_opNav_chol = chol(R(2:3, 2:3), 'lower');
    end
elseif measNum == 5 % relativ erange and bearing and feature tracking 
    % range, relative bearing 
    R = diag([sigma_r, R_opNav', R_opNav'].^2);
    R_rangeRangeRate_chol = sqrt(R(1, 1));
    R_rangeRangeRate = R(1, 1);
    R_opNav_chol = chol(R(2:3, 2:3), 'lower');
end

meas = NaN(measNum, N);
meas_truth = NaN(measNum, N);

R_far = R;
R_near = diag([0.01/L_m, R_opNav'].^2);

R_all = zeros(3, 3, N);
% R = R*100000000;

% Assume that ACAF and ACI are coaligned at the start time
% Assume that ACAF rotates as a constant rate of omegaErros about the z
% axis
% Assume both s/c are keeping boresight approximately towards the other s/c
% groundFeature(x, RCF_N, RN_rot, L_m)
N_feature = 1;
RMoon = 1737400 / L;  % m/m [-]
% [L_set] = genFeatures(REarth, N_feature);
L_set = [1 - mu;0;0]; % if we pretend that the feature is at the center of the moon


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
        [meas_i, H_i] = relNavRngRngRate_opTarget_feature_meas_H(truthState(:, i), RCF1_Body2_N(:, :, i), RCF2_Body2_N(:, :, i), eye(3), L_set(:, 1), [9, 9]);
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
                % pass the meas model a new DCM
                NR = nonDimRotational_dimInertial_conversion(truthState(:, i), time(i), L_m, n, tCharacteristic, 'DCM_R2N');
                RCF1_Body2_R = RCF1_Body2_N(:, :, i);%*NR;
                [meas_i, H_i] = relNavRng_opTarget_meas_H(truthState(:, i), RCF1_Body2_R, [9, 9]);
                meas_truth(:, i) = meas_i;
                meas(1, i) = meas_i(1) + R_rangeRangeRate_chol*randn(1, 1);
                meas(2:3, i) = meas_i(2:3) + R_opNav_chol*randn(2, 1);
                H(:, :, i) = H_i;
                % if  j < 400 % meas(1, i) > 100/L_m &&
                %     R_all(:, :, i) = R_far;
                % else
                %     R_all(:, :, i) = R_near;
                % end
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

% R = R_all;
% figure(317)
% for iiPlotMeas = 1:3
% subplot(3, 1, iiPlotMeas)
% if iiPlotMeas > 1
%     plot(time, meas_truth(iiPlotMeas, :)*180/pi)
% else
%     plot(time, meas_truth(iiPlotMeas, :))
% end
% end

% figure(318)
% for iiPlotMeas = 1:2
%     subplot(2, 1, iiPlotMeas)
%     plot(time, meas_truth(iiPlotMeas+3, :)*180/pi)
% end


% meas = meas_truth;

%% Save Off Info
if ~exist(fileBase, 'dir')
    mkdir(fileBase)
end

if ~contains(flagScenario, 'Rel')
    truthState = truthState_N;
end



save(fullfile(fileBase, strcat('truthState_dual.mat')), 'truthState_dual')
save(fullfile(fileBase, strcat('truthState.mat')), 'truthState')
save(fullfile(fileBase, 'meas_.mat'), 'meas')
save(fullfile(fileBase, 'truemeas.mat'), 'meas_truth')
save(fullfile(fileBase, 'groundFeatIndex.mat'), 'featureIndex')
save(fullfile(fileBase, 'acaf2aci.mat'), 'RN_rot')
save(fullfile(fileBase, 'aci2body1CF1.mat'), 'RCF1_Body1_N')
save(fullfile(fileBase, 'aci2body1CF2.mat'), 'RCF2_Body1_N')
save(fullfile(fileBase, 'aci2body2CF1.mat'), 'RCF1_Body2_N')
save(fullfile(fileBase, 'aci2body2CF2.mat'), 'RCF2_Body2_N')
% save(fullfile(fileBase, 'body12N.mat'), 'RBody1_N')
% save(fullfile(fileBase, 'body22N.mat'), 'RBody2_N')
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