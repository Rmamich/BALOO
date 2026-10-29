%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dual Burn Rendezvous Scenario -- J2-J8
% Scratch Script to determine burn times and initial offset
%
%
% Created: 11/11/23
% Last Modified: 11/11/23
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('flagSetFlags', 'var')
    plotFigsScratch = false;
    runInitialProp = false;
    rng(97123 + iiTruth+1)
else
    clear all
    close all
    clc
    plotFigsScratch = true;
    runInitialProp = true;
    rng(97123)
    beta = [0; 0; 0];
    gate = 1e-7;
    eomTruth = @tbpJ2_J6_accelerometerGating_rpod_lidarBias;
    mu = 3.986e5*(1e3)^3; % m^3/s^2
    J2 = .0010826269;
    J3 = -0.0000025323;
    J4 = -0.0000016204;
    J5 = -2.27753590730836e-7;
    J6 = 5.40666576283813e-7;
    J7 = -3.50551795713742e-7;
    J8 = -2.03993125929884e-7;
    RBody = 6378*1e3; % m

    gravParamTruth = [mu, RBody, J2, J3, J4, J5, J6, 0, 0, 0];

end

fileBaseBurns = fullfile('.', 'MAT', 'RendezvousScenario');


%% Design analytically
% Convert to inertial elements
mu = 3.986e5*(1e3)^3; % m^3/s^2


% rp/(1-ecIntial) = acInitial;
rp = (6378+500)*1e3;
eTargetInitial = 0.01;
% Initial Chaser Keplarian Elements
aTargetInitial = rp/(1-eTargetInitial); % m
incTarget = 30*pi/180;  % rad
OmegaTarget = 0; % rad RAAN
omegaTarget = 0; % rad % aop
McTarget = 2*pi-0.1; % rad mean anomaly
eFinal = 0.0;
incFinal = incTarget;% + 0.1*pi/180;
nTargetInitial = sqrt(mu/aTargetInitial^3); % rad/s

dt = 1;
P = 2*pi*sqrt(aTargetInitial^3/mu);
time = 0:dt:(3500);

dMTarget = 2*pi - McTarget;
time2perapsis = dMTarget/nTargetInitial; % rad/(rad)*s = s -- start of sim to first burn


aFinal = rp/(1-eFinal);
PFinal = 2*pi*sqrt(aFinal^3/mu); % s

inputFlag = "kep";
outputFlag = "cartesian";

[ xcBeforeBurn ] = orbitalElmConverter( inputFlag, [aTargetInitial; eTargetInitial; incTarget*180/pi; OmegaTarget*180/pi; omegaTarget*180/pi; 0], 0, outputFlag, mu );
[ xcAfterBurn ] = orbitalElmConverter( inputFlag, [aFinal; eFinal; (incFinal*180/pi); OmegaTarget*180/pi; omegaTarget*180/pi; 0], 0, outputFlag, mu );

% Chaser Burn Information
burnTimeTarget = 60; % s
burnTimeInitTarget = round(time2perapsis - burnTimeTarget);
burnTimeEndTarget = burnTimeTarget +burnTimeInitTarget;%+ round(time2perapsis);
dv1 = xcAfterBurn(4:6) - xcBeforeBurn(4:6); % km/s
dv1Finite = dv1/burnTimeTarget;
burnDirTarget = dv1Finite/norm(dv1Finite);
burnMagTarget = norm(dv1Finite);






% Initial Target Keplarian Elements
eChaserInitial = 0.000001; % m
ra = rp; % ra = atInitial*(1+ecInitial); % aFinal 
aChaserInitial = ra/(1 + eChaserInitial);%  = atInitial*(1+ecInitial)
% incc = 30*pi/180;  % rad
% Omegac = 0; % rad RAAN
% omegac = 0; % rad % aop

nChaserInitial = sqrt(mu/aChaserInitial^3); % rad/s
PChaserInitial = 2*pi/nChaserInitial;


% Target will need to travel nInitial*PFinal radians to burn to rendezvous
% at periapsis
timeSinceInit2Burn = time2perapsis + PFinal*1/2; % total time since beginning of sim to target burn
dRad = nChaserInitial*timeSinceInit2Burn; % rad travled from beginning
Mt0 = 2*pi- dRad; % mean anomaly of target initially


inputFlag = "kep";
outputFlag = "cartesian";

[ xtBeforeBurn ] = orbitalElmConverter( inputFlag, [aChaserInitial; eChaserInitial; incTarget*180/pi; OmegaTarget*180/pi; omegaTarget*180/pi; pi], 0, outputFlag, mu );
[ xtAfterBurn ] = orbitalElmConverter( inputFlag, [aFinal; eFinal; (incFinal*180/pi); OmegaTarget*180/pi; omegaTarget*180/pi; pi], 0, outputFlag, mu );

% Chaser Burn Information
burnTimeChaser = 30; % s
burnTimeInitChaser = round(timeSinceInit2Burn - burnTimeChaser);
burnTimeEndChaser = round(timeSinceInit2Burn);
dv2 = xtAfterBurn(4:6) - xtBeforeBurn(4:6); % km/s
dv2Finite = dv2/burnTimeChaser;
burnDirChaser = dv2Finite/norm(dv2Finite);
burnMagChaser = norm(dv2Finite);



%% Back Prop from final conditions
time2End = time(end)-burnTimeEndChaser;
[ xcEnd ] = orbitalElmConverter( inputFlag, [aFinal; eFinal; (incFinal*180/pi); OmegaTarget*180/pi; omegaTarget*180/pi; pi], time2End, outputFlag, mu );
[ xdEnd ] = orbitalElmConverter( inputFlag, [aFinal; eFinal; (incFinal*180/pi); OmegaTarget*180/pi; omegaTarget*180/pi; pi+0.0001],time2End, outputFlag, mu );
xcf = xcEnd; % [-6824796.32950823; 1120046.79669146; 646664.452253465; -1410.5469571719; -6445.50890025666; -3721.33264744889];
xdf = xdEnd; % [-6824796.32950823; 1120046.79669146; 646664.452253465; -1410.5469571719; -6445.50890025666; -3721.33264744889];
dt = 1;

burnStart = {round(burnTimeInitTarget/dt), round(burnTimeInitChaser/dt)}; % time step 803 {150/dt, 150/dt}; %
burnEnd = {burnStart{1} + burnTimeTarget/dt, burnStart{2} + burnTimeChaser/dt}; % time step 806
burn = {burnMagTarget, burnMagChaser};
burnMagSigma = 0.15;
burnDirSigma = 0.5; % deg
burnDirection = {burnDirTarget, burnDirChaser}; % in inertial frame

n_burn = [length(burnStart{1}), length(burnStart{2})];
% eomTruth = @tbpJ2_J8_accelerometerGating_rpod_lidarBias;
% mu = 3.986e5*(1e3)^3; % m^3/s^2
% J2 = .0010826269;
% J3 = -0.0000025323;
% J4 = -0.0000016204;
% J5 = -2.27753590730836e-7;
% J6 = 5.40666576283813e-7;
% J7 = -3.50551795713742e-7;
% J8 = -2.03993125929884e-7;
% RBody = 6378*1e3; % m
% 
% gravParamTruth = [mu, RBody, J2, J3, J4, J5, J6, 0, 0, 0];



N = length(time);
% beta = [0; 0; 0];
% gate = 1e-7;

iplus1 = N;
ii = N-1;
dtProp = 1;
truthState_N = NaN(19, N);
truthState_N(:, N) = [xdf; zeros(3, 1);xcf; zeros(4, 1)];
% 
RBody1_N = NaN(3, 3, N);
RBody2_N = NaN(3, 3, N);
pointNoise = 0;
RBody1_N(:, :, iplus1) = cameraPoint(truthState_N(1:6, iplus1), truthState_N(10:15, iplus1), pointNoise);
RBody2_N(:, :, iplus1) = cameraPoint(truthState_N(10:15, iplus1), truthState_N(1:6, iplus1), pointNoise);

for i = 1:(N-1)
    timeInt = time(ii:iplus1);
    burnVecActual = zeros(6, 1);
    burnVecCommand = zeros(6, 1);

    if timeInt(1) < burnTimeEndChaser && timeInt(2)>= burnTimeInitChaser
        burnVecActual(4:6) = burnDirChaser*burnMagChaser;
        burnVecCommand(4:6) = burnDirChaser*burnMagChaser;
        % fprintf('Chaser Burining %d, %d\n', timeInt)
    else
        fjwip = 0;
    end

    if timeInt(1) < burnTimeEndTarget && timeInt(2) >= burnTimeInitTarget
        burnVecActual(1:3) = burnDirTarget*burnMagTarget;
        burnVecCommand(1:3) = burnDirTarget*burnMagTarget;
        % fprintf('Target Burining %d, %d\n', timeInt)
    end

    burnHistActual(:, i) = burnVecActual;
    burnHistCommand(:, i) = burnVecCommand;
    timeProp = flip(timeInt(1):dtProp:timeInt(2));
    Nprop = length(timeProp)-1;
    initialState = truthState_N(:, iplus1);
    initialState(7:9) = burnVecActual(1:3);
    accelerometer(:, i) = burnVecActual(4:6) + initialState(16:18); % + Qchol_c*randn(3, 1);
    for j = 1:Nprop
        state = rk4(initialState, timeProp(j:j+1), eomTruth, gravParamTruth, beta, accelerometer(:, i), gate, eye(3), eye(3), [9, 9], {});
        initialState = state;
    end
    truthState_N(:, ii) = state;
    RBody1_N(:, :, ii) = cameraPoint(truthState_N(1:6, ii), truthState_N(10:15, ii), pointNoise);
    RBody2_N(:, :, ii) = cameraPoint(truthState_N(10:15, ii), truthState_N(1:6, ii), pointNoise);
    iplus1 = N-i;
    ii = iplus1-1;
    % if mod(i, 1000) == 0
    %     fprintf('%d\n', i)
    % end
end

trim = 0;

truthState_N(:, 1:trim) = [];
time((N-trim+1):end) = [];
N = length(time);
burnTimeInitTarget = burnTimeInitTarget-trim;
burnTimeInitChaser = burnTimeInitChaser-trim;
burnTimeEndTarget = burnTimeEndTarget-trim;
burnTimeEndChaser = burnTimeEndChaser-trim;
burnStart = {round((burnTimeInitTarget)/dt), round((burnTimeInitChaser)/dt)}; % time step 803 {150/dt, 150/dt}; %
burnEnd = {burnStart{1} + burnTimeTarget/dt, burnStart{2} + burnTimeChaser/dt}; % time step 806


if plotFigsScratch
    figure(201)
    hold on
    plot3(truthState_N(1, :), truthState_N(2, :), truthState_N(3, :))
    grid minor
    plot3(truthState_N(10, :), truthState_N(11, :), truthState_N(12, :))

    figure(203)
    hold on
    for k = 1:3
        subplot(3, 1, k)
        hold on
        grid minor
        plot(time, truthState_N(k, :))
        plot(time, truthState_N(k+9, :))
        title('Both Plotted')
    end

    figure(204)
    hold on
    for k = 1:3
        subplot(3, 1, k)
        hold on
        grid minor
        plot(time, truthState_N(k, :)-truthState_N(k+9, :))
        title('Err Plot')
    end
end



outputFlag = "kep";
inputFlag = "cartesian";

[ oec ] = orbitalElmConverter( inputFlag, truthState_N(1:6, 1), 0, outputFlag, mu );
[ oed ] = orbitalElmConverter( inputFlag, truthState_N(10:15, 1), 0, outputFlag, mu );



% fprintf('\n\n\n')

iminus1 = 1;
dtProp = 1;
truthState_N_new = NaN(19, N);
truthState_N_new(:, 1) = truthState_N(:, 1);
for i = 2:(N)
    timeInt = time(iminus1:i);
    burnVecActual = zeros(6, 1);
    burnVecCommand = zeros(6, 1);

    if timeInt(1) < burnTimeEndChaser && timeInt(2)>= burnTimeInitChaser
        burnVecActual(4:6) = burnDirChaser*burnMagChaser;
        burnVecCommand(4:6) = burnDirChaser*burnMagChaser;
        % fprintf('Chaser Burining %d, %d\n', timeInt)
    else
        fjwip = 0;
    end

    if timeInt(1) < burnTimeEndTarget && timeInt(2) >= burnTimeInitTarget
        burnVecActual(1:3) = burnDirTarget*burnMagTarget;
        burnVecCommand(1:3) = burnDirTarget*burnMagTarget;
        % fprintf('Target Burining %d, %d\n', timeInt)
    end

    burnHistActual(:, i) = burnVecActual;
    burnHistCommand(:, i) = burnVecCommand;
    timeProp = timeInt(1):dtProp:timeInt(2);
    Nprop = length(timeProp)-1;
    initialState = truthState_N_new(:, iminus1);
    accelerometer(:, i) = burnVecActual(4:6) + initialState(16:18); % + Qchol_c*randn(3, 1);
    initialState(7:9) = burnVecActual(1:3);
    for j = 1:Nprop
        state = rk4(initialState, timeProp(j:j+1), eomTruth, gravParamTruth, beta, accelerometer(:, i), gate, eye(3), eye(3), [9, 9], {});
        initialState = state;
    end
    truthState_N_new(:, i) = state;
    iminus1 = i;
    % if mod(i, 1000) == 0
    %     fprintf('%d\n', i)
    % end
end


if plotFigsScratch
    figure(201)
    hold on
    plot3(truthState_N_new(1, :), truthState_N_new(2, :), truthState_N_new(3, :), '--')
    plot3(truthState_N_new(10, :), truthState_N_new(11, :), truthState_N_new(12, :), '--')

    figure(203)
    hold on
    for k = 1:3
        subplot(3, 1, k)
        hold on
        grid minor
        plot(time, truthState_N_new(k, :), '--')
        plot(time, truthState_N_new(k+9, :), '--')
    end

    figure(204)
    hold on
    for k = 1:3
        subplot(3, 1, k)
        hold on
        grid minor
        plot(time, truthState_N_new(k, :)-truthState_N_new(k+9, :), '--')
    end

     rangeMeas = NaN(1, N);
    for k = 1:N
        rangeMeas(k) = norm(truthState_N(1:3, k)-truthState_N(10:12, k));
    end

    figure(208)
    hold on
    grid minor 
    plot(time, rangeMeas, 'b')
end
save(fullfile(fileBaseBurns, 'truth_forwardBackward.mat'), 'truthState_N_new');
save(fullfile(fileBaseBurns, 'burnMagRendezvous.mat'), 'burn')
save(fullfile(fileBaseBurns, 'burnDirRendezvous.mat'), "burnDirection")
save(fullfile(fileBaseBurns, 'burnStartRendezvous.mat'), "burnStart")
save(fullfile(fileBaseBurns, 'burnEndRendezvous.mat'), "burnEnd")
save(fullfile(fileBaseBurns, 'timeRendezvous.mat'), "time")
save(fullfile(fileBaseBurns, 'rbody1_n.mat'), "RBody1_N")
save(fullfile(fileBaseBurns, 'rbody2_n.mat'), "RBody2_N")

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