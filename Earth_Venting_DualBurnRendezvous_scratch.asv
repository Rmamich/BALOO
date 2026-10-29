%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Venting Rendezvous Scenario
% Scratch Script to determine burn times and initial offset
%
%
% Created: 03/05/24
% Last Modified: 03/05/24
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('flagSetFlags', 'var')
    if plotCheck
        plotFigsScratch = true;
    else
        plotFigsScratch = false;
    end
    runInitialProp = false;
    rng(97123 + iiTruth+1)
else
    clear all
    close all
    clc
    plotFigsScratch = true;
    runInitialProp = true;
    rng(97123+9)
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

fileBaseBurns = fullfile('.', 'MAT', 'VentingScenario');

%% Design analytically -- maneuver goals 
% Convert to inertial elements
mu = 3.986e5*(1e3)^3; % m^3/s^2


% rp/(1-ecIntial) = acInitial;
rp = (6378+500)*1e3;
eTargetInitial = 0.001;
% Initial Chaser Keplarian Elements
aTargetInitial = rp/(1-eTargetInitial); % m
incTarget = 30*pi/180;  % rad
OmegaTarget = 0; % rad RAAN
omegaTarget = 0; % rad % aop
McTarget = 2*pi-0.2; % rad mean anomaly
eFinal = 0.0;
incFinal = incTarget;% + 0.1*pi/180;
nTargetInitial = sqrt(mu/aTargetInitial^3); % rad/s

dt = 1;
P = 2*pi*sqrt(aTargetInitial^3/mu);
time = 0:dt:(3000);

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
burnTimeInitt = round(time2perapsis - burnTimeTarget);
burnTimeEndt = round(time2perapsis);
dv1 = xcAfterBurn(4:6) - xcBeforeBurn(4:6); % km/s
dv1Finite = dv1/burnTimeTarget;
burnDirt = dv1Finite/norm(dv1Finite);
burnMagt = norm(dv1Finite);






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
burnTimeInitc = round(timeSinceInit2Burn - burnTimeChaser);
burnTimeEndc = round(timeSinceInit2Burn);
dv2 = xtAfterBurn(4:6) - xtBeforeBurn(4:6); % km/s
dv2Finite = dv2/burnTimeChaser;
burnDirc = dv2Finite/norm(dv2Finite);
burnMagc = norm(dv2Finite);



%% Back Prop from final conditions
time2End = time(end)-burnTimeEndc;
[ xcEnd ] = orbitalElmConverter( inputFlag, [aFinal; eFinal; (incFinal*180/pi); OmegaTarget*180/pi; omegaTarget*180/pi; pi], time2End, outputFlag, mu );
[ xdEnd ] = orbitalElmConverter( inputFlag, [aFinal; eFinal; (incFinal*180/pi); OmegaTarget*180/pi; omegaTarget*180/pi; pi+0.00001],time2End, outputFlag, mu );
xcf = xcEnd; % [-6824796.32950823; 1120046.79669146; 646664.452253465; -1410.5469571719; -6445.50890025666; -3721.33264744889];
xdf = xdEnd; % [-6824796.32950823; 1120046.79669146; 646664.452253465; -1410.5469571719; -6445.50890025666; -3721.33264744889];
dt = 1;


 % norm(xcf(1:3) - xdf(1:3))  


%% Venting Info 
% Convert to inertial elements
% The following info is from Power and Prop Element Spacecraft Reference 
% Traj Doc

mass = 39000; % kg -- Power and Prop Element Spacecraft Refereence Traj Doc
CO2mag = 20; % kg*m/s -- Power and Prop Element Spacecraft Refereence Traj Doc
wasteWaterMag = 132; % kg*m/s -- Power and Prop Element Spacecraft Refereence Traj Doc  0.311; % kg m/s^2%

ventMagCO2 = CO2mag/mass/3; % m/s -- 1 sigma
ventMagwasteWater = wasteWaterMag/mass/3; % m/s -- 1 sigma PPE
% ventMagwasteWater = wasteWaterMag/mass; % m/s^2 -- 1 sigma -- ammonia info from dad 

CO2Freq = 10*60; % every 10  min -- we'll say this takes 10 s to do? 
wasteWaterFreq = 6*3600; % every 6 hours -- we'll say this takes 1 min?


CO2Duration = 10; % s
wasteWaterDuration = 60; % s

ventMagCO2 = CO2mag/mass/3/CO2Duration; % m/s^2 -- 1 sigma
ventMagwasteWater = wasteWaterMag/mass/3/wasteWaterDuration; % m/s^2 -- 1 sigma

ventDirMag = norm([-0.5; -0.866; 0]);
CO2ventDir_Body1 = [-0.5; -0.866; 0]/ventDirMag;
wasteWaterDir_Body1 = [-0.5; 0; -0.866]/ventDirMag;


% ventStart = time(1)+100+unidrnd(200);
ventStart = time(1)+unidrnd(300);
CO2VentStart = (ventStart:CO2Freq:time(end))';
wasteWaterVentStart = ((ventStart+1400):wasteWaterFreq:time(end))';

nCO2 = length(CO2VentStart);
nwasteWater = length(wasteWaterVentStart);

CO2VentInfo = [CO2VentStart, CO2VentStart+CO2Duration, abs(randn(nCO2, 1)*ventMagCO2), repmat(CO2ventDir_Body1', [nCO2, 1])];
wasteWaterVentInfo = [wasteWaterVentStart, wasteWaterVentStart+wasteWaterDuration, abs(randn(nwasteWater, 1)*ventMagwasteWater), repmat(wasteWaterDir_Body1', [nwasteWater, 1])];

ventAll = [CO2VentInfo; wasteWaterVentInfo]; % round(P), round(P)+75, 700e-4, -1/sqrt(2), -1/sqrt(2), 0];
ventAll = sortrows(ventAll);


n_vent = length(ventAll(:, 1));
flagDoNotAppend = false;
iiVent = 1;
while iiVent <= n_vent
    if ~flagDoNotAppend
        % Will go here if we have not had overlap between burn and vent
        if (ventAll(iiVent, 1) >= burnTimeInitt &&  ventAll(iiVent, 1) <= burnTimeEndt) || (ventAll(iiVent, 2) >= burnTimeInitt &&  ventAll(iiVent, 2) <= burnTimeEndt)
                % Will go here when there is an overlap on ithe iVent-th vent
                ventAll(iiVent, :) = [burnTimeInitt, burnTimeEndt, burnMagt, burnDirt'];
                flagDoNotAppend = true;
                % TODO: Maybe come back here and combine vent with maneuver?
                % Seems like a lot of work for not much benefit though
        end
    else
        % Will go here after we have already dealt with the overlap
        if (ventAll(iiVent, 1) >= burnTimeInitt &&  ventAll(iiVent, 1) <= burnTimeEndt) || (ventAll(iiVent, 2) >= burnTimeInitt &&  ventAll(iiVent, 2) <= burnTimeEndt)
            % Should never get here... but if two vents happen in the
            % window of the burn it will erase the  vent
            ventAll(iiVent, :) = [];
            n_vent = length(ventAll(:, 1));
        end
    end
    iiVent = iiVent+1;
end

if ~flagDoNotAppend
    % No overlapping Vents
    ventAll = [ventAll; burnTimeInitt, burnTimeEndt, burnMagt, burnDirt'];
end
ventAll = sortrows(ventAll);







%% Back Prop from final conditions
burnStart = {ventAll(:, 1)/dt, burnTimeInitc}; % time step 803 {150/dt, 150/dt}; %
burnEnd = {ventAll(:, 2)/dt; burnTimeEndc}; % time step 806
burn = {ventAll(:, 3); burnMagc};
burnDirection_body = {ventAll(:, 4:6)', burnDirc};


n_burn = [length(burnStart{1}), length(burnStart{2})];



N = length(time);

accelerometer = zeros(3, N);

iplus1 = N;
ii = N-1;
dtProp = 1;
truthState_N = NaN(19, N);
truthState_N(:, N) = [xdEnd; zeros(3, 1);xcEnd; zeros(4, 1)];


ventCount = size(ventAll, 1);

RBody1_N = NaN(3, 3, N);
RBody2_N = NaN(3, 3, N);
pointNoise = 0;
RBody1_N(:, :, iplus1) = cameraPoint(truthState_N(1:6, iplus1), truthState_N(10:15, iplus1), pointNoise);
RBody2_N(:, :, iplus1) = cameraPoint(truthState_N(10:15, iplus1), truthState_N(1:6, iplus1), pointNoise);

for i = 1:(N-1)
    if iplus1 == 181
        hwuhdwe = 0;
    end
    timeInt = time(ii:iplus1);
    burnVecActual = zeros(6, 1);
    burnVecCommand = zeros(6, 1);

    if timeInt(2) < ventAll(ventCount, 2) && timeInt(1) >= ventAll(ventCount, 1)
        % timeInt(1) < ventAll(ventCount, 2) && timeInt(2) >= ventAll(ventCount, 1)
        burnVecActual(1:3) = RBody1_N(:, :, iplus1)'*ventAll(ventCount, 4:6)'*ventAll(ventCount, 3);
        burnVecCommand(1:3) = RBody1_N(:, :, iplus1)'*ventAll(ventCount, 4:6)'*ventAll(ventCount, 3);
        % if iplus1 <= 200 && iplus1 >= 110
            % fprintf('Weird thing happening Target Burining %d, %d\n', timeInt)
        % end
    elseif ventCount > 1 && timeInt(2) < ventAll(ventCount, 1)
        ventCount = ventCount - 1;
        % fprintf('\n')
    end

    if timeInt(2) < burnTimeEndc && timeInt(1) >= burnTimeInitc
        burnVecActual(4:6) = RBody2_N(:, :, iplus1)*burnDirc*burnMagc;
        burnVecCommand(4:6) = RBody2_N(:, :, iplus1)*burnDirc*burnMagc;
        % fprintf('Chaser Burining %d, %d\n', timeInt)
    end

    burnHistActual(:, i) = burnVecActual;
    burnHistCommand(:, i) = burnVecCommand;
    timeProp = flip(timeInt(1):dtProp:timeInt(2));
    Nprop = length(timeProp)-1;
    initialState = truthState_N(:, iplus1);
    initialState(7:9) = burnVecActual(1:3);
    accelerometer(:, iplus1) = burnVecActual(4:6) + initialState(16:18); % + Qchol_c*randn(3, 1);
    for j = 1:Nprop
        state = rk4(initialState, timeProp(j:j+1), eomTruth, gravParamTruth, [0, 0, 0], accelerometer(:, iplus1), 0, RBody2_N(:, :, iplus1), eye(3), [9, 9], {});
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
    % 
    % figure(204)
    % hold on
    % % for k = 1:3
    % %     subplot(3, 1, k)
    %     hold on
    %     grid minor
    %     % plot(time, truthState_N(k, :)-truthState_N(k+9, :))
    %     title('Err Plot')
    % % end

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
ventCount = 1;
allVent = size(ventAll, 1);
for i = 2:N
    if i == 178
        fnejkwfwhe = 0;
    end

    timeInt = time(iminus1:i);
    burnVecActual = zeros(6, 1);
    burnVecCommand = zeros(6, 1);

    if timeInt(1) < burnTimeEndc && timeInt(2)>= burnTimeInitc
        burnVecActual(4:6) = RBody2_N(:, :, i)*burnDirc*burnMagc;
        burnVecCommand(4:6) = RBody2_N(:, :, i)*burnDirc*burnMagc;
        % fprintf('Chaser Burining %d, %d\n', timeInt)
    else
        fjwip = 0;
    end

    if timeInt(1) >= ventAll(ventCount, 1) && timeInt(2) < ventAll(ventCount, 2)
        burnVecActual(1:3) = RBody1_N(:, :, i)'*ventAll(ventCount, 4:6)'*ventAll(ventCount, 3);
        burnVecCommand(1:3) = RBody1_N(:, :, i)'*ventAll(ventCount, 4:6)'*ventAll(ventCount, 3);
        % fprintf('Target Burining %d, %d\n', timeInt)
        % if i <= 200 && i >= 119
        %     fprintf('Target Burining %d, %d\n', timeInt)
        %     burnVecActual(1:3)
        % end
    elseif ventCount < allVent && timeInt(1) > ventAll(ventCount, 2) 
        ventCount = ventCount + 1;
        % fprintf('\n')
    end

    burnHistActual(:, i) = burnVecActual;
    burnHistCommand(:, i) = burnVecCommand;
    timeProp = timeInt(1):dtProp:timeInt(2);
    Nprop = length(timeProp)-1;
    initialState = truthState_N_new(:, iminus1);
    % accelerometer(:, iminus1) = burnVecActual(4:6) + initialState(16:18); % + Qchol_c*randn(3, 1);
    initialState(7:9) = burnVecActual(1:3);
    for j = 1:Nprop
        state = rk4(initialState, timeProp(j:j+1), eomTruth, gravParamTruth, [0, 0, 0], accelerometer(:, i), 0, RBody2_N(:, :, i), eye(3), [9, 9], {});
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

    targetMinusChaser = truthState_N(1:6, :) - truthState_N(10:15, :); 
    targetMinusChaser_new = truthState_N_new(1:6, :) - truthState_N_new(10:15, :); 

    rangeMeas = NaN(1, N);
    for k = 1:N
        rangeMeas(k) = norm(truthState_N(1:3, k)-truthState_N(10:12, k));
    end

    figure(208)
    hold on
    grid minor 
    plot(time, rangeMeas, 'b')

    figure(204)
    hold on
    for k = 1:3
        subplot(3, 1, k)
        hold on
        grid minor
        plot(time, targetMinusChaser(k, :))
        plot(time, targetMinusChaser_new(k, :), '--')
        title('Err Plot')
    end

    figure(205)
    hold on 
    for k = 1:3
        subplot(3, 1, k)
        hold on
        % grid minor
        plot((time), (truthState_N_new(k, :)-truthState_N(k, :)), '--r')
    end

    figure(206)
    hold on 
    for k = 1:3
        subplot(3, 1, k)
        hold on
        % grid minor
        plot((time), (truthState_N_new(k+3, :)-truthState_N(k+3, :)), '--r')
    end
end

burnHistActual_scratch = burnHistActual;

save(fullfile(fileBaseBurns, 'truth_forwardBackward.mat'), 'truthState_N_new');
save(fullfile(fileBaseBurns, 'burnMagRendezvous.mat'), 'burn')
save(fullfile(fileBaseBurns, 'burnDirRendezvous.mat'), "burnDirection_body")
save(fullfile(fileBaseBurns, 'burnStartRendezvous.mat'), "burnStart")
save(fullfile(fileBaseBurns, 'burnEndRendezvous.mat'), "burnEnd")
save(fullfile(fileBaseBurns, 'timeRendezvous.mat'), "time")
save(fullfile(fileBaseBurns, 'rbody1_n.mat'), 'RBody1_N')
save(fullfile(fileBaseBurns, 'rbody2_n.mat'), 'RBody2_N')






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