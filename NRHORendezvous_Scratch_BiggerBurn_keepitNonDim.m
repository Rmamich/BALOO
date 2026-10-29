%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Venting Rendezvous Scenario
% Scratch Script to determine burn times and initial offset
% Setting everything in dimensional, inertial space and converting back and
% forth before and after prop
%
%
% Created: 03/05/24
% Last Modified: 08/06/24
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
    beta = zeros(3, 1);
    L = 384400; % km -- characteristic length

    tCharacteristic = 3.751902619517e5; % s -- characteristic time from paper
    beta = [0; 0; 0]; % 1/s * s
    beta = beta * tCharacteristic; % nondim
    gate = 1e-4; % m/s/s
    gate = gate/(L*1000)*tCharacteristic*tCharacteristic; % [-]
    fontSize = 25;                                              % Sets font size for plotting script plots
    line_width = 4;                                             % Sets line width for plotting script plots
    textInterpreter = 'latex';                                  % Sets text interpreter for plotting script plots
    axesLineWidth = 2;                                          % Sets axes width for plotting script plots
    gridLineWidth = 1;                                          % Sets grid line width for plotting script plots
    markerSize = 15;                                             % Sets marker size for plotting script plots
    set(0, 'DefaultTextInterpreter', 'latex')                   % Sets text interpreter for plotting script plots
    set(0, 'DefaultLegendInterpreter', 'latex')                 % Sets text interpreter for plotting script plots
    set(0, 'DefaultTextFontSize', fontSize)                     % Sets font size for plotting script plots
    set(0, 'DefaultLineLineWidth', line_width)                  % Sets line width for plotting script plots
    set(0, 'DefaultAxesTickLabelInterpreter', textInterpreter)
end

%% Commented out all vents and burns -- test with perfect models 8/26/24 - 4:15pm
% Changed final spread to be further apart to keep covarianve normal


fileBaseBurns = fullfile('.', 'MAT', 'NRHOVentingScenario_nonDim');




%% Forward Prop -- get final approximate states for target/chaser
% From Lorin's AAS GNC Paper 2024
mu = 0.012150584270; % mass ratio from paper
L = 384400; % km -- characteristic length
L_m = L*1000; % m 
% gate = 0;

rNondimChaser = [1.022023976774, 0, -0.182098475077]'; % nonimensaional initial conditions from paper
vNondimChaser = [0, -0.103261718478, 0]'; % nondimensional initial conditions from paper

rDimChaser = L*rNondimChaser; % makign position dimensional -- r nondimensional * characteristic length

PNondim = 1.511143593137; % period non dimensional from paper
tCharacteristic = 3.751902619517e5; % s -- characteristic time from paper

tDim = PNondim*tCharacteristic; % dimensional period -- convert as nondimensional period * characteristic time
P = tDim; %
n = 2*pi/P; % 1/s
% 0.001 rad/ n

vDim = vNondimChaser*L*n; % m/s -- convert velocity to dimensional

tDim_prop = 3000; % my desired propagation time for MILTON
tNondim = tDim_prop/tCharacteristic; % convert to nondimensional

% RK4 Test -- testing to see if fixed step rk4 even works here
dt = 1; % dimensional dt
dtNonDim = dt/tCharacteristic;





% timeNonDim = 0:dtNonDim:tNondim;
% time = timeNonDim*tCharacteristic;

time = 0:dt:P;
timeNonDim = time/tCharacteristic;
N = length(timeNonDim);
initialState = [rNondimChaser; vNondimChaser];
eomProp = @CRTBP_targetChaser;
gravParamTruth = [mu, L, n, tCharacteristic];


if runInitialProp
    nonDimStateTimeHist = zeros(19, N);
    nonDimStateTimeHist(1:6, 1) = initialState;
    nonDimStateTimeHist(10:15, 1) = initialState;

    initialState = nonDimStateTimeHist(:, 1);
    for j = 2:N
        stateRK4prop = rk4(initialState, timeNonDim(j-1:j), eomProp, gravParamTruth, beta, zeros(3, 1), 0, eye(3), eye(3), [9, 9], {});
        nonDimStateTimeHist(:, j) = stateRK4prop;
        initialState = stateRK4prop;
    end





    dimStateTimeHist = zeros(19, N);
    nonDimStateTimeHist_goingBack = zeros(19, N);
    % convert to inertial and back
    for j = 1:N
        dimStateTimeHist(:, j) = nonDimRotational_dimInertial_conversion(nonDimStateTimeHist(:, j), time(j), L_m, n, tCharacteristic, 'fullState_nondimR2dimN');
        nonDimStateTimeHist_goingBack(:, j) = nonDimRotational_dimInertial_conversion(dimStateTimeHist(:, j), time(j), L_m, n, tCharacteristic, 'fullState_dimN2nondimR');
    end

    if plotFigsScratch
        
        figure(1)
        plot3(nonDimStateTimeHist(1, :), nonDimStateTimeHist(2, :),nonDimStateTimeHist(3, :))
        hold on
        ax = gca;
        % ax.LineWidth = 4;
        ax.GridLineWidth = 1;
        ax.TickLabelInterpreter = 'latex';
        set(gca, 'FontSize', fontSize)
        % plot3(nonDimStateTimeHist(10, :), nonDimStateTimeHist(11, :), nonDimStateTimeHist(12, :))
        % plot3(nonDimStateTimeHist_goingBack(1, :), nonDimStateTimeHist_goingBack(2, :),nonDimStateTimeHist_goingBack(3, :))
        % plot3(nonDimStateTimeHist_goingBack(10, :), nonDimStateTimeHist_goingBack(11, :), nonDimStateTimeHist_goingBack(12, :))
        % title(strcat('Plot NRHO in the Rotating Frame'))
        % xlabel('Rotating X [-]', 'Interpreter','latex', 'FontSize',fontSize)
        % ylabel('Rotating Y [-]', 'Interpreter','latex', 'FontSize',fontSize)
        % zlabel('Rotating Z [-]', 'Interpreter','latex', 'FontSize',fontSize)
        axis equal
        grid minor
        [X, Y, Z] = sphere(100);
        radL = 1740/L;
        X = radL*X;
        X = X + 1 - mu;
        Y = radL*Y;
        Z = radL*Z;
        % plot3(-mu, 0, 0, '.', 'MarkerSize', 10) % Earth
        texture = imread('Milton_Moon_Texture.jpg');
        % texture = imread('texture.jpg');
        surface(X,Y,-Z,'FaceColor','texturemap',... % nice moon surface
        'EdgeColor','none', 'Cdata',texture)
        % plot3(1 - mu, 0, 0, '.', 'MarkerSize', 10) % moon 
        % legend('S/C', 'Moon')
        % legend('S/C', 'Earth', 'Moon')
        % exportgraphics(figure(1), fullfile('Figures' ,'NRHO.png'), 'resolution', 300)
        % exportgraphics(figure(1), fullfile('Figures' ,'NRHO.eps'), 'ContentType', 'vector')





        figure(28)
        plot3(nonDimStateTimeHist(1, :), nonDimStateTimeHist(3, :),nonDimStateTimeHist(2, :))
        hold on
        ax = gca;
        % ax.LineWidth = 4;
        ax.GridLineWidth = 1;
        ax.TickLabelInterpreter = 'latex';
        set(gca, 'FontSize', fontSize)
        % plot3(nonDimStateTimeHist(10, :), nonDimStateTimeHist(11, :), nonDimStateTimeHist(12, :))
        % plot3(nonDimStateTimeHist_goingBack(1, :), nonDimStateTimeHist_goingBack(2, :),nonDimStateTimeHist_goingBack(3, :))
        % plot3(nonDimStateTimeHist_goingBack(10, :), nonDimStateTimeHist_goingBack(11, :), nonDimStateTimeHist_goingBack(12, :))
        % title(strcat('Plot NRHO in the Rotating Frame'))
        % xlabel('Rotating X [-]', 'Interpreter','latex', 'FontSize',fontSize)
        % ylabel('Rotating Y [-]', 'Interpreter','latex', 'FontSize',fontSize)
        % zlabel('Rotating Z [-]', 'Interpreter','latex', 'FontSize',fontSize)
        axis equal
        grid minor
        [X, Y, Z] = sphere(100);
        radL = 1740/L;
        X = radL*X;
        X = X + 1 - mu;
        Y = radL*Y;
        Z = radL*Z;
        % plot3(-mu, 0, 0, '.', 'MarkerSize', 10) % Earth
        texture = imread('Milton_Moon_Texture.jpg');
        % texture = imread('texture.jpg');
        surface(X,-Z,Y,'FaceColor','texturemap',... % nice moon surface
        'EdgeColor','none', 'Cdata',texture)
        % plot3(1 - mu, 0, 0, '.', 'MarkerSize', 10) % moon 
        % legend('S/C', 'Moon')
        % legend('S/C', 'Earth', 'Moon')


        figure(27)
        hold on 
        axis equal
        
        surface(X,Y,-Z,'FaceColor','texturemap',... % nice moon surface
        'EdgeColor','none', 'Cdata',texture)



        figure(2)
        plot3(dimStateTimeHist(1, :), dimStateTimeHist(2, :), dimStateTimeHist(3, :))
        hold on
        plot3(dimStateTimeHist(10, :), dimStateTimeHist(11, :), dimStateTimeHist(12, :))
        title(strcat('Plot of NRHO in the Inertial Frame'))
        xlabel('Inertial X [km]')
        ylabel('Inertial Y [km]')
        zlabel('Inertial Z [km]')
        plot3(0, 0, 0, '.', 'MarkerSize', 10)
        axis equal
        grid minor
        legend('Target', 'Chaser','Earth')
    end

    save(fullfile(fileBaseBurns, 'nonDimStateTimeHist.mat'), 'nonDimStateTimeHist');
else
    load(fullfile(fileBaseBurns, 'nonDimStateTimeHist.mat'));
end

finalSpreadIndex = 200;
idxTargetEnd = 284792+tDim_prop/2;

finalState_target_nonDim_R = nonDimStateTimeHist(1:6, idxTargetEnd); % in "front"
finalState_chaser_nonDim_R = nonDimStateTimeHist(1:6, idxTargetEnd-finalSpreadIndex); % trailing behind


[finalState_target_dim_N] = nonDimRotational_dimInertial_conversion(finalState_target_nonDim_R, timeNonDim(idxTargetEnd), L_m, n, tCharacteristic,'R2N');
[finalState_chaser_dim_N] = nonDimRotational_dimInertial_conversion(finalState_chaser_nonDim_R, timeNonDim(idxTargetEnd), L_m, n, tCharacteristic,'R2N');

finalState_chaser_dim_N = finalState_target_dim_N + [ones(3, 1)*2; zeros(3, 1)];

%% Venting Information/ Chaser burn
% The following info is from Power and Prop Element Spacecraft Reference
% Traj Doc

% time = 0:dt:tDim_prop;

% time = time(idxTargetEnd):dt:P+time(idxTargetEnd);
time = (time(idxTargetEnd)-tDim_prop):dt:time(idxTargetEnd);
timeNonDim = time/tCharacteristic;



mass = 39000; % kg -- Power and Prop Element Spacecraft Refereence Traj Doc
CO2mag = 20*1e-3; % kg*km/s -- Power and Prop Element Spacecraft Refereence Traj Doc
wasteWaterMag = 132*1e-3; % kg*km/s -- Power and Prop Element Spacecraft Refereence Traj Doc  0.311; % kg m/s^2%

ventMagCO2 = CO2mag/mass/3; % lm/s -- 1 sigma
ventMagwasteWater = wasteWaterMag/mass/3; % km/s -- 1 sigma PPE
% ventMagwasteWater = wasteWaterMag/mass; % m/s^2 -- 1 sigma -- ammonia info from dad

CO2Freq = 10*60; % every 10  min -- we'll say this takes 10 s to do?
wasteWaterFreq = 6*3600; % every 6 hours -- we'll say this takes 1 min?


CO2Duration = 10; % [s]
wasteWaterDuration = 60; % [s]

ventMagCO2 = CO2mag/mass/3/CO2Duration/(L)*tCharacteristic^2; % [-] -- kg*km/s / kg / s 1 sigma
ventMagwasteWater = wasteWaterMag/mass/3/wasteWaterDuration/(L)*tCharacteristic^2; % [-] -- 1 sigma

CO2Duration = 10/tCharacteristic; % [-]
wasteWaterDuration = 60/tCharacteristic; % [-]

ventDirMag = norm([-0.5; -0.866; 0]);
CO2ventDir_Body1 = [-0.5; -0.866; 0]/ventDirMag;
wasteWaterDir_Body1 = [-0.5; 0; -0.866]/ventDirMag;


ventStart = time(1)+1+unidrnd(30);
CO2VentStart = ((ventStart):CO2Freq:time(end))'/tCharacteristic; % [-]
wasteWaterVentStart = ((ventStart+1400):wasteWaterFreq:time(end))'/tCharacteristic; % [-]

nCO2 = length(CO2VentStart);
nwasteWater = length(wasteWaterVentStart);

CO2VentInfo = [CO2VentStart, (CO2VentStart+CO2Duration), abs(randn(nCO2, 1)*ventMagCO2), repmat(CO2ventDir_Body1', [nCO2, 1])];
wasteWaterVentInfo = [wasteWaterVentStart, wasteWaterVentStart+wasteWaterDuration, abs(randn(nwasteWater, 1)*ventMagwasteWater), repmat(wasteWaterDir_Body1', [nwasteWater, 1])];

ventAll = [CO2VentInfo; wasteWaterVentInfo]; % round(P), round(P)+75, 700e-4, -1/sqrt(2), -1/sqrt(2), 0];
ventAll = sortrows(ventAll);

%% Chaser and Target Burn
burnTimet = 30/tCharacteristic; % [-]
burnTimeInitt = (time(1)+350)/tCharacteristic; % [-]
burnTimeEndt = (burnTimeInitt+burnTimet); % [-]
% burnDirc = [-1/sqrt(2); -1/sqrt(2); 0];
burnDirt = [1; 1; 1]/sqrt(3);
burnMagt = 0.1/(L_m)*tCharacteristic^2;


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



% Chaser Burn Information -- may need to adjust this
burnTimec = 30/tCharacteristic; % [-]
burnTimeInitc = (time(end)-150)/tCharacteristic; % [-]
burnTimeEndc = (burnTimeInitc+burnTimec); % [-]
% burnDirc = [-1/sqrt(2); -1/sqrt(2); 0];
burnDirc = [1; 1; 1]/sqrt(3);
burnMagc = 1e-3/(L_m)*tCharacteristic^2;



%% Back Prop from final conditions
% ventAll(:, 3) = 0;
burnStart = {ventAll(:, 1), burnTimeInitc}; % nondim time to start
burnEnd = {ventAll(:, 2); burnTimeEndc}; % nondim time to end
burn = {ventAll(:, 3); burnMagc}; % nondim burn magnitude
burnDirection_body = {ventAll(:, 4:6)', burnDirc}; % direction of uburn in body frame



n_burn = [length(burnStart{1}), length(burnStart{2})];

time = timeNonDim;

N = length(time);
% beta = [0; 0; 0]; % 1/s * s
% beta = beta * tCharacteristic; % nondim
% gate = 1e-4*0; % m/s/s
% gate = gate/(L_m)*tCharacteristic^2; % [-]

iplus1 = N;
ii = N-1;
dtProp = 1;
dtProp_nonDim = dtProp/tCharacteristic;
truthState_N = NaN(19, N);
truthState_N(:, N) = [finalState_target_dim_N; zeros(3, 1);finalState_chaser_dim_N; zeros(4, 1)];
truthState_R = NaN(19, N);
truthState_R(:, N) = nonDimRotational_dimInertial_conversion(truthState_N(:, N), time(end)*tCharacteristic, L_m, n, tCharacteristic, 'fullState_dimN2nondimR');
accelerometer = zeros(3, N);

ventCount = size(ventAll, 1);

RBody1_R = NaN(3, 3, N);
RBody2_R = NaN(3, 3, N);
pointNoise = 0;
RBody1_R(:, :, end) = cameraPoint(truthState_R(1:6, iplus1), truthState_R(10:15, iplus1), pointNoise);
RBody2_R(:, :, end) = cameraPoint(truthState_R(10:15, iplus1), truthState_R(1:6, iplus1), pointNoise);


initialState_nonDimR = nonDimRotational_dimInertial_conversion(truthState_N(:, N), time(end)*tCharacteristic, L_m, n, tCharacteristic, 'fullState_dimN2nondimR');

% fprintf('Initial diff: %0.6e, %0.6e, %0.6e, %0.6e, %0.6e, %0.6e\n', initialState_nonDimR(1:6) - finalState_target_nonDim_R)
% fprintf('Initial diff: %0.6e, %0.6e, %0.6e, %0.6e, %0.6e, %0.6e\n', initialState_nonDimR(10:15) - finalState_target_nonDim_R)

% initialState_nonDimR(7:9) = 0;
initialState = truthState_N(:, N);
for i = 1:(N-1)
    timeInt = timeNonDim(ii:iplus1); % dimensional
    burnVecActual = zeros(6, 1);
    burnVecCommand = zeros(6, 1);

    if timeInt(1) < ventAll(ventCount, 2) && timeInt(2)>= ventAll(ventCount, 1)
        % burnVecActual(1:3) = RBody1_R(:, :, iplus1)'*ventAll(ventCount, 4:6)'*ventAll(ventCount, 3); %convert fron body to inertial dimensional
        % burnVecCommand(1:3) = RBody1_R(:, :, iplus1)'*ventAll(ventCount, 4:6)'*ventAll(ventCount, 3); %convert fron body to inertial dimensional
        burnVecActual(1:3) = ventAll(ventCount, 4:6)'*ventAll(ventCount, 3); %convert fron body to inertial dimensional
        burnVecCommand(1:3) = ventAll(ventCount, 4:6)'*ventAll(ventCount, 3); %convert fron body to inertial dimensional
        % fprintf('Target Burining %d, %d\n', timeInt)
    elseif ventCount > 1 && timeInt(2) < ventAll(ventCount, 1)
        ventCount = ventCount - 1;
        % fprintf('\n')
    end

    if timeInt(1) < burnTimeEndc && timeInt(2) >= burnTimeInitc
        burnVecActual(4:6) = RBody2_R(:, :, iplus1)*burnDirc*burnMagc; % from inertial to body frame dimensional
        burnVecCommand(4:6) = RBody2_R(:, :, iplus1)*burnDirc*burnMagc; % from inertial to body frame dimensional
        % burnVecActual(4:6) = burnDirc*burnMagc; % from inertial to body frame dimensional
        % burnVecCommand(4:6) = burnDirc*burnMagc; % from inertial to body frame dimensional
        % fprintf('Chaser Burining %d, %d\n', (timeInt-time(1))*tCharacteristic)
    end

    burnHistActual(:, i) = burnVecActual;
    burnHistCommand(:, i) = burnVecCommand;
    timeProp = flip(timeInt(1):dtProp_nonDim:timeInt(2));
    Nprop = length(timeProp)-1;
    initialState = truthState_R(:, iplus1);
    initialState(7:9) = burnVecActual(1:3);
    accelerometer(:, ii) = burnVecActual(4:6) + initialState(16:18);
    % Convert prop values to nondimensional
    % [finalState_chaser_dim_N] = nonDimRotational_dimInertial_conversion(finalState_chaser_nonDim_R, timeNonDim(1), L_m, n, tCharacteristic,'R2N');
    % timeProp_nonDim = timeProp/tCharacteristic;
    % initialState_nonDimR = truthState_R(:, iplus1);
    % initialState_nonDimR = nonDimRotational_dimInertial_conversion(initialState, timeProp(1), L_m, n, tCharacteristic, 'fullState_dimN2nondimR');
    % initialState_goingBack = nonDimRotational_dimInertial_conversion(initialState_nonDimR, timeProp(1), L_m, n, tCharacteristic, 'fullState_nondimR2dimN');
    % initialState_c_nonDimR = nonDimRotational_dimInertial_conversion(initialState(10:15), timeProp_nonDim(1), L_m, n, tCharacteristic,'N2R');
    % accelerometer_nonDim = nonDimRotational_dimInertial_conversion(accelerometer(:, i), timeProp(1), L_m, n, tCharacteristic, 'accelerometer');
    for j = 1:Nprop
        state = rk4(initialState, timeProp(j:j+1), eomProp, gravParamTruth, beta, accelerometer(:, ii), 0, RBody2_R(:, :, iplus1), eye(3), [9, 9], {});
        initialState = state;
    end
    truthState_R(:, ii) = state;
    truthState_N(:, ii) = nonDimRotational_dimInertial_conversion(truthState_R(:, ii), timeProp(end)*tCharacteristic, L_m, n, tCharacteristic, 'fullState_nondimR2dimN');
    RBody1_R(:, :, ii) = cameraPoint(truthState_R(1:6, ii), truthState_R(10:15, ii), pointNoise);
    RBody2_R(:, :, ii) = cameraPoint(truthState_R(10:15, ii), truthState_R(1:6, ii), pointNoise);
    iplus1 = N-i;
    ii = iplus1-1;
    if mod(i, 100000) == 0
        % fprintf('step %d of %d\n', i, N)
    end
    % if i == 1
    %     accelerometer(:, N) = burnVecActual(4:6) + state(16:18);
    % end

end


% % convert to inertial and back
% for j = 1:N
%     truthState_N(:, j) = nonDimRotational_dimInertial_conversion(nonDimStateTimeHist(:, j), time(j), L_m, n, tCharacteristic, 'fullState_nondimR2dimN');
%     % nonDimStateTimeHist_goingBack(:, j) = nonDimRotational_dimInertial_conversion(dimStateTimeHist(:, j), time(j), L_m, n, tCharacteristic, 'fullState_dimN2nondimR');
% end

if plotFigsScratch
    figure(3)
    plot3(truthState_R(1, :), truthState_R(2, :), truthState_R(3, :))
    hold on
    plot3(truthState_R(10, :), truthState_R(11, :), truthState_R(12, :))
    title(strcat('Plot NRHO in the Rotating Frame'))
    xlabel('Rotating X []')
    ylabel('Rotating Y []')
    zlabel('Rotating Z []')
    axis equal
    grid minor
    plot3(-mu, 0, 0, '.', 'MarkerSize', 10)
    plot3(1 - mu, 0, 0, '.', 'MarkerSize', 10)
    legend('Target', 'Chaser', 'Earth', 'Moon')


    figure(4)
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
    % legend('complexProp Target', 'complexProp Chaser', 'simpleProp Target', 'simpleProp Chaser', 'Earth')÷
    legend('complexProp Target', 'complexProp Chaser', 'Earth')


    figure(203)
    hold on
    for k = 1:3
        subplot(3, 1, k)
        hold on
        grid minor
        plot((time-time(1))*tCharacteristic, truthState_N(k, :), 'b')
        plot((time-time(1))*tCharacteristic, truthState_N(k+9, :), 'r')
        title('Both Plotted')
    end

    figure(204)
    hold on
    for k = 1:3
        subplot(3, 1, k)
        hold on
        grid minor
        plot((time-time(1))*tCharacteristic, (truthState_N(k, :)-truthState_N(k+9, :)), 'b')
        title('Err Plot')
    end
    outputFlag = "kep";
    inputFlag = "cartesian";
    
    rangeMeas = NaN(1, N);
    for k = 1:N
        rangeMeas(k) = norm(truthState_N(1:3, k)-truthState_N(10:12, k));
    end

    figure(208)
    hold on
    grid minor 
    plot((time-time(1))*tCharacteristic, rangeMeas, 'b')
end



% fprintf('\n\n\n')

iminus1 = 1;
dtProp = 1;
truthState_N_new = NaN(19, N);
truthState_N_new(:, 1) =   truthState_N(:, 1); % dimStateTimeHist(:, 1); %

truthState_R_new = NaN(19, N);
truthState_R_new(:, 1) =   truthState_R(:, 1); % nonDimStateTimeHist(:, 1);
ventCount = 1;
initialState_nonDimR = truthState_R_new(:, 1);
allVent = size(ventAll, 1);

for i = 2:N
    timeInt = timeNonDim(iminus1:i);
    burnVecActual = zeros(6, 1);
    burnVecCommand = zeros(6, 1);


    if timeInt(1) < burnTimeEndc && timeInt(2)>= burnTimeInitc
        % Need to fix that the EOM think this will be in inertial -- not
        % rotating frame
        burnVecActual(4:6) = RBody2_R(:, :, i)*burnDirc*burnMagc;
        burnVecCommand(4:6) = RBody2_R(:, :, i)*burnDirc*burnMagc;
        % burnVecActual(4:6) = burnDirc*burnMagc;
        % burnVecCommand(4:6) = burnDirc*burnMagc;
        % fprintf('Chaser Burining %d, %d\n', iminus1:i)
    else
        fjwip = 0;
    end

    if timeInt(1) < ventAll(ventCount, 2) && timeInt(2) >= ventAll(ventCount, 1)
        % Need to fix that the EOM think this will be in inertial -- not
        % rotating frame
        % burnVecActual(1:3) = RBody1_R(:, :, i)'*ventAll(ventCount, 4:6)'*ventAll(ventCount, 3);
        % burnVecCommand(1:3) = RBody1_R(:, :, i)'*ventAll(ventCount, 4:6)'*ventAll(ventCount, 3);
        burnVecActual(1:3) = ventAll(ventCount, 4:6)'*ventAll(ventCount, 3);
        burnVecCommand(1:3) = ventAll(ventCount, 4:6)'*ventAll(ventCount, 3);
        % fprintf('Target Burining %d, %d\n', timeInt)
    elseif ventCount < allVent && timeInt(1) > ventAll(ventCount, 2)
        ventCount = ventCount + 1;
        % fprintf('\n')
    end

    burnHistActual(:, i) = burnVecActual;
    burnHistCommand(:, i) = burnVecCommand;
    timeProp = timeInt(1):dtProp_nonDim:timeInt(2);
    Nprop = length(timeProp)-1;
    initialState = truthState_R_new(:, iminus1);
    % accelerometer(:, iminus1) = burnVecActual(4:6) + initialState(16:18); % + Qchol_c*randn(3, 1);
    % accelerometer(:, iminus1) = burnVecActual(4:6) + initialState(16:18);
    initialState(7:9) = burnVecActual(1:3);

    % timeProp_nonDim = timeProp/tCharacteristic;
    % initialState_nonDimR_new = nonDimRotational_dimInertial_conversion(initialState, timeProp(1), L_m, n, tCharacteristic, 'fullState_dimN2nondimR');
    % initialState_nonDimR = initialState_nonDimR_new;
    % accelerometer_nonDim = nonDimRotational_dimInertial_conversion(accelerometer(:, i), timeProp(1), L_m, n, tCharacteristic, 'accelerometer');
    % if norm(initialState_nonDimR_new - initialState_nonDimR) > 1e-15
    %     cdwo = 0;
    % end
    for j = 1:Nprop
        state = rk4(initialState, timeProp(j:j+1), eomProp, gravParamTruth, beta, accelerometer(:, i), 0, RBody2_R(:, :, i), eye(3), [9, 9], {});
        initialState = state;
    end
    truthState_R_new(:, i) = initialState;
    truthState_N_new(:, i) = nonDimRotational_dimInertial_conversion(truthState_R_new(:, i), timeProp(end)*tCharacteristic, L_m, n, tCharacteristic, 'fullState_nondimR2dimN');
    iminus1 = i;
    if mod(i, 100000) == 0
        % fprintf('%d\n', i)
    end
end




if plotFigsScratch
    figure(5)
    plot3(truthState_R_new(1, :), truthState_R_new(2, :), truthState_R_new(3, :))
    hold on
    plot3(truthState_R_new(10, :), truthState_R_new(11, :), truthState_R_new(12, :))
    title(strcat('Plot NRHO in the Rotating Frame'))
    xlabel('Rotating X []')
    ylabel('Rotating Y []')
    zlabel('Rotating Z []')
    axis equal
    grid minor
    plot3(-mu, 0, 0, '.', 'MarkerSize', 10)
    plot3(1 - mu, 0, 0, '.', 'MarkerSize', 10)
    legend('Target', 'Chaser', 'Earth', 'Moon')


    figure(6)
    plot3(truthState_N_new(1, :), truthState_N_new(2, :), truthState_N_new(3, :))
    hold on
    plot3(truthState_N_new(10, :), truthState_N_new(11, :), truthState_N_new(12, :))
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
    % legend('complexProp Target', 'complexProp Chaser', 'simpleProp Target', 'simpleProp Chaser', 'Earth')÷
    legend('complexProp Target', 'complexProp Chaser', 'Earth')




    figure(201)
    hold on
    plot3(truthState_N_new(1, :), truthState_N_new(2, :), truthState_N_new(3, :), '--')
    plot3(truthState_N_new(10, :), truthState_N_new(11, :), truthState_N_new(12, :), '--')

    figure(203)
    hold on
    for k = 1:3
        subplot(3, 1, k)
        hold on
        % grid minor
        plot((time-time(1))*tCharacteristic, truthState_N_new(k, :), '--c')
        plot((time-time(1))*tCharacteristic, truthState_N_new(k+9, :), '--m')
    end

    figure(204)
    hold on
    for k = 1:3
        subplot(3, 1, k)
        hold on
        % grid minor
        plot((time-time(1))*tCharacteristic, (truthState_N_new(k, :)-truthState_N_new(k+9, :)), '--r')
    end


    figure(205)
    hold on 
    for k = 1:3
        subplot(3, 1, k)
        hold on
        % grid minor
        plot((time-time(1))*tCharacteristic, (truthState_N_new(k, :)-truthState_N(k, :)), '--r')
    end

    figure(206)
    hold on 
    for k = 1:3
        subplot(3, 1, k)
        hold on
        % grid minor
        plot((time-time(1))*tCharacteristic, (truthState_N_new(k+3, :)-truthState_N(k+3, :)), '--r')
    end
end


%% Convert outputs to the expected names for the rest of the sim
truthState_N_new = truthState_R_new; % return nondim and rotating states
RBody2_N = RBody2_R;
RBody1_N = RBody1_R;
% Leaving frames defined as body to inertial -- may come back and want to
% consider trying to define body to rotating instead

% figure(206)
%     for k = 1:3
%         subplot(3, 1, k)
%         hold on
%         grid minor
%         plot(time, truthState_N_new(k+6, :))
%     end
%     title('Acceleration')

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
vc = xc(4:6);
rd = xd(1:3);
vd = xd(4:6);

yaxis = rc - rd;
yaxis = yaxis/norm(yaxis);

zaxis = vc - vd;
zaxis = zaxis/norm(zaxis);
zaxis = cross(zaxis, yaxis);
zaxis = zaxis/norm(zaxis);

% rdUnit_N = rd/norm(rd);
% zaxis = cross(rdUnit_N, yaxis);
% zaxis = zaxis/norm(zaxis);

% zaxis = vd/norm(vd);
% zaxis = cross(zaxis, yaxis);
% zaxis = zaxis/norm(zaxis);


xaxis = cross(yaxis, zaxis);
xaxis = xaxis/norm(xaxis);

if norm(pointNoise) == 0
    noise = eye(3);
else
    noise = pointNoise*randn(3, 1);
    noise = PRV2C(noise);
end
RBody_N = [xaxis, yaxis, zaxis]'*noise;

end