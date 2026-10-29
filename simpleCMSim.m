%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple Simulation Truth and Measurment Gen
%
%
% Created: 09/07/22
% Last Modified: 11/14/22
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% house keeping
if ~exist('flagSetFlags', 'var')
    clear
    close all
    flagScenario = 'Erros';
    fileBase = fullfile('.', 'MAT', flagScenario);
    measNum = 5;
    flagGravModelTruth = 'j2On';                      % Sets label for saving data with the truth data model
    plotCheck = true;
    beta = [0, 1/3600]; % time constan -- 1 hour
    bias = [0; 0; 0; 1e-9; 1e-9; 1e-9]; % initial bias
    gate = 1e-7;
    eom = @sphericalHarmonics_accelerometerGating;
end


%% Initialization

% Initial Chief Spacecraft Keplarian Elements
ac = 40e3; % m
ec = 0.01;
incc = 95*pi/180;  % rad
Omegac = 0; % rad
omegac = 0; % rad
Mc = 0; % rad
uc = Mc + omegac;


% Initial Deputy Spacecraft Keplarian Elements

doe = [0; 5; 0; 2; 0; 2]*1e3/ac;
ad = doe(1)*ac + ac; % m
incd = doe(5) + incc; % rad
Omegad = doe(6)/sin(incc) + Omegac; %rad
ud = doe(2) + uc - (Omegad - Omegac)*cos(incc);

b = doe(3:4) + [ec*cos(omegac); ec*sin(omegac)];
omegad = atan2(b(2), b(1)); % rad
Md = ud - omegad; % rad
ed = b(2)/cos(omegad); %


% Convert to inertial elements
% 433 Eros Info
mu = 4.4627e5; % m^3/s^2 -- From: A Global Solution for the Gravity Field, Rotation, Landmarks, and Ephemeris of Eros
J2 = 0.052478; % From: Determination of Shape, Gravity, and Rotational State of Asteroid 433 Eros
J3 = -0.0014; 
J4 = -0.012134;
J5 = 0;
J6 = 0;
J7 = 0; 
J8 = 0;
RBody = 16e3; % m % From: Determination of Shape, Gravity, and Rotational State of Asteroid 433 Eros

if strcmp(flagGravModelTruth, '8x8On')
    gravParamTruth = [mu, RBody, J2, J3, J4, J5, J6, J7, J8, 1];
elseif strcmp(flagGravModelTruth, 'j2On')
    gravParamTruth = [mu, RBody, J2, 0, 0, 0, 0, 0, 0, 0];
elseif strcmp(flagGravModelTruth, 'j2j3On')
    gravParamTruth = [mu, RBody, J2, J3, 0, 0, 0, 0, 0, 0];
elseif strcmp(flagGravModelTruth, 'j2j3j4On')
    gravParamTruth = [mu, RBody, J2, J3, J4, 0, 0, 0, 0, 0];
else
    fprintf('Only have J2-J4 Saved, Propagating with J2-J4 On\n')
    gravParamTruth = [mu, RBody, J2, J3, J4, 0, 0, 0, 0, 0];
end

P = 2*pi*sqrt(ac^3/mu);
dt = 2.5*3600/30;
N = ceil(4*P/dt);

truthState_N = zeros(18, N);
time = 0:dt:P*4;

inputFlag = "kep";
outputFlag = "cartesian";

[ xc ] = orbitalElmConverter( inputFlag, [ac; ec; incc*180/pi; Omegac*180/pi; omegac*180/pi; Mc*180/pi], 0, outputFlag, mu );
[ xd ] = orbitalElmConverter( inputFlag, [ad; ed; incd*180/pi; Omegad*180/pi; omegad*180/pi; Md*180/pi], 0, outputFlag, mu );

opts = odeset('RelTol',1e-13,'AbsTol',1e-13);
burnStart = round(3.2*P/dt); % time step 803
burnEnd = burnStart + 15*60/dt; % time step 806
burn = [0; 720e-6];
burnMagSigma = 0.15*burn;
burnDirSigma = 0.5; % deg




burnHistCommand = zeros(6, N);
burnHistActual = zeros(6, N);
accelerometer = zeros(6, N);


calcBurn = false;

initialState = [xc; bias(1:3); xd; bias(4:6)];
truthState_N(:, 1) = initialState';
iminus1 = 1;

omegaErros = 360 /(5.27025547*3600); % deg/s
lam0 = 0;
RBody1_N = NaN(3, 3, N);
RBody2_N = NaN(3, 3, N);
RCF1_Body1_N = NaN(3, 3, N);
RCF2_Body1_N = NaN(3, 3, N);
RCF1_Body2_N = NaN(3, 3, N);
RCF2_Body2_N = NaN(3, 3, N);

RN_rot = NaN(3, 3, N);

RN_rot(:, :, 1) = eye(3);

pointNoise = eye(3)*(1e-6);
pointNoise_chol = chol(pointNoise);




if dt > 1
   dtProp = 1;
else
    dtProp = dt;
end

RBody1_N(:, :, 1) = cameraPoint(truthState_N(1:6, 1), truthState_N(10:15, 1), pointNoise);
RBody2_N(:, :, 1) = cameraPoint(truthState_N(10:15, 1), truthState_N(1:6, 1), pointNoise);
for i = 2:N
    RN_rot(:, :, i) = rotz(lam0 + omegaErros*time(i));
    RBody1_N(:, :, i) = cameraPoint(truthState_N(1:6, iminus1), truthState_N(10:15, iminus1), pointNoise);
    RBody2_N(:, :, i) = cameraPoint(truthState_N(10:15, iminus1), truthState_N(1:6, iminus1), pointNoise);
    
    timeInt = time(iminus1:i);
    if i > burnStart && i < burnEnd
        if ~calcBurn
            hVec = cross(initialState(1:3), initialState(4:6));
            hUnit = hVec/norm(hVec);
            burnUnitNoisy = -rotz(randn(1)*burnDirSigma)*roty(randn(1)*burnDirSigma)*rotx(randn(1)*burnDirSigma)*hUnit;
            burnVecCommand = -burn(2)*hUnit';
            burnVecActual = burn(2)*burnUnitNoisy*randn(1)*burnMagSigma(2);
            burnVecCommand = [zeros(3, 1); burnVecCommand'];
            burnVecActual = [zeros(3, 1); burnVecActual] + burnVecCommand;
            calcBurn = true;
        end
        burnHistActual(:, i) = burnVecActual;
        burnHistCommand(:, i) = burnVecCommand;
        timeProp = timeInt(1):dtProp:timeInt(2);
        Nprop = length(timeProp)-1;
        initialState = truthState_N(:, iminus1);
        accelerometer(:, i) = burnVecActual + initialState([7:9, 16:18]);
        Rbody_N = zeros(3, 3, 2);
        Rbody_N(:, :, 1) = RBody1_N(:, :, i);
        Rbody_N(:, :, 2) = RBody2_N(:, :, i);
        
        for j = 1:Nprop
            state = rk4(initialState, timeProp(j:j+1), eom, gravParamTruth, beta, accelerometer(:, i), gate, Rbody_N, RN_rot);
            initialState = state;
        end
        truthState_N(:, i) = state;
    else
        timeProp = timeInt(1):dtProp:timeInt(2);
        Nprop = length(timeProp)-1;
        initialState = truthState_N(:, iminus1);
        Rbody_N = zeros(3, 3, 2);
        Rbody_N(:, :, 1) = RBody1_N(:, :, i);
        Rbody_N(:, :, 2) = RBody2_N(:, :, i);
        for j = 1:Nprop
            state = rk4(initialState, timeProp(j:j+1), eom, gravParamTruth, beta, accelerometer(:, i), gate, Rbody_N, RN_rot);
            initialState = state;
        end
        truthState_N(:, i) = state;
    end
    RBody1_N(:, :, i) = cameraPoint(truthState_N(1:6, i), truthState_N(10:15, i), pointNoise);
    RBody2_N(:, :, i) = cameraPoint(truthState_N(10:15, i), truthState_N(1:6, i), pointNoise);
    iminus1 = i;
end


truthState = truthState_N;
truthState(1:9, :) = truthState(1:9, :) - truthState(10:18, :);

%% Verify Traj
if plotCheck
    figure
    plot3(truthState_N(1, :), truthState_N(2, :), truthState_N(3, :))
    hold on
    grid minor
    plot3(truthState_N(10, :), truthState_N(11, :), truthState_N(12, :))
end

outputFlag = "kep";
inputFlag = "cartesian";

[ oec_1 ] = orbitalElmConverter( inputFlag, truthState(1:6, burnStart-1)', 0, outputFlag, mu );
[ oed_1 ] = orbitalElmConverter( inputFlag, truthState(10:15, burnStart-1)', 0, outputFlag, mu );


[ oec ] = orbitalElmConverter( inputFlag, truthState(1:6, burnEnd+1)', 0, outputFlag, mu );
[ oed ] = orbitalElmConverter( inputFlag, truthState(10:15, burnEnd+1)', 0, outputFlag, mu );

acdiy = angdiff(oed(4)*pi/180, oec(4)*pi/180)*sind(oec(3))*oec(1);



%% Generate Measuremnts
focalLength = 7.6e-3; % m
fov_rad = [44; 32]*pi/180; % rad
numPixel = [2592; 1944]; % -

R_opNav = fov_rad./numPixel/2;
if measNum == 6
    R = diag([0.01, 0.001, R_opNav', R_opNav'].^2);
    R_rangeRangeRate_chol = chol(R(1:2, 1:2), 'lower');
    R_rangeRangeRate = R(1:2, 1:2);
    R_opNav_chol = chol(R(3:4, 3:4), 'lower');
elseif measNum == 4
    R = diag([0.01, 0.001, R_opNav'].^2);
    R_rangeRangeRate_chol = chol(R(1:2, 1:2), 'lower');
    R_rangeRangeRate = R(1:2, 1:2);
    R_opNav_chol = chol(R(3:4, 3:4), 'lower');
elseif measNum == 3
    % range, relative bearing 
    R = diag([0.01, R_opNav'].^2);
    R_rangeRangeRate_chol = sqrt(R(1, 1));
    R_rangeRangeRate = R(1, 1);
    R_opNav_chol = chol(R(2:3, 2:3), 'lower');
elseif measNum == 5
    % range, relative bearing 
    R = diag([0.01, R_opNav', R_opNav'].^2);
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
% load('featureSet.mat')
% [L_set] = genFeatures(RBody, N_feature);
L_set = [0; 0; 0];

% Assumes boresight is in the z direction of the camera frame
RCF1_Body1 = rotx(-90); % looking at body 2
RCF2_Body1 = roty(90); % looking toward the central body (ish)
RCF1_Body2 = rotx(-90); % looking at body 1
RCF2_Body2 = roty(90); % looking toward the central body (ish)


H = NaN(measNum, 18, N);

featureIndex = NaN(1, N);


for i = 1:N
    RCF1_Body1_N(:, :, i) = RCF1_Body1*RBody1_N(:, :, i);
    RCF2_Body1_N(:, :, i) = RCF2_Body1*RBody1_N(:, :, i);
    RCF1_Body2_N(:, :, i) = RCF1_Body2*RBody2_N(:, :, i);
    RCF2_Body2_N(:, :, i) = RCF2_Body2*RBody2_N(:, :, i); 
    
    if measNum == 4
        [meas_i, H_i] = rngRngRate_opTarget_meas_H(truthState(:, i), RCF1_Body2_N(:, :, i), true);
        meas_truth(:, i) = meas_i;
        meas(1:2, i) = meas_i(1:2) + R_rangeRangeRate_chol*randn(2, 1);
        meas(3:4, i) = meas_i(3:4) + R_opNav_chol*randn(2, 1);
        H(:, :, i) = H_i;
    elseif measNum == 5
        j = 1;
        while j < N_feature
            feature_i = inFOV(squeeze(RN_rot(:, :, i)), L_set(:, j), truthState(:, i), squeeze(RCF2_Body2_N(:, :, i)), fov_rad);
            if  ~isnan(feature_i(1))
                featureIndex(i) = j;
                RCF2_Body2_N(:, :, i) = cameraPoint_atFeature(truthState(10:18, i), pointNoise, L_set(:, j), RN_rot(:, :, i));
                [meas_i, H_i] = relNavRngRngRate_opTarget_feature_meas_H(truthState(:, i), RCF1_Body2_N(:, :, i), RCF2_Body2_N(:, :, i), RN_rot(:, :, i), L_set(:, 1), true);
                meas_truth(:, i) = meas_i;
                meas(1, i) = meas_i(1) + R_rangeRangeRate_chol*randn(1, 1);
                meas(2:3, i) = meas_i(2:3) + R_opNav_chol*randn(2, 1);
                meas(4:5, i) = meas_i(4:5) + R_opNav_chol*randn(2, 1);
                H(:, :, i) = H_i;
                j = N_feature;
            end
            j = j + 1;
        end
        if isnan(meas_truth(:, i))
            [meas_i, H_i] = relNavRngRngRate_opTarget_feature_meas_H(truthState(:, i), RCF1_Body2_N(:, :, i), RCF2_Body2_N(:, :, i), RN_rot(:, :, i), L_set(:, 1), true);
            meas_truth(:, i) = meas_i;
            meas(1, i) = meas_i(1) + R_rangeRangeRate_chol*randn(1, 1);
            meas(2:3, i) = meas_i(2:3) + R_opNav_chol*randn(2, 1);
            meas(4:5, i) = meas_i(4:5) + R_opNav_chol*randn(2, 1);
            H(:, :, i) = H_i;
        end
    elseif measNum == 3
        if strcmp(flagScenario, 'ErrosRel')
            [meas_i, H_i] = relNavRng_opTarget_meas_H(truthState(:, i), RCF1_Body2_N(:, :, i), true);
            meas_truth(:, i) = meas_i;
            meas(1, i) = meas_i(1) + R_rangeRangeRate_chol*randn(1, 1);
            meas(2:3, i) = meas_i(2:3) + R_opNav_chol*randn(2, 1);
            H(:, :, i) = H_i;
        else 
            [meas_i, H_i] = rng_opTarget_meas_H(truthState_N(:, i), RCF1_Body2_N(:, :, i), true);
            meas_truth(:, i) = meas_i;
            meas(1, i) = meas_i(1) + R_rangeRangeRate_chol*randn(1, 1);
            meas(2:3, i) = meas_i(2:3) + R_opNav_chol*randn(2, 1);
            H(:, :, i) = H_i;
        end
    end
    

end



%% Save Off Info
if ~exist(fileBase, 'dir')
    mkdir(fileBase)
end

if ~contains(flagScenario, 'Rel')
    truthState = truthState_N;
end

save(fullfile(fileBase, 'truthState.mat'), 'truthState')
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

noise = pointNoise*randn(3, 1);
noise = PRV2C(noise);

RBody_N = [xaxis, yaxis, zaxis]'; % *noise;

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

noise = pointNoise*randn(3, 1);
noise = PRV2C(noise);

RCF_N = [xaxis, yaxis, zaxis]'; % *noise;

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