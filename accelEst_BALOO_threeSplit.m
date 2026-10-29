function [x_2, P_2, w_2] = accelEst_BALOO_threeSplit(kk, x, P, w, ri, y, time, bufferSize, RCF1_Body2_N, eom, gravParam, beta_markov, accelerometer, gate, Rbody_N, RN_rot,relNavFlag, accelBias, n_elm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Takes in information on a node to be split into two and produces
% the two new nodes (1 unchnaged other than its weight, one new node with
% an updated unmodeled target acceleration)
% 
% Inputs:
%   x       Time history of the state of the node to be split 
%   P       Time history of the covaraince of the nodes to be split
%   w       Current weight of the nodes to be split 
%   ri      Prefit measurment residuals over a window in the previous time
%           history
% 
% Outputs:
%   x_2     New nodes' time history of the state (NaN for new node until
%           current time)
%   P_2     New nodes' time history of the state covariance (NaN for new 
%           node until current time)
%   w_2     New nodes' weights 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eom, gravParam, beta_markov, accelerometer, gate, Rbody_N, RN_rot,relNavFlag, accelBias, n_elm,
[n, ~, ~] = size(x);
x_2 = NaN(n, 3);
P_2 = NaN(n, n, 3);
w_2 = zeros(3, 1);


x_2(:, 1) = x(:, 1, end);
P_2(:, :, 1) = P(:, :, 1, end);
w_2(1) = w(1);

x_2(:, 2) = x(:, :, end);
P_2(:, :, 2) = P(:, :, :, end);
w_2(2) = w(1);

x_2(:, 3) = x(:, :, end);
P_2(:, :, 3) = P(:, :, :, end);
w_2(3) = w(1);

[mbp1, ~] = size(ri);
m = mbp1/(bufferSize);

y = y(~isnan(ri)); % remove NaN values if present 
ri = ri(~isnan(ri)); % remove NaN values if present 


[mbp1, ~] = size(ri);
bufferSizep1 = mbp1/m;



measEstRhoVec_CF = zeros(3, 1, bufferSizep1);
y = reshape(y, m, 1, bufferSizep1);


% Create rho vec in camera frame 
cosAtilde2 = cos(y(2, :, :)).^2;
tanBtilde = tan(y(3, :, :));
denominator = (1 + pagemtimes(cosAtilde2, tanBtilde.^2));
numerator = pagemtimes(y(1, :, :).^2, cosAtilde2);
measEstRhoVec_CF(3, :, :) = sqrt(pagemtimes(numerator, 1  ./ denominator));
measEstRhoVec_CF(1, :, :) = pagemtimes(measEstRhoVec_CF(3, :, :), tan(y(2, :, :)));
measEstRhoVec_CF(2, :, :) = pagemtimes(measEstRhoVec_CF(3, :, :), tanBtilde);

x_flip = flip(x, 3);
RCF1_Body2_N_flip = flip(RCF1_Body2_N, 3);

% Convert rho vec from measurments to inertial frame 
measEstRhoVec_N = pagemtimes(RCF1_Body2_N_flip(:, :, 1:bufferSizep1), 'transpose', measEstRhoVec_CF, 'none');

% Difference of actual and estimated 
diffMeasEst_state = measEstRhoVec_N - x_flip(1:3, :, 1:bufferSizep1);


% Weight more recent measurments more
b = 2;
s0 = (b-1)/(b^(bufferSizep1) - 1);
w_measAcc_favorNewer = reshape(flip(s0*b.^(0:(bufferSizep1-1))), 1, 1, bufferSizep1);
% w_measAcc = w_measAcc_favorNewer;


% Equally weighted 
w_measAcc =  ones(1, 1, bufferSizep1)/bufferSizep1;

diffMeasEst_state_weightedMeas = pagemtimes(diffMeasEst_state, w_measAcc);


dt = diff(time);
diffMeasEst_state_weightedMeas = pagemtimes(diffMeasEst_state_weightedMeas,  reshape(1 ./ (dt(1:bufferSizep1).^2), 1, 1, []));
diffMeasEst_state_weightedMeas = 2*sum(diffMeasEst_state_weightedMeas, 3);

aEst_1 = diffMeasEst_state_weightedMeas;

% 0.063423146260032

% diffMeasEst_state = reshape(diffMeasEst_state, [], 1);
% aEst = (timeMat'*timeMat)\eye(3)*timeMat'*diffMeasEst_state;

% aEst = (timeMat'*timeMat)\eye(3)*timeMat'*diffMeasEst_state_weightedMeas;




dt_0 = (time - time(1)).^2;
dt_0 = flip(dt_0(2:end-1));
dt_0 = dt_0(1:bufferSizep1);
dt = ones(size(dt_0));
w2divdt2_0 = w_measAcc_favorNewer.*reshape(2./dt_0, 1, 1, []);
% w2divdt2_0 = w_measAcc_favorNewer.*reshape(2./dt_0, 1, 1, []);
w2divdt2_0 = pagemtimes(w2divdt2_0,  repmat(eye(3), 1, 1, bufferSizep1));
w2divdt2_0 = reshape(w2divdt2_0, 3, []);


daEst = w2divdt2_0*reshape(diffMeasEst_state, [], 1);


% What if we use the last half of measurments to get accel est?
numBack = round(bufferSizep1/2);
dt_0_new = (time - time(1)).^2;
dt_0_new = flip(dt_0_new(2:numBack+1));
dt = ones(size(dt_0_new));

b = 2;
s0 = (b-1)/(b^(numBack) - 1);
w_measAcc_favorNewer = reshape(flip(s0*b.^(0:(numBack-1))), 1, 1, numBack);

w_measAcc =  ones(1, 1, numBack)/numBack;
w2divdt2_0 = w_measAcc.*reshape(2./dt, 1, 1, []);
w2divdt2_0 = pagemtimes(w2divdt2_0,  repmat(eye(3), 1, 1, numBack));
w2divdt2_0 = reshape(w2divdt2_0, 3, []);

daEst_new = w2divdt2_0*reshape(diffMeasEst_state(3*(bufferSizep1-numBack)+1:3*bufferSizep1), [], 1);

% Apply maneuver half way back 
% drEst = 1/2*daEst*(dt_0(round(bufferSizep1/2)));
% dvEst = daEst*sqrt(dt_0(round(bufferSizep1/2)));
% 
% drEst_new = 1/2*daEst_new*(dt_0_new(round(bufferSizep1/2)));
% dvEst_new = daEst_new*sqrt(dt_0_new(round(bufferSizep1/2)));

% % Apply maneuver at end of buffer
drEst = 1/2*daEst*(dt_0(end));
dvEst = daEst*sqrt(dt_0(end));

drEst_new = 1/2*daEst_new*(dt_0(end));
dvEst_new = daEst_new*sqrt(dt_0(end));

% % Apply maneuver at beginning of buffer
% drEst = 1/2*daEst*(dt_0(1));
% dvEst = daEst*sqrt(dt_0(1));
% 
% drEst_new = 1/2*daEst_new*(dt_0_new(1));
% dvEst_new = daEst_new*sqrt(dt_0_new(1));


%% TODO: Try taking the estimate of acceleration, then taking xminus from 
% the state buffer and using the full nonlinear prop to get the new 
% relative target state 


% drEst = 1/2*daEst*(dt_0(1));
% dvEst = daEst*sqrt(dt_0(1));



% 0; -0.054926055849122; -0.031711573130019

% 0; -0.1099; -0.0634




a_0 = x_2(7:9, 2) + daEst;


if relNavFlag 
    % If Conversion needed
    xr = x(1:n_elm(1), :, end-numBack);
    xc = x(n_elm(1) + 1:end, :, end-numBack);
    xt = xr;
    xt(1:6, :, :) = xc(1:6, :) + xr(1:6, :);
    initialState = [xt; xc];
    % M(1:6, (n_elm(1)+(1:6))) = -eye(6);
else
    initialState = x;
end

initialState_zeroBurn = initialState;

initialState(7:9) = a_0;
initialState_zeroBurn(7:9) = zeros(3, 1);


N = length(time);

for j = (N-numBack):N-1
    times = time(j:j+1);

    state = rk4(initialState, times, eom, gravParam, beta_markov, accelerometer, gate, Rbody_N, RN_rot, n_elm, accelBias);
    initialState = state;

    state_zeroBurn = rk4(initialState_zeroBurn, times, eom, gravParam, beta_markov, accelerometer, gate, Rbody_N, RN_rot, n_elm, accelBias);
    initialState_zeroBurn = state_zeroBurn;

end

if relNavFlag
    xt = state(1:n_elm(1));
    xc = state(n_elm(1)+1:end);
    xr = xt;
    xr(1:6) = xt(1:6) - xc(1:6);
    
    xminus = [xr; xc];

    xt = state_zeroBurn(1:n_elm(1));
    xc = state_zeroBurn(n_elm(1)+1:end);
    xr = xt;
    xr(1:6) = xt(1:6) - xc(1:6);
    
    xminus_zeroState = [xr; xc];
else
    xminus = state;
    xminus_zeroState = state_zeroBurn;
end



% lifted from AGEIS algorithm from Orbital Debris
% % alphaTilde3 = [0.709303325662705,   0.290696674337295];
% alphaT
% ilde3 = [0.9,   0.1];
% alphaTilde3 = [0.2,   0.8];
alphaTilde3 = [1/3,   1/3, 1/3];
% mTilde3 = [0, 1.0575154615];
% sigmaTilde3 = [1.007349433050000, 1.007349433050000];
sigmaTilde3 = [1.00, 1.00, 1.00];


% x_2(7:9, 2) = x_2(7:9, 2) + aEst_1;
% x_2(7:9, 2) = x_2(7:9, 2) + daEst;
% x_2(1:9, 2) = x_2(1:9, 2) + [drEst; dvEst; daEst];
x_2(1:9, 2) = xminus(1:9);
x_2(1:9, 3) = xminus_zeroState(1:9);

% x_2(1:9, 2) = x_2(1:9, 2) + [drEst_new; dvEst_new; daEst_new];
w_2 = alphaTilde3'.*w_2;
P_2 = sigmaTilde3(1)*P_2;



end 