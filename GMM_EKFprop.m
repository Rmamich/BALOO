function [xminus, Pminus, F] = GMM_EKFprop(times, x, P, eom, Qhat, gravParam, beta_markov, accelerometer, gate, Rbody_N, RN_rot,relNavFlag, accelBias, n_elm, accelVRW, flagUseTrueJacobians, Phi_true, CRTBPFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Propagates the state x forward in inertial space using the given
% eom 
% 
% Inputs: 
%   times               times for propagation, expected length of 2 [s]
%   x                   input state:
%                           x(1:3)     target or relative position in the 
%                                       inertial frame [m] 
%                           x(4:6)     target or relative velocity in the 
%                                       inertial frame [m/s]
%                           x(7:9)     target or relative unmodeled 
%                                        acceleration [m/s^2]
%                                       (may not be present)
%                           x(10:12)   chaser position [m]
%                           x(13:15)   chaser velocity [m/s]
%                           x(16:18)   chaser unmodeled acceleration est 
%                                       [m/s^2] (may not be present)
%   P
%   eom
%   burn
%   Qhat
%   gravParam
%   beta_markov
%   frameConvertFlag    flag to determine if the first 1-(6 or 9)
%                           components are relative inertial or full
%                           inertial 
% 
% Outputs:
%   xminus
%   Pminus
%   Phi
% 
% Author: Rachel Mamich 
% Date Created: 10/24/2022
% Date Modified: 11/02/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = times(2) - times(1);
M = eye(size(x, 1)); % identity allocation

if relNavFlag 
    % If Conversion needed
    xr = x(1:n_elm(1), :);
    xc = x(n_elm(1) + 1:end, :);
    xt = xr;
    xt(1:6, :, :) = xc(1:6, :) + xr(1:6, :);
    initialState = [xt; xc];
    % M(1:6, (n_elm(1)+(1:6))) = -eye(6);
else
    initialState = x;
end


if nargin > 17
    [F_j] = buildPhi_page(initialState, gravParam, dt, beta_markov, n_elm,relNavFlag, CRTBPFlag, times(1));
    % buildPhi(      initialState_nonDimR, gravParamTruth, dt, beta_nonDim, [9, 9], true, true, timeProp(j));
    F = F_j;
else
    [F_j] = buildPhi_page(initialState, gravParam, dt, beta_markov, n_elm,relNavFlag);
    % [F_j_check] = buildPhi(initialState, gravParam, dt, beta_markov, n_elm,relNavFlag);
    F = F_j;
end

if flagUseTrueJacobians
    F = Phi_true;
end


Pminus = pagemtimes(F, pagemtimes(P, 'none', F, 'transpose'));

% TODO: Add this back in later
% if norm(accelerometer) > gate
%     QhatAccel = continuious2discreteProcNoise(dt, beta_markov, accelVRW, false, false, false); % 12x12
%     if relNavFlag
%         Qhat(1:6, 1:6) = Qhat(1:6, 1:6) + QhatAccel(1:6, 1:6);
%     end
%     if n_elm(1) > 6
%         Qhat(10:15, 10:15) = Qhat(10:15, 10:15) + QhatAccel(7:12, 7:12);
%     else
%         Qhat(7:12, 7:12) = Qhat(7:12, 7:12) + QhatAccel(7:12, 7:12);
%     end
% end

Pminus = Pminus + pagemtimes(M,pagemtimes(Qhat, 'none', M, 'transpose'));


state = rk4(initialState, times, eom, gravParam, beta_markov, accelerometer, gate, Rbody_N, RN_rot, n_elm, accelBias);

if relNavFlag
    xt = state(1:n_elm(1), :);
    xc = state(n_elm(1)+1:end, :);
    xr = xt;
    xr(1:6, :) = xt(1:6, :) - xc(1:6, :);
    
    xminus = [xr; xc];
else
    xminus = state;
end



end

