function [ri, PriMinus, H] = GMM_EKFmeas_prefitOutput(x, P, meas, R, measFunc, L, RCF1_N_i, RCF2_N_i, RN_rot_i, n_elm, measFlag, flagUseTrueJacobians, H_true)
[n, n_nodes] = size(x);
nonNaN = ~isnan(meas);
if strcmp(measFlag, 'relRngRelAzEl') % This case works 
    % if nargin > 19
    %     [yh, H] = measFunc(x, RCF1_N_i, n_elm);
    %     % H(2:3, :) = H(2:3, :)/L_char;
    % else
        [yh, H] = measFunc(x, RCF1_N_i, n_elm);
    % end
    m = 3;
elseif strcmp(measFlag, 'relRngRelAzElFeatAzEl') % have NOT touched this
    [yh, H] = measFunc(x, RCF1_N_i, RCF2_N_i, eye(3), L, n_elm);
    m = 5;
elseif strcmp(measFlag, 'direct')% have NOT touched this
    [yh, H] = measFunc(x);
    m = 3;
else % have NOT touched this
    [yh, H] = measFunc(x, RCF1_N_i, RCF2_N_i, RN_rot_i, L, n_elm);
end

if flagUseTrueJacobians
    H = repmat(H_true, 1, 1, n_nodes);
end

H = H(nonNaN, :, :); % correct shape
Ri = R(nonNaN, nonNaN);
if strcmp(measFlag, 'direct')
    ri = meas-yh;
else
    ri = diffMeas(meas, yh, nonNaN);
end

Pxy = pagemtimes(P,'none', H, 'transpose');
HPHT = pagemtimes(H, Pxy);
Pyy = (HPHT + Ri);
PriMinus(nonNaN, nonNaN, :) = Pyy(nonNaN, nonNaN, :); % correct shape







ri = ri(nonNaN, :, :);
ri = reshape(ri, [], n_nodes); % correct shape




end



function [difference] = diffMeas(yi, yh, nonNan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differences measurments when angular measurments are used
% 
% Inputs: 
%   yi      recieved measurment 
%   yh      estimated measurment 
%   nonNan  measurment values that are not NaN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = length(nonNan);
[~, ~, n_nodes] = size(yh);
difference = NaN(m, 1, n_nodes);
difference(1, :, :) = yi(1, :, :) - yh(1, :, :);
for i = 2:m
    if nonNan(i)
        for j = 1:n_nodes
            difference(i, :, j) = angdiff(yh(i, :, j), yi(i));
        end
    end
end
difference = difference(nonNan, :, :);
end
