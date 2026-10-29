function [x, P, w_nodesPlus, riPlus, ri, Pri, PriMinus, H] = GMM_EKFmeas(x, P, w_nodes, meas, R, measFunc, L, RCF1_N_i, RCF2_N_i, RN_rot_i, considerParam, n_elm, measFlag, underweight, accelerometer, gate, accelNoise, flagUnderweight, flagUseTrueJacobians, H_true, flagInteracting)
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

if n_elm(2) == 9 
    [rowH, ~] = size(H);
    if norm(accelerometer) < gate
        % When the acceleration is below the threshold to sense a maneuver
        H = [H; zeros(3, n-3, n_nodes), repmat(eye(3), 1, 1, n_nodes)]; % Add in accelerometer to measurements
        ri = [ri; reshape(accelerometer - x((n_elm(1)+7):(n_elm(1) + 9), :), [], 1, n_nodes)]; % Add in measurment to the measurements
        Ri = blkdiag(Ri, accelNoise);
        m = m+3;
        appendFlag = true;
        % ri = ri*0;
    else
        appendFlag = false;
        % fprintf('Chaser burning, no accelerometer bias update\n')
    end
end

Pxy = pagemtimes(P,'none', H, 'transpose');
HPHT = pagemtimes(H, Pxy);

if strcmp(flagUnderweight, 'variable') % have NOT touched this
    if trace(P(1:6, 1:6)) > 1e-3 % have NOT touched this
        ck = 3/yh(1)^2; % m -- assumes that the first measurement it range
        HPHT_rel = H(1:rowH, 1:6)*P(1:6, 1:6)*H(1:rowH, 1:6)';
        underweight = ck * trace(P(1:6, 1:6))^2/(2 * trace(HPHT_rel));
    else
        underweight = 0;
    end
elseif strcmp(flagUnderweight, 'switchOff')% have NOT touched this
    if trace(P(1:6, 1:6)) < 1
        underweight = 0;
    end
end


B = underweight*HPHT;
Pyy = (HPHT + Ri);
Pyyinv = pagemldivide(Pyy, eye(m));
Ki = pagemtimes(Pxy, pagemldivide((Pyy + B), eye(m)));
Ki(considerParam, :, :) = 0;
x = x + reshape(pagemtimes(Ki,ri), n, n_nodes); % correct shape
riPlusTemp = pagemtimes(Ri, pagemtimes(Pyyinv, ri));
riPlus(nonNaN, :, :) = riPlusTemp(nonNaN, :, :);

ImKH = eye(n) - pagemtimes(Ki,H);
P = pagemtimes(ImKH, pagemtimes(P, 'none', ImKH, 'transpose')) + pagemtimes(Ki, pagemtimes(Ri,'none', Ki, 'transpose')); % correct shape
PriMinus(nonNaN, nonNaN, :) = Pyy(nonNaN, nonNaN, :); % correct shape
PriTemp = pagemtimes(Ri, pagemtimes(Pyyinv,'none', Ri, 'transpose')); % correct shape
Pri(nonNaN, nonNaN, :) = PriTemp(nonNaN, nonNaN, :);



% Calculating New Weights 
detPyy = prod(pageeig(Pyy(1:rowH, 1:rowH, :), 'vector'));
Pyyinv_mxm = pagemldivide(Pyy(1:rowH, 1:rowH, :), eye(rowH));
qyi =  squeeze(1./sqrt(2*pi*detPyy).*exp(-1/2*pagemtimes(ri(1:rowH, :, :), 'transpose', pagemtimes(Pyyinv_mxm, ri(1:rowH, :, :)), 'none'))); % mvnpdf(meas, yh, Pyy);




if n_nodes > 1 && any(isinf(qyi))
    djwiodhwioq = 0;
    qyi(isinf(qyi)) = 0;
end

if ~flagInteracting
    w_nodesPlus = w_nodes.*qyi;
else
    w_nodesPlus = qyi;
end
if n_nodes == 1
    w_nodesPlus = 1;
else
    if sum(w_nodesPlus) == 0
        fprintf('No weight, this is bad!\n')
    end
    w_nodesPlus = w_nodesPlus/sum(w_nodesPlus);

end


ri = ri(nonNaN, :, :);
riPlus = reshape(riPlus, [], n_nodes); % correct shape
ri = reshape(ri, [], n_nodes); % correct shape


% x, P, riPlus, ri, Pri, PriMinus, H

% if appendFlag
%     H = H(1:rowH, :);
% end

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
