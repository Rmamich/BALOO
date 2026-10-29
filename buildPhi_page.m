function [F, A] = buildPhi_page(x, gravParam, dt, beta, n_elm, relNavFlag, CRTBPFlag, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%
% Inputs:
%   x           x(1:n_elm(1)) target states in inertial space
%               x(n_elm(1)+1, end) chaser states in inertial space
%   gravParam   gravitational parameter array for central body
%   dt          delte time [s]
%   beta        time constants for gauss markov processes [1/s]
%   n_elm       number of elements for target and chaser states 
%   relNavFlag  true when the target states in the stm are for the relative
%                   inertial target 
%
% Outputs:
%   F           discrete time stm 
%   A           continuious time jacobian of dxdot/ dx
% 
% Author: Rachel Mamich
% Date Modified: 04/26/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n, n_bank] = size(x);
xt = x(1:n_elm(1), :);
xc = x((1+n_elm(1)):n_elm(1)+n_elm(2), :);
eye_page = repmat(eye(3), 1, 1, n_bank);

if nargin > 6
    if ~CRTBPFlag
        % Normal TBP + whatever
        % Chaser
        Ac = zeros(n_elm(2), n_elm(2), n_bank);
        Ac(1:3, 4:6, :) = eye_page;
        [dadr] = gravJac_page(xc(1:3, :), gravParam);
        Ac(4:6, 1:3, :) = dadr;
        if n_elm(2) == 9
            %     Ac(4:6, 7:9) = eye(3);
            Ac(7:9, 7:9, :) = -eye_page*beta(2); 
        end

        % Target
        At = zeros(n_elm(1), n_elm(1), n_bank);
        At(1:3, 4:6, :) = eye_page;
        [dadr] = gravJac_page(xt(1:3, :), gravParam);
        At(4:6, 1:3, :) = dadr;
        if n_elm(1) == 9
            At(4:6, 7:9, :) = eye_page;
            At(7:9, 7:9, :) = -eye_page*beta(1); % Assumes constant burn
        end
    else
        % CRTBP -- Assumes that the state passes in is in nondim,
        % rotational frame, also assumes that beta is nondim 
        % Chaser
        Ac = zeros(n_elm(2), n_elm(2), n_bank);
        Ac(1:3, 4:6, :) = eye_page;
        [dadr] = gravJac_CRTBP_page(xc(1:3, :), gravParam);
        Ac(4:6, 1:3, :) = dadr;
        Ac(4:6, 4:6, :) = repmat([0, 2, 0; -2, 0, 0; 0, 0, 0], 1, 1, n_bank);
        if n_elm(2) == 9
            %     Ac(4:6, 7:9) = eye(3);
            Ac(7:9, 7:9, :) = -eye_page*beta(2); % Assumes that beta is nondim
        end

        % dX_nonDimR_dX_dimN_c = dX_nonDimRot_dX_dimN(xc, t, gravParam);
        % % Ac is now wrt dimensional and inertial frame
        % Ac = dX_nonDimR_dX_dimN_c*Ac*dX_nonDimR_dX_dimN_c';
        % Target

        At = zeros(n_elm(1), n_elm(1), n_bank);
        At(1:3, 4:6, :) = eye_page;
        [dadr] = gravJac_CRTBP_page(xt(1:3, :), gravParam);
        At(4:6, 1:3, :) = dadr;
        At(4:6, 4:6, :) = repmat([0, 2, 0; -2, 0, 0; 0, 0, 0], 1, 1, n_bank);
        if n_elm(1) == 9
            At(4:6, 7:9, :) = eye_page;
            At(7:9, 7:9, :) = -eye_page*beta(1); % Assumes that beta is nondim
        end

        % Transform A by doing:
        % dX_dimN/dX_nonDimRot * A_nonDimRot * dX_dimN/dX_nonDimRot'
        % dX_nonDimR_dX_dimN_t = dX_nonDimRot_dX_dimN(xt, t, gravParam);
        % % At is now wrt dimensional and inertial frame
        % At = dX_nonDimR_dX_dimN_t*At*dX_nonDimR_dX_dimN_t';

    end
else
    % Normal TBP + whatever
    % Chaser
    Ac = zeros(n_elm(2), n_elm(2), n_bank);
    Ac(1:3, 4:6, :) = eye_page;
    [dadr] = gravJac_page(xc(1:3, :), gravParam);
    Ac(4:6, 1:3, :) = dadr;
    if n_elm(2) == 9
        %     Ac(4:6, 7:9) = eye(3);
        Ac(7:9, 7:9, :) = -eye_page*beta(2);
    end

    % Target
    At = zeros(n_elm(1), n_elm(1), n_bank);
    At(1:3, 4:6, :) = eye_page;
    [dadr] = gravJac_page(xt(1:3, :), gravParam);
    At(4:6, 1:3, :) = dadr;
    if n_elm(1) == 9
        At(4:6, 7:9, :) = eye_page;
        At(7:9, 7:9, :) = -eye_page*beta(1); % Assumes constant burn
    end

end


if ~relNavFlag
    % Full inertial
    A = [At, zeros(n_elm(1), n_elm(2), n_bank); zeros(n_elm(2), n_elm(1), n_bank), Ac];
else
    % Rel Nav with different configs
    if n_elm(1) == 6 && n_elm(2) == 9
        A = [At, At-Ac(1:6, 1:6, :), zeros(n_elm(1), 3, n_bank); zeros(n_elm(2), n_elm(1), n_bank), Ac];
    elseif n_elm(1) == 6 && n_elm(2) == 6
        A = [At, At-Ac;... 
            zeros(6, 6, n_bank), Ac];
    elseif n_elm(1) == 9 && n_elm(2) == 9
        A = [At, [At(1:6, 1:6, :)-Ac(1:6, 1:6, :), zeros(6, 3, n_bank); zeros(3, 9, n_bank)];...
            zeros(9, 9, n_bank), Ac];
    elseif n_elm(1) == 9 && n_elm(2) == 6
        A = [At, [At(1:6, 1:6, :)-Ac(1:6, 1:6, :), zeros(6, 3, n_bank); zeros(3, 9, n_bank)]; zeros(n_elm(2), n_elm(1), n_bank), Ac];
    end
end


% Convert to discrete time 
% Eigen Decomp to get F = expm(A*dt); 
Adt = A*dt;
F = zeros(n,n,n_bank);
for j = 1:n_bank
   F(:, :, j) = expm(Adt(:, :, j)); 
end
% [V, lambda_ii] = pageeig(Adt, 'vector');
% expLam = exp(lambda_ii);
% F = pagemtimes(V.*reshape(expLam, 1, n_elm(1) + n_elm(2), []), pagemldivide(V, eye(n_elm(1) + n_elm(2))));


if n > n_elm(1) + n_elm(2)
    % Lidar bias present 
    F_new = repmat(eye(n), 1, 1, n_bank);
    F_new(1:(n_elm(1) + n_elm(2)), 1:(n_elm(1) + n_elm(2)), :) = F;
    F = F_new;
end
end



