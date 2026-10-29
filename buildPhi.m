function [F, A] = buildPhi(x, gravParam, dt, beta, n_elm, relNavFlag, CRTBPFlag, t)
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

xt = x(1:n_elm(1));
xc = x((1+n_elm(1)):n_elm(1)+n_elm(2));

if nargin > 6
    if ~CRTBPFlag
        % Normal TBP + whatever
        % Chaser
        Ac = zeros(n_elm(2));
        Ac(1:3, 4:6) = eye(3);
        [dadr] = gravJac(xc(1:3), gravParam);
        Ac(4:6, 1:3) = dadr;
        if n_elm(2) == 9
            %     Ac(4:6, 7:9) = eye(3);
            Ac(7:9, 7:9) = -eye(3)*beta(2); 
        end

        % Target
        At = zeros(n_elm(1));
        At(1:3, 4:6) = eye(3);
        [dadr] = gravJac(xt(1:3), gravParam);
        At(4:6, 1:3) = dadr;
        if n_elm(1) == 9
            At(4:6, 7:9) = eye(3);
            At(7:9, 7:9) = -eye(3)*beta(1); % Assumes constant burn
        end
    else
        % CRTBP -- Assumes that the state passes in is in nondim,
        % rotational frame, also assumes that beta is nondim 
        % Chaser
        Ac = zeros(n_elm(2));
        Ac(1:3, 4:6) = eye(3);
        [dadr] = gravJac_CRTBP(xc(1:3), gravParam);
        Ac(4:6, 1:3) = dadr;
        Ac(4:6, 4:6) = [0, 2, 0; -2, 0, 0; 0, 0, 0];
        if n_elm(2) == 9
            %     Ac(4:6, 7:9) = eye(3);
            Ac(7:9, 7:9) = -eye(3)*beta(2); % Assumes that beta is nondim
        end

        % dX_nonDimR_dX_dimN_c = dX_nonDimRot_dX_dimN(xc, t, gravParam);
        % % Ac is now wrt dimensional and inertial frame
        % Ac = dX_nonDimR_dX_dimN_c*Ac*dX_nonDimR_dX_dimN_c';
        % Target

        At = zeros(n_elm(1));
        At(1:3, 4:6) = eye(3);
        [dadr] = gravJac_CRTBP(xt(1:3), gravParam);
        At(4:6, 1:3) = dadr;
        At(4:6, 4:6) = [0, 2, 0; -2, 0, 0; 0, 0, 0];
        if n_elm(1) == 9
            At(4:6, 7:9) = eye(3);
            At(7:9, 7:9) = -eye(3)*beta(1); % Assumes that beta is nondim
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
    Ac = zeros(n_elm(2));
    Ac(1:3, 4:6) = eye(3);
    [dadr] = gravJac(xc(1:3), gravParam);
    Ac(4:6, 1:3) = dadr;
    if n_elm(2) == 9
        %     Ac(4:6, 7:9) = eye(3);
        Ac(7:9, 7:9) = -eye(3)*beta(2);
    end

    % Target
    At = zeros(n_elm(1));
    At(1:3, 4:6) = eye(3);
    [dadr] = gravJac(xt(1:3), gravParam);
    At(4:6, 1:3) = dadr;
    if n_elm(1) == 9
        At(4:6, 7:9) = eye(3);
        At(7:9, 7:9) = -eye(3)*beta(1); % Assumes constant burn
    end

end


if ~relNavFlag
    % Full inertial
    A = [At, zeros(n_elm(1), n_elm(2)); zeros(n_elm(2), n_elm(1)), Ac];
else
    % Rel Nav with different configs
    if n_elm(1) == 6 && n_elm(2) == 9
        A = [At, At-Ac(1:6, 1:6), zeros(n_elm(1), 3); zeros(n_elm(2), n_elm(1)), Ac];
    elseif n_elm(1) == 6 && n_elm(2) == 6
        A = [At, At-Ac;... 
            zeros(6), Ac];
    elseif n_elm(1) == 9 && n_elm(2) == 9
        A = [At, [At(1:6, 1:6)-Ac(1:6, 1:6), zeros(6, 3); zeros(3, 9)];...
            zeros(9), Ac];
    elseif n_elm(1) == 9 && n_elm(2) == 6
        A = [At, [At(1:6, 1:6)-Ac(1:6, 1:6), zeros(6, 3); zeros(3, 9)]; zeros(n_elm(2), n_elm(1)), Ac];
    end
end


% Convert to discrete time 
F = expm(A*dt);

if length(x) > n_elm(1) + n_elm(2)
    % Lidar bias present 
    F = blkdiag(F, 1);
end
end



