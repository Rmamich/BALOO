function [out] = nonDimRotational_dimInertial_conversion(in, t, L, n, tCharacteristic, inputFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Convert between nondimensional vectors in the rotational frame
% and dimensional vectors in the inertial frame
% 
% Inputs: 
%   in      input vector 
%   t       time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nStates = length(in);
if strcmp(inputFlag, 'R2N') 
    % Converts from nondimensional, rotational frame used in CRTBP to 
    % dimensional and inertial  
    r_R = in(1:3)*L; % km
    v_R = in(4:6)*L*n; % km/s
    if nStates == 9
        a_R = in(7:9)*(L)/tCharacteristic/tCharacteristic;  %km/s^2
    end

    omega = [0;0;n]; % 1/s

    % rotational to inertial DCM
    NR = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1];
    r_N = NR * r_R; % km

    % accounting for angular velocity of the frame
    % assumes that rotation is omega = [0; 0; n]
    % so NR*v_r + cross(omega, r_N);
    v_N = NR*v_R + cross(omega, r_N); % km/s
    if nStates == 9
        a_N = NR*(a_R);
        % km/s^2 + 1/s* km/s                 + 1/s * 1/s * km
        out = [r_N; v_N; a_N];
    else
        out = [r_N; v_N];
    end

    
elseif strcmp(inputFlag, 'N2R')
    % Convert from the dimensional, inertial frame to the nondimensional
    % rotational frame for CRTBP 
    r_N = in(1:3); % km
    v_N = in(4:6); % km/s
    t = t/tCharacteristic; % []
    if nStates == 9
        a_N = in(7:9); % km/s^2
    end

    omega = [0;0;n]; % 1/s

    % rotational to inertial DCM
    RN = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1]';

    r_R = RN*r_N;

    v_R = RN*(v_N - cross(omega, r_N)); % km/s

    if nStates == 9
        a_R = RN*a_N;
        r_R = r_R/L;
        v_R = v_R/L/n;
        a_R = a_R/(L)*tCharacteristic*tCharacteristic;
        out = [r_R; v_R; a_R];
    else
        r_R = r_R/L;
        v_R = v_R/L/n;
        out = [r_R; v_R];
    end
elseif strcmp(inputFlag, 'accelN2R_posVelR')
    % Mostly for prop script, this case is for when you have a
    % nondimensional position and velocity, but the acceleration is in
    % nondimensional and in the inertial frame
    a_N = in(7:9); % km/s^2
    % rotational to inertial DCM
    RN = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1]';
    % a_R = -2*cross(omega, v_R) - cross(omega, cross(omega, r_R)) + 
    a_R = RN*a_N;
    out = a_R;
elseif strcmp(inputFlag, 'fullState_dimN2nondimR')
    nState = length(in);
    r_t_N = in(1:3); % m
    v_t_N = in(4:6); % m/s
    a_t_N = in(7:9); % m/s/s 
    t = t/tCharacteristic; % [] % Assume I am getting dimensaional time
    r_c_N = in(10:12); % m
    v_c_N = in(13:15); % m/s
    b_c_B = in(16:18); % m/s/s
    if nState == 19
        l_B = in(19); % m
    end

    omega = [0;0;n]; % 1/s

    % rotational to inertial DCM
    RN = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1]';

    r_t_R = RN*r_t_N;  % km

    v_t_R = RN*(v_t_N - cross(omega, r_t_N)); % km/s
    % v_t_R = RN*(v_t_N + [n*r_t_N(2); -n*r_t_N(1); 0]);
    
    v_t_NinR = RN*v_t_N;
    % -2*cross(omega, v_R) = [2*omega*v_t_NinR(2); -2*omega*v_t_NinR(1); 0];
    % - cross(omega, cross(omega, r_t_R)) = [omega^2*r_t_R(1); omega^2*r_t_R(2); 0];
    % if norm(a_t_N) >= 0.0001
    %     a_t_R = [2*n*v_t_NinR(2); -2*n*v_t_NinR(1); 0] + [n^2*r_t_R(1); n^2*r_t_R(2); 0] + RN*a_t_N; % km/s/s
    % else
    %     a_t_R = zeros(3, 1);
    % end
    a_t_R = RN*a_t_N;
    r_t_R = r_t_R/L; % [-]
    v_t_R = v_t_R/L/n; % [-]
    a_t_R = a_t_R/(L)*tCharacteristic*tCharacteristic; % [-]

    r_c_R = RN*r_c_N; % km
    v_c_R = RN*(v_c_N - cross(omega, r_c_N)); % km/s

    
    r_c_R = r_c_R/L; % [-]
    v_c_R = v_c_R/L/n; % [-]
    b_c_B = b_c_B/(L)*tCharacteristic*tCharacteristic; % [-]

    

    if nState == 19
        l_B = l_B/L; % [-]
        out = [r_t_R; v_t_R; a_t_R; r_c_R; v_c_R; b_c_B; l_B]; % [-]
    else
        out = [r_t_R; v_t_R; a_t_R; r_c_R; v_c_R; b_c_B]; % [-]
    end
elseif strcmp(inputFlag, 'accelerometer')
    out = in/(1000*L)*tCharacteristic*tCharacteristic; % dimensional to non dimensional 
elseif strcmp(inputFlag, 'fullState_nondimR2dimN')
    % Converts from nondimensional, rotational frame used in CRTBP to 
    % dimensional and inertial  
    nState = length(in);
    r_t_R = in(1:3)*L; % m
    v_t_R = in(4:6)*L*n; % m/s
    a_t_R = in(7:9)*(L)/tCharacteristic/tCharacteristic; % m/s/s 
    t = t/tCharacteristic; % [] % Assume I am getting dimensaional time
    r_c_R = in(10:12)*L; % m
    v_c_R = in(13:15)*L*n; % m/s
    b_c_B = in(16:18)*(L)/tCharacteristic/tCharacteristic; % m/s/s
    if nState == 19
        l_B = in(19)*L; % m
    end

    omega = [0;0;n]; % 1/s

    % rotational to inertial DCM
    NR = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1];
    r_t_N = NR * r_t_R; % m

    % accounting for angular velocity of the frame
    % assumes that rotation is omega = [0; 0; n]
    % so NR*v_r + cross(omega, r_N);
    v_t_N = NR*v_t_R + cross(omega, r_t_N); % km/s

    v_t_NinR = NR'*v_t_N;
    % 2*cross(omega, v_R) = -[2*omega*v_t_NinR(2); -2*omega*v_t_NinR(1); 0];
    % cross(omega, cross(omega, r_t_R)) = -[omega^2*r_t_R(1); omega^2*r_t_R(2); 0];
    % if norm(a_t_R) >= 0.0001
    %     a_t_N = NR*(a_t_R - [2*n*v_t_NinR(2); -2*n*v_t_NinR(1); 0] - [n^2*r_t_R(1); n^2*r_t_R(2); 0]);
    % else 
    %     a_t_N = zeros(3, 1);
    % end
    a_t_N = NR*a_t_R;
    % km/s^2 + 1/s* km/s                 + 1/s * 1/s * km

    r_c_N = NR * r_c_R; % m

    % accounting for angular velocity of the frame
    % assumes that rotation is omega = [0; 0; n]
    % so NR*v_r + cross(omega, r_N);
    v_c_N = NR*v_c_R + cross(omega, r_c_N); % km/s
    
    
    if nState == 19
        out = [r_t_N; v_t_N; a_t_N; r_c_N; v_c_N; b_c_B; l_B]; % [-]
    else
        out = [r_t_N; v_t_N; a_t_N; r_c_N; v_c_N; b_c_B]; % [-]
    end

elseif strcmp(inputFlag, 'DCM_R2N')
    t = t/tCharacteristic; % [] % Assume I am getting dimensaional time
    NR = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1];

    out = NR;


elseif strcmp(inputFlag, 'fullStateDot_nondimR2dimN')
    % Converts from nondimensional, rotational frame used in CRTBP to 
    % dimensional and inertial  
    nState = length(in);
    v_t_R = in(1:3)*L*n; % m/s
    a_t_R = in(4:6)*L*n*n; % m/s/s
    j_t_R = in(7:9)*(L)/tCharacteristic/tCharacteristic/tCharacteristic; % m/s/s/s 
    t = t/tCharacteristic; % [] % Assume I am getting dimensaional time
    v_c_R = in(10:12)*L*n; % m/s
    a_c_R = in(13:15)*L*n*n; % m/s/s
    bDot_c_B = in(16:18)*(L)/tCharacteristic/tCharacteristic/tCharacteristic; % m/s/s/s
    
    r_t_R = in(19:21)*L; % m  
    r_c_R = in(22:24)*L; % m 

    omega = [0;0;n]; % 1/s

    % rotational to inertial DCM
    NR = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1];
    r_t_N = NR * r_t_R; % m

    % accounting for angular velocity of the frame
    % assumes that rotation is omega = [0; 0; n]
    % so NR*v_r + cross(omega, r_N);
    v_t_N = NR*v_t_R + cross(omega, r_t_N); % km/s

    v_t_NinR = NR'*v_t_N;
    % 2*cross(omega, v_R) = -[2*omega*v_t_NinR(2); -2*omega*v_t_NinR(1); 0];
    % cross(omega, cross(omega, r_t_R)) = -[omega^2*r_t_R(1); omega^2*r_t_R(2); 0];
    % if norm(a_t_R) >= 0.0001
    a_t_N = NR*(a_t_R - [2*n*v_t_NinR(2); -2*n*v_t_NinR(1); 0] - [n^2*r_t_R(1); n^2*r_t_R(2); 0]);
    % else 
    %     a_t_N = zeros(3, 1);
    % end
    j_t_N = NR*j_t_R;
    % km/s^2 + 1/s* km/s                 + 1/s * 1/s * km

    r_c_N = NR * r_c_R; % m

    % accounting for angular velocity of the frame
    % assumes that rotation is omega = [0; 0; n]
    % so NR*v_r + cross(omega, r_N);
    v_c_N = NR*v_c_R + cross(omega, r_c_N); % km/s
    v_c_NinR = NR'*v_c_N;
    % 2*cross(omega, v_R) = -[2*omega*v_t_NinR(2); -2*omega*v_t_NinR(1); 0];
    % cross(omega, cross(omega, r_t_R)) = -[omega^2*r_t_R(1); omega^2*r_t_R(2); 0];
    % if norm(a_t_R) >= 0.0001
    a_c_N = NR*(a_c_R - [2*n*v_c_NinR(2); -2*n*v_c_NinR(1); 0] - [n^2*r_c_R(1); n^2*r_c_R(2); 0]);
    
    out = [v_t_N; a_t_N; j_t_N; v_c_N; a_c_N; bDot_c_B]; % [-]

% elseif trcmp(inputFlag, 'fullStateCovariance_nondimR2dimN_relNav_9x9')
%     % CRTBP -- Assumes that the state passes in is in nondim,
%     % rotational frame, also assumes that beta is nondim
%     % Chaser
%     dX_nonDimR_dX_dimN_c = dX_nonDimRot_dX_dimN(0, t, [mu, L, n, tCharacteristic]);
% 
%     % Target
%     dX_nonDimR_dX_dimN_t = dX_nonDimRot_dX_dimN(0, t, gravParam);
% 
% 
%     A = [At, [dX_nonDimR_dX_dimN_t(1:6, 1:6)-dX_nonDimR_dX_dimN_c(1:6, 1:6), zeros(6, 3); zeros(3, 9)];...
%         zeros(9), dX_nonDimR_dX_dimN_c];
        



else
    fprintf('Input flag does not match expected\n')
    out = in;
end

end

