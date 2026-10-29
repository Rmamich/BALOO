function [A, B, Mw] = constructAB_page(F, H, bufferSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Constructs observability gramians, A, B, and stm's necessary to
% perform correlation method process nosie adaptation
% 
% Inputs:
%   F           time history of stms assumes size (nA, nA, buffer+1, n_nodes) where
%                   each 3D layer props the previous time step to the 
%                   "current" one
%   H           time history of measurement mapping matrices 
%                   (nM*(buffer+1), nA, n_nodes)
%   bufferSize  size of buffer considered   
% 
% Outputs:
%   A           A coefficient matrices (nA, nM*buffer+1, n_nodes)
%   B           B coefficient matrices (nA, nAbuffer+1, n_nodes)
% 
% Author: Rachel Mamich 
% Date Created: 03/11/23
% Date Modified: 04/04/25
% Source Material: Adaptive Algorithms for Identification of Symmetric and
% Positive Definite Matrices Ch4 and Design of measurement difference
% autocovariance method for estimation of process anf measuremebnt noise
% covariances by J. Dunik O. Kost, and O. Straka 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m, n, n_nodes] = size(H);
m = m/(bufferSize+1);


eye_nxn_nodes = repmat(eye(n), 1, 1, n_nodes); % useful in several spots

% Preallocte memory for Phi stack 
Phi_bufferHist_jmbufferp1 = NaN([n*bufferSize, n, n_nodes]);
Phi_bufferHist_jmbuffer = NaN([n*bufferSize, n, n_nodes]); % map from initial time in buffer to the following times j - bufferSize


% Initialize last Phi in each stack to be eye(n);
Phi_bufferHist_jmbufferp1((1 + n*(bufferSize-1)):(n*bufferSize), 1:n, :) = eye_nxn_nodes; % map from initial time in buffer to the following times j - bufferSize + 1
Phi_bufferHist_jmbuffer((1 + n*(bufferSize-1)):(n*bufferSize), 1:n, :) = eye_nxn_nodes;
% Populate Phi stacks 
for kk = 2:(bufferSize)
    idx = bufferSize - kk + 1;
    idxPrev = bufferSize - kk + 2;
    Phi_bufferHist_jmbufferp1((1 + n*(idx-1)):(n*idx), :, :) = pagemtimes(F(:, :, :, idx), Phi_bufferHist_jmbufferp1((1 + n*(idxPrev-1)):(n*idxPrev), :, :));
    Phi_bufferHist_jmbuffer((1 + n*(idx-1)):(n*idx), :, :) = pagemtimes(F(:, :, :, idxPrev), Phi_bufferHist_jmbuffer((1 + n*(idxPrev-1)):(n*idxPrev), :, :));
end



% Preallocate memory for blk diag H, used to get Observability matrix 
H_blkDiag = zeros(m*(bufferSize+1), n*(bufferSize+1), n_nodes);
% Populate H blk diag 
for jj = 1:(bufferSize+1)
    H_blkDiag((1 + m*(jj-1)):(jj*m), (1+ n*(jj-1)):(jj*n), :) = H((1+ m*(jj-1)):(jj*m), :, :);
end

%Calculate Observability matrix and gramian 
Omatj_jmbuffp1 = pagemtimes(H_blkDiag(1:(bufferSize*m), 1:(bufferSize*n), :), Phi_bufferHist_jmbufferp1);
Mj_jmbuffp1 = pagemtimes(Omatj_jmbuffp1, 'transpose', Omatj_jmbuffp1, 'none');

Omatjm1_jmbuff = pagemtimes(H_blkDiag((1+m):end, (1 + n):end, :), Phi_bufferHist_jmbuffer);
Mjm1_jmbuff = pagemtimes(Omatjm1_jmbuff, 'transpose', Omatjm1_jmbuff, 'none');


% Invert observability gramian 
Mj_jmbuffp1INV = pagemldivide(Mj_jmbuffp1,  eye(n));
Mjm1_jmbuffINV = pagemldivide(Mjm1_jmbuff, eye(n));


% Get A matrix components for each buffer 
A1 = pagemtimes(Mj_jmbuffp1INV, 'none', Omatj_jmbuffp1, 'transpose');
A2 = pagemtimes(F(:, :, :, bufferSize), pagemtimes(Mjm1_jmbuffINV, 'none', Omatjm1_jmbuff, 'transpose'));

% Combine 
A = [A1, zeros(n, m, n_nodes)] - [zeros(n, m, n_nodes), A2];

% Preallocate Mw matrix for each buffer and for all components 
Mwjm1_jmbuffp1 = zeros(m*(bufferSize), n*(bufferSize-1), n_nodes);
Mwjm2_jmbuff = zeros(m*(bufferSize), n*(bufferSize-1), n_nodes);
Mw = zeros(m*(bufferSize+1), n*(bufferSize), n_nodes);

for kk = 1:(bufferSize+2)
    Phi_temp = eye_nxn_nodes;
    for ii = kk:(bufferSize)
        HPhi = pagemtimes(H((1+m*(kk-1)):(m*kk), :, :), Phi_temp);
        if kk <= bufferSize && ii <= (bufferSize-1)
            Mwjm1_jmbuffp1(1+(m*(kk-1)):(m*kk), 1+(n*(ii-1)):(n*ii), :) = HPhi;
        end
        if kk >= 2 
            kkminus1 = kk-1;
            iiminus1 = ii - 1;
            Mwjm2_jmbuff(1+(m*(kkminus1-1)):(m*kkminus1), 1+(n*(iiminus1-1)):(n*iiminus1), :) = HPhi;
        end
        Mw(1+(m*(kk-1)):(m*kk), 1+(n*(ii-1)):(n*ii), :) = HPhi;
        Phi_temp = pagemtimes(Phi_temp,F(:, :, :, ii));
    end
end

% Get A matrix components for each buffer 
B1 = pagemtimes(A1,Mwjm1_jmbuffp1);
B2 = pagemtimes(A2, Mwjm2_jmbuff);

% Calculate B
B = [B1, eye_nxn_nodes] - [zeros(n, n, n_nodes), B2];


end

