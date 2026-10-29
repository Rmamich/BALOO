function [A, B, Mw] = constructAB(F, H, bufferSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Constructs observability gramians, A, B, and stm's necessary to
% perform correlation method process nosie adaptation
% 
% Inputs:
%   F           time history of stms assumes size (nA, nA, buffer+1) where
%                   each 3D layer props the previous time step to the 
%                   "current" one
%   H           time history of measurement mapping matrices 
%                   (nM, nA, buffer+1)
%   bufferSize  size of buffer considered   
% 
% Outputs:
%   A           A coefficient matrices (nA, nM, buffer+1)
%   B           B coefficient matrices (nA, nA, buffer)
% 
% Author: Rachel Mamich 
% Date Created: 03/11/23
% Date Modified: 03/11/23
% Source Material: Adaptive Algorithms for Identification of Symmetric and
% Positive Definite Matrices Ch4 and Design of measurement difference
% autocovariance method for estimation of process anf measuremebnt noise
% covariances by J. Dunik O. Kost, and O. Straka 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhiSize = size(F(:, :, 1), 1);
m = size(H, 1);
Mj_jmbuffp1 = zeros(PhiSize);
Mjm1_jmbuff = zeros(PhiSize);
Mwjm1_jmbuff = zeros(bufferSize*m, PhiSize*(bufferSize-1));
Mwj_jmbuffp1 = zeros(bufferSize*m, PhiSize*(bufferSize-1));
Mw = zeros(bufferSize*m, PhiSize*bufferSize);

Phi_kk_jmbufferSizep1 = eye(size(Mj_jmbuffp1));
Phi_kk_jmbufferSize = eye(size(Mj_jmbuffp1));
Phi_bufferHist_jmbufferp1 = NaN([size(Mj_jmbuffp1), bufferSize]);

% build Mj_jmbuffp1
for kk = 0:(bufferSize-1)
    % Includes time step j, excludes j-buffer
    Mj_jmbuffp1 = Mj_jmbuffp1 + Phi_kk_jmbufferSizep1'*H(:, :, bufferSize-kk)'*H(:, :, bufferSize-kk)*Phi_kk_jmbufferSizep1;
    if kk < (bufferSize-1)
        Phi_kk_jmbufferSizep1  = F(:, :, bufferSize-kk - 1) * Phi_kk_jmbufferSizep1;
    end
end

% build Mjm1_jmbuff
for kk = 0:(bufferSize-1)
    % Includes time step j-buffer, excludes step j
    Mjm1_jmbuff = Mjm1_jmbuff + Phi_kk_jmbufferSize'*H(:, :, bufferSize-kk+1)'*H(:, :, bufferSize-kk+1)*Phi_kk_jmbufferSize;
    Phi_kk_jmbufferSize  = F(:, :, bufferSize-kk) * Phi_kk_jmbufferSize;
end

Phi_bufferHist_jmbufferp1(:, :, bufferSize) = eye(PhiSize); % map from initial time in buffer to the following times
for kk = 1:(bufferSize-1)
    Phi_bufferHist_jmbufferp1(:, :, bufferSize-kk) = F(:, :, bufferSize-kk) * Phi_bufferHist_jmbufferp1(:, :, bufferSize-kk+1);
end


% Build A coeffs
A = NaN(PhiSize, m, bufferSize+1);
Mj_jmbuffp1INV = Mj_jmbuffp1\eye(PhiSize);
Mjm1_jmbuffINV = Mjm1_jmbuff\eye(PhiSize);

A(:, :, 1) = Mj_jmbuffp1INV*Phi_bufferHist_jmbufferp1(:, :, 1)'*H(:, :, 1)';  % A0
for kk = 2:(bufferSize)
    % Ai
    A(:, :, kk) = Mj_jmbuffp1INV*Phi_bufferHist_jmbufferp1(:, :, kk)'*H(:, :, kk)' - F(:, :, bufferSize)*Mjm1_jmbuffINV*F(:, :, bufferSize)'*Phi_bufferHist_jmbufferp1(:, :, kk)'* H(:, :, kk)';
end
% Am 
A(:, :, bufferSize+1) = -F(:, :, bufferSize)*Mjm1_jmbuffINV*H(:, :, bufferSize+1)';


% build Mw -- k-3 -- k -- previous back -- green
for kk = 1:(bufferSize-1)
    Mwjm1_jmbuff(((kk-1)*m+1):kk*m, ((kk-1)*PhiSize+1):kk*PhiSize) = H(:, :, kk+1);
    for jj = (kk+1):(bufferSize-1)
        jjminus1 = jj-1;
        Mwjm1_jmbuff(((kk-1)*m+1):kk*m, ((jj-1)*PhiSize+1):jj*PhiSize) = Mwjm1_jmbuff(((kk-1)*m+1):kk*m, ((jjminus1-1)*PhiSize+1):jjminus1*PhiSize)*F(:, :, jj);
    end
end

populateMat = zeros(4,3 );
% build Mw -- k-2 -- k+1 -- current back -- orange
for kk = 1:(bufferSize-1)
    Mwj_jmbuffp1(((kk-1)*m+1):kk*m, ((kk-1)*PhiSize+1):kk*PhiSize) = H(:, :, kk);
    populateMat(kk, kk) = 1;
    for jj = (kk+1):(bufferSize-1)
        jjminus1 = jj-1;
        Mwj_jmbuffp1(((kk-1)*m+1):kk*m, ((jj-1)*PhiSize+1):jj*PhiSize) = Mwj_jmbuffp1(((kk-1)*m+1):kk*m, ((jjminus1-1)*PhiSize+1):jjminus1*PhiSize)*F(:, :, jjminus1);
        populateMat(kk, jj) = 1;
    end
end
Mw(1:bufferSize*m, 1:PhiSize*(bufferSize-1)) = Mwj_jmbuffp1;
Mw(1:m, PhiSize*(bufferSize-1)+1:PhiSize*bufferSize) = Mwj_jmbuffp1(1:m, PhiSize*(bufferSize-2)+1:(bufferSize-1)*PhiSize)*F(:, :, bufferSize-1);
Mw(m+1:bufferSize*m, PhiSize+1:PhiSize*bufferSize) = Mwjm1_jmbuff(1:(bufferSize-1)*m, :);

% build O^T -- k-3 -- k -- previous back -- green
Omatjm1_jmbuff = zeros(m*bufferSize, PhiSize);
for kk = 1:(bufferSize-1)
    Omatjm1_jmbuff((kk-1)*m+1:kk*m, :) = H(:, :, kk+1)*Phi_bufferHist_jmbufferp1(:, :, kk+1)*F(:, :, bufferSize);
end
Omatjm1_jmbuff((bufferSize-1)*m+1:bufferSize*m, :) = H(:, :, bufferSize+1);

% build O^T -- k-2 -- k+1 -- current back -- orange
Omatj_jmbuffp1 = zeros(m*bufferSize, PhiSize);
for kk = 1:(bufferSize-1)
    Omatj_jmbuffp1((kk-1)*m+1:kk*m, :) = H(:, :, kk)*Phi_bufferHist_jmbufferp1(:, :, kk);
end
Omatj_jmbuffp1((bufferSize-1)*m+1:bufferSize*m, :) = H(:, :, bufferSize);

OTMwjm1_jmbuff = Omatjm1_jmbuff'*Mwjm1_jmbuff; % pervious
OTMwj_jmbuffp1 = Omatj_jmbuffp1'*Mwj_jmbuffp1; % current 

A_firstBracket = [Mj_jmbuffp1INV*Omatj_jmbuffp1', zeros(PhiSize, m)]; % current  
A_secondBracket = [zeros(PhiSize, m),  F(:, :, bufferSize)*Mjm1_jmbuffINV*Omatjm1_jmbuff']; % previous 

A_test =  A_firstBracket-A_secondBracket;

A_3D = reshape(A_test, PhiSize, m, []);

B_firstBracket = [Mj_jmbuffp1INV*OTMwj_jmbuffp1, eye(PhiSize)]; % current  
B_secondBracket = [zeros(PhiSize),  F(:, :, bufferSize)*Mjm1_jmbuffINV*OTMwjm1_jmbuff]; % previous 

B_test = B_firstBracket-B_secondBracket;

B = reshape(B_test, PhiSize, PhiSize, []);



% B = NaN(PhiSize, PhiSize);
% 
% B(:, :, 1) = Mj_jmbuffp1INV*Phi_bufferHist_jmbufferp1(:, :, 1)'*(H(:, :, 1)'*H(:, :, 1));
% 
% for ii = 2:(bufferSize-1)
%     Mj_jmiip1 = zeros(PhiSize);
%     Phikk_jmiip1 = eye(PhiSize);
%     for kk = 0:(ii-1)
%         Mj_jmiip1 = Mj_jmiip1 + Phikk_jmiip1'*H(:, :, ii - kk)'*H(:, :, ii - kk)*Phikk_jmiip1;
%         if kk < (ii-1)
%             Phikk_jmiip1  = F(:, :, ii-kk-1) * Phikk_jmiip1;
%         end
%     end
%     Mjm1_jmiip1 = zeros(PhiSize);
%     Phikk_jmiip1 = eye(PhiSize);
%     for kk = 0:(ii-2)
%         Mjm1_jmiip1 = Mjm1_jmiip1 + Phikk_jmiip1'*H(:, :, ii - kk)'*H(:, :, ii - kk)*Phikk_jmiip1;
%         Phikk_jmiip1  = F(:, :, ii-kk-1) * Phikk_jmiip1;
%     end
%     Phijmiip1_jmbuff = Phi_bufferHist_jmbufferp1(:, :, ii)*F(:, :, bufferSize);
% 
%     B(:, :, ii) = Mj_jmbuffp1INV*Phi_bufferHist_jmbufferp1(:, :, ii)'*Mj_jmiip1 - F(:, :, bufferSize)*Mjm1_jmbuffINV*Phijmiip1_jmbuff'*Mjm1_jmiip1;
% end
% 
% Mjm1_jmbuffp1 = zeros(PhiSize);
% Phikk_jmbuffp1 = eye(PhiSize);
% for kk = 0:(bufferSize-2)
%     Mjm1_jmbuffp1 = Mjm1_jmbuffp1 + Phikk_jmbuffp1'*H(:, :, bufferSize - kk)'*H(:, :, bufferSize - kk)*Phikk_jmbuffp1;
%     Phikk_jmbuffp1  = F(:, :, bufferSize-kk-1) * Phikk_jmbuffp1;
% end
% 
% B(:, :, bufferSize) = eye(PhiSize) - F(:, :, bufferSize)*Mjm1_jmbuffINV*F(:, :, bufferSize)'*Mjm1_jmbuffp1;


end

