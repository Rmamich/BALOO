function [H, ABlock] = formMeasMatrix_variableBufferSize_page(Ahist, Bhist, bufferSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Creates the "measurement mapping" matrix from the vectorized 
% process noise matrix to the autocovariance of the measurement residuals.
%
% Inputs: 
%   Ahist   history of the A matrices over the buffer window 
%           (nA*(buffer+1), nM*(buffer+1), n_nodes)
%   Bhist   history of the B matrices over the buffer window
%           (nA*buffer, nA*buffer, n_nodes)
% 
% Outputs: 
%   H       "Measurement mapping matrix" for least square est of process
%               noise psd 
%   ABlock  stacked A matrix that can be used to remove determinsitic
%               knowledge of measurement noise from autocovariance
%               bk - ABlock*vech(R) = BBlock*vech(Q)
%               (n*(n+1)/2 + n^2*buffSize x nM*(nM+1)/2)
% Author: Rachel Mamich 
% Date Created: 03/17/23
% Date Modified: 04/07/25
% Modified from code produced by Rahul Moghe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[nA, nM, n_nodes] = size(Ahist);
nA = nA/(bufferSize+1);
nM = nM/(bufferSize+1);

ut1_nA = nA*(nA+1)/2;
ut_nM = nM*(nM+1)/2;

nA2 = nA^2;

ABlock = zeros((bufferSize)*nA2+ ut1_nA, ut_nM, n_nodes);
BBlock = zeros((bufferSize)*nA2+ ut1_nA, ut1_nA, n_nodes);
 
% A is relevant for estimating measuremnt noise matrix 
% First Entry A
Atemp = zeros(ut1_nA, ut_nM, n_nodes);
for i = 1:(bufferSize+1)
    Atemp = Atemp + kronh_page(Ahist(1:nA, (1+nM*(i-1)):(nM*(i)),:)); % All A at k
end

ABlock(1:ut1_nA, :, :) = Atemp;

%   Ahist   history of the A matrices over the buffer window 
%           (nA, nM, buffer+1, buffer+1)
%   Bhist   history of the B matrices over the buffer window
%           (nA, nA, buffer, buffer)

% Loop through to populate the rest of A 
for kk = 2:(bufferSize+1)
    Atemp = zeros(nA2, ut_nM, n_nodes);
    AkmpmeasIdx = 1;
    for i = kk:(bufferSize +1)
        % Atemp = Atemp + kronu(Ahist(:,:,i+(kk-1), 1), Ahist(:,:,i,kk));
        Atemp = Atemp + kronu_page(Ahist((1+nA*(kk-1)):(nA*(kk)), (1+nM*(AkmpmeasIdx-1)):(nM*(AkmpmeasIdx)), :), Ahist(1:nA, (1+nM*(i-1)):(nM*(i)), :));
        AkmpmeasIdx = AkmpmeasIdx + 1;
    end
    ABlock((nA2*(kk-2)+1+ ut1_nA):(nA2*(kk-1) + ut1_nA), :, :) = Atemp;
end

% First Entry B 
Btemp = zeros(ut1_nA, ut1_nA, n_nodes);
for i = 1:bufferSize
    Btemp = Btemp + kronh_page(Bhist(1:nA, (1+nA*(i-1)):(nA*(i)),:));
end

BBlock(1:ut1_nA, :, :) = Btemp;

% Loop through to populate the rest of B
for kk = 2:bufferSize
    Btemp = zeros(nA2, ut1_nA, n_nodes);
    BkmpmeasIdx = 1;
    for i = kk:bufferSize
        % Btemp = Btemp + kronu(Bhist(:,:,i+kk-1, 1), Bhist(:,:,i,kk));
        Btemp = Btemp + kronu_page(Bhist((1+nA*(kk-1)):(nA*(kk)), (1+nA*(BkmpmeasIdx-1)):(nA*(BkmpmeasIdx)),:), Bhist(1:nA, (1+nA*(i-1)):(nA*(i)),:, :));
        BkmpmeasIdx = BkmpmeasIdx + 1;
    end
    BBlock((nA2*(kk-2)+1+ ut1_nA):(nA2*(kk-1) + ut1_nA), :, :) = Btemp;
end

 
H =  BBlock;


end
