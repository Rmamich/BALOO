function [kron_hBB] = kronh_page(B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Calculates the Kronecker product taking two sided symmetry into
% account. 
%           vech(BXB') = (B kron_h B)vech(X)
% Assumes: X is symmetric (m x m)
% 
% Inputs: 
%   B           matrix with size (n x m) 
%   B           matrix with size (n x m)
% 
% Outputs:
%   kron_uBA    Kronecker product levereging single sided symmetry
% 
% Author: Rachel Mamich 
% Date Created: 03/27/23
% Date Modified: 04/07/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [n, m, n_nodes] = size(B);
    kronBB = pagekron(B, B);

    % kronBB = kronBB(:, :, 1);
    % n_nodes = 1;
    
    indexVec = 1:m^2;
    indexMat = reshape(indexVec, m, m)';
    mmp1div2 = m*(m+1)/2;
    n2 = n^2;
    kron_hBB = zeros(n2, mmp1div2, n_nodes)*B(1, 1);
    count = 1;
    % Combine columns for symmetry 
    for i = 1:m
        kron_hBB(:, count,  :) = kronBB(:, indexMat(i, i), :);
        count = count + 1;
        for j = i+1:m
            kron_hBB(:, count, :) = kronBB(:, indexMat(i, j), 1) + kronBB(:, indexMat(j, i), :);
            count = count +1;
        end
    end
    
    % Reduce Duplicate Rows
    indexVec = 1:n^2;
    indexMat = reshape(indexVec, n, n);
    rowKeep = vech(indexMat);
    kron_hBB = kron_hBB(rowKeep, :, :);
end