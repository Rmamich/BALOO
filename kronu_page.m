function [kron_uBA] = kronu_page(B, A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Calculates the Kronecker product taking one sided symmetry into
% account. 
%           (AXB')(:) = (B kron_u A)vech(X)
% Assumes: X is symmetric (m x m)
% 
% Inputs: 
%   B           matrix with size (n x m) 
%   A           matrix with size (n x m)
% 
% Outputs:
%   kron_uBA    Kronecker product levereging single sided symmetry
% 
% Author: Rachel Mamich 
% Date Created: 03/27/23
% Date Modified: 04/07/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [n, m, n_nodes] = size(B);
    kronBA = pagekron(B, A);


    
    indexVec = 1:m^2;
    indexMat = reshape(indexVec, m, m)';
    
    kron_uBA = zeros(n^2, m*(m+1)/2, n_nodes);

    count = 1;
    
    % Combine columns for symmetry
    for i = 1:m
        kron_uBA(:, count, :) = kronBA(:, indexMat(i, i), :);
        count = count + 1;
        for j = i+1:m
            kron_uBA(:, count, :) = kronBA(:, indexMat(i, j), :) + kronBA(:, indexMat(j, i), :);
            count = count + 1;
        end
    end
end