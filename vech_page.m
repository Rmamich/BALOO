function vech = vech_page(matA)
% Function return elements from and below the main diagonal, then
% stacking by column to have a vector of size K*(K+1)/2.
% For example:  matA = [b11 b12 b13;b21 b22 b23;b31 b32 b33] is (3x3)
%               vech(matA) returns: [b11 b12 b13 b22 b23 b33]' is (6x1)
% Author: Binh T. Pham. Date: 2018-02-15.
[M,N, n_nodes] = size(matA);

if (M == N)
    vech  = zeros(M*(M+1)/2, n_nodes);
    for ii=1:n_nodes
        UT = reshape(triu(matA(:, :, ii)')', [] , 1);
        vech(:, ii) = reshape(UT(UT~=0), [], 1);
    end
else
     error('Input must be a symmetric matrix.')
end
end