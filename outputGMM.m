function [x, P] = outputGMM(x_nodes, P_nodes, w_nodes, n_nodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Combine Gaussian nodes in to single GMM estimate of mean and
% covariance 
% 
% Inputs: 
%   x_nodes     Array of Gaussian means for each node (n, n_nodes)
%   P_nodes     Array of Gaussian covariance matrices for each node (n, n,
%               n_nodes)
%   w_nodes     Array of weights for each Gaussian node (n_nodes, 1)
%   n_nodes     Number of nodes in GMM 
% 
% Output:
%   x           GMM mean estimate
%   P           GMM covariance estiamte 
% 
% Author: Rachel Mamich 
% Date Created: 04/01/25
% Date Modified: 04/02/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = x_nodes*w_nodes;
e_minus = x_nodes - x;
e_minus = reshape(e_minus, [], 1, n_nodes);
P = sum(pagemtimes(reshape(w_nodes, 1, 1, []), P_nodes + pagemtimes(e_minus, 'none', e_minus, 'transpose')), 3);
end