function [transitionProbMat] = transitionalProbabilityMatIMM(n_nodes, staysTheSame )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Creates a transitional probability matrix for use in an
% interacting multiple model filter (IMM).
% 
% Inputs:
%   n_nodes             number of nodes in the system
%   staysTheSame        probability that the mode does not change (diagonal
%                           values in the matrix. If no input given, will 
%                           default to 2/3. 
% 
% Outputs:          
%   transitionProbMat   transitional probability matrix 
%
% Note to self: There are probably more creative options here, but this was
% the most obvious approach. If things are super weird, try something more
% intellegent... maybe favor earlier nodes in the off diagonal elements? 
% 
% Created: 04/03/2025
% Modified: 04/03/2025
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n_nodes >1
    if nargin <2
        staysTheSame = 2/3;
    end
    otherValues = (1-staysTheSame)/(n_nodes-1);
    transitionProbMat = ones(n_nodes)*otherValues;

    idxMat = reshape(1:(n_nodes^2), n_nodes, n_nodes);
    idxDiag = diag(idxMat);

    transitionProbMat(idxDiag) = transitionProbMat(idxDiag)/otherValues*staysTheSame;
else
    transitionProbMat = 1;
end

end