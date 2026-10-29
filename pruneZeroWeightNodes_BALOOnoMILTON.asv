function [x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, n_nodes, nodeIDs, splitWait] = pruneZeroWeightNodes_BALOOnoMILTON(epsPrune, x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, nodeIDs, splitWait, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Pruning function for BALOO. Nodes are removed when their weights
% fall below the threshold epsPrune.
%
% Inputs:
%   
%
% Outputs:
%
%
% Author: Rachel Mamich
% Date Created: 04/15/25
% Date Modified: 04/15/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk = find(~isnan(x_nodes(1, 1, :)), 1, 'last');
[~, n_nodes, ~] = size(x_nodes);

zeroNodes = w_nodes(:, kk) < epsPrune;

% zeroNodes_2 = w_nodes(:, kk) < 1e-10;
% 
% if any(zeroNodes_2)
%     w_nodes(:, kk)
% end

nTrunk = 10; 
n_nodes = length(w_nodes(:, kk));

idList = [w_nodes(:, kk), [1:n_nodes]'];
idListSort = sortrows(idList, 'descend');

if n_nodes > nTrunk
    zeroNodes = idListSort(nTrunk+1:end, 2);
    zeroNodes = sort(zeroNodes);
end




if any(zeroNodes) 
    nodeList = 1:(length(zeroNodes));
    fprintf('Killing off at time step %d\n', kk)
    fprintf('Killing off %d\n', zeroNodes)
    x_nodes(:, zeroNodes, :) = [];
    P_nodes(:, :, zeroNodes, :) = [];
    w_nodes(zeroNodes, :) = [];
    w_nodes(:, kk) = w_nodes(:, kk) / sum(w_nodes(:, kk));
    Pminus_nodes(:, :, zeroNodes, :) = [];
    Phi(:, :, zeroNodes, :) = [];
    H_hist(:, :, zeroNodes, :) = [];
    riPlus_nodes(:, zeroNodes, :) = []; 
    riMinus_nodes(:, zeroNodes, :) = [];
    Pri_nodes(:, :, zeroNodes, :) = [];
    PriMinus_nodes(:, :, zeroNodes, :) = [];
    rWindow_nodes(:, zeroNodes) = [];
    Fwindow_nodes(:, :, zeroNodes, :) = [];
    Hwindow_nodes(:, :, zeroNodes) = [];
    
    splitWait(zeroNodes) = [];
    % fullEnoughSplit(zeroNodes) = [];
    % fullEnoughAdapt(zeroNodes) = [];
    
    
    
    [~, n_nodes, ~] = size(x_nodes);
    nodeIDs(zeroNodes) = [];
   
    
    if isempty(n_nodes)
        fprintf('Hey Silly Goose, you killed off all of your nodes\n')
    end
end
end