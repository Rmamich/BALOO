function [x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi, H_hist, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, n_nodes, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait] = BALOOnoMILTON_splitCheck_NRHO(resConsider, epsilon_split, bufferSize, time, x_nodes, P_nodes, Pminus_nodes, w_nodes, Phi,  yWindow, H_hist, nodeIDs, x_nodes_all, P_nodes_all, Pminus_nodes_all, Phi_all, w_nodes_all, splitWait, RCF1_Body2_N, riPlus_nodes, riMinus_nodes, Pri_nodes, PriMinus_nodes, rWindow_nodes, Fwindow_nodes, Hwindow_nodes, fullEnough, eom, gravParam, beta_markov, accelerometer, gate, Rbody_N, RN_rot,relNavFlag, accelBias, n_elm, R)
% qkminus, Pkminus,  Awindow, Bwindow, Zwindow, vBarHist, uWindow_nodes,
% Xk_nodes, Qhat  =  qkminus, Pkminus, qk0, Pk0, Awindow, Bwindow, Zwindow, vBarHist, uWindow_nodes, Xk, Qhat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Splitting function for BALOO. Split is performed when the mean
% of the prefit measurment residual buffer exceedes the split limit,
% epsilon_split. The nodes will be split into 3, splitting only on the
% unmodeeld target acceleration.
%
% Inputs:
%   epsilon_split   The user defined limit for splitting nodes
%   n_nodes         The number of nodes in the system
%   rWindow_nodes   Time history of prefit measurement residuals for each
%                       node (nM, buffer+1, n_nodes)
%   Fwindow_nodes   Time history of stms assumes size (nA, nA, buffer+1, n_nodes)
%                   where each 3D layer props the previous time step to the
%                   "current" one
%   Hwindow_nodes   Time history of measurement mapping matrices
%                   (nM, nA, buffer+1, n_nodes)
%
% Outputs:
%
%
% Author: Rachel Mamich
% Date Created: 02/26/25
% Date Modified: 04/15/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nMmbp1, nA, n_nodes] = size(Hwindow_nodes);
splitWait = splitWait - 1;


if any(fullEnough)
    kk = find(~isnan(x_nodes(1, 1, :)), 1, 'last');
    buffCalcSplit = resConsider;
    if kk <= bufferSize-1
        timeStart = 1;
    else
        timeStart = kk-buffCalcSplit-1;
    end
    
    n_nodes_viable = fullEnough;

    
    num_nodes_viable = sum(fullEnough);

    nM = nMmbp1/(bufferSize+1);

    b = 2;
    s0 = (b-1)/(b^(resConsider) - 1);
    w_meas_favorNewer =      diag(ones(1, resConsider))/resConsider; %diag(flip(s0*b.^(0:(resConsider-1))));%

    % % nA = nAmbp1/(bufferSize);
    % PriMinus_nodes_buff = flip(PriMinus_nodes(:, :, n_nodes_viable, (timeStart+1):kk), 4);
    % PriMinus_nodes_buff = PriMinus_nodes_buff(:, :, :, 1:resConsider);
    % Pyyinv = zeros(nM*(resConsider), nM*(resConsider), num_nodes_viable);
    % for i = 1:(resConsider)
    %     Pyyinv((1+nM*(i-1)):(nM*i), (1+nM*(i-1)):(nM*i), :) = w_meas_favorNewer(i, i)*pagemldivide(PriMinus_nodes_buff(:, :, :, i), eye(nM));%w_meas_favorNewer(i, i)*;
    % end

    rWindow_nodes_splitViable = reshape(rWindow_nodes(:, n_nodes_viable), [], 1, num_nodes_viable);

   
    Rinv = repmat(R\eye(nM), resConsider, resConsider, num_nodes_viable);

    % SNIS = squeeze(pagemtimes(rWindow_nodes_splitViable((1:resConsider*nM), :, :), 'transpose', pagemtimes(Pyyinv, rWindow_nodes_splitViable(1:(resConsider*nM), :, :)),'none')/nM);
    SNIS = squeeze(pagemtimes(rWindow_nodes_splitViable((1:resConsider*nM), :, :), 'transpose', pagemtimes(Rinv, rWindow_nodes_splitViable(1:(resConsider*nM), :, :)),'none')/nM);
    mean_buff = sum(pagemtimes(reshape(rWindow_nodes_splitViable(1:resConsider*nM, :), nM, [], num_nodes_viable), w_meas_favorNewer), 2, 'omitnan'); %mean(rWindow_nodes_splitViable(:, 1:resConsider), 2, 'omitnan');
    norm_mean_buff = reshape(pagenorm(mean_buff), 1, []);
    % nodes_to_split = (norm_mean_buff) > epsilon_split;
    nodes_to_split = (SNIS) > epsilon_split;
    num_nodes_to_split = sum(nodes_to_split);

    % if  kk < 300
    %     % if n_nodes_viable > 1
    %     % % fprintf('norm mean buff %0.3e on %d \n', norm_mean_buff, repmat(kk, 1, n_nodes_viable))
    % if kk > 975 && kk < 1200
    %     fprintf('NIS %0.3e on %d \n', SNIS, repmat(kk, 1, n_nodes_viable))
    % 
    % end
    %     % else
    %     %     % fprintf('norm mean buff %0.3e on %d \n', norm_mean_buff, kk)
    %     %     fprintf('NIS %0.3e on %d \n', SNIS, kk)
    %     % end
    % 
    %     % if SNIS > 3.578130795349372/2*1.2
    %     %     fprintf('NIS > 1.5*bound on %d\n', kk)
    %     % end
    % end


    if num_nodes_to_split > 0
        % if  kk < 300

                fprintf('splitting on %d \n', kk)
                fprintf('SNIS %0.3e \n', SNIS)
      
        % end
        [~, nodes_all, ~] = size(x_nodes_all);
        [n, ~, N] = size(x_nodes);

        x_nodes_all_new = NaN(n, nodes_all, N);
        P_nodes_all_new = NaN(n, n, nodes_all, N);
        w_nodes_all_new = zeros(nodes_all, N);
        Pminus_nodes_all_new = NaN(n, n, nodes_all, N);
        Phi_all_new = NaN(n, n, nodes_all, N);


        
        new_nodes = num_nodes_to_split; % Create 1 new node
        n_nodes_new = n_nodes + new_nodes;

        rWindow_nodes_new = NaN(nMmbp1, n_nodes_new);
        Fwindow_nodes_new = zeros(nA, nA, n_nodes_new, bufferSize);
        Hwindow_nodes_new = zeros(nMmbp1, nA, n_nodes_new);
        
        
        
        x_nodes_new = NaN(n, n_nodes_new, N);
        P_nodes_new = NaN(n, n, n_nodes_new, N);
        Pminus_nodes_new = NaN(n, n, n_nodes, N);
        w_nodes_new = zeros(n_nodes_new, N);
        Phi_new = NaN(n, n, n_nodes_new, N);

        H_hist_new = NaN(nM, n, n_nodes_new, N);
        riPlus_nodes_new = NaN(nM, n_nodes_new, N);
        riMinus_nodes_new = NaN(nM, n_nodes_new, N);
        Pri_nodes_new  = NaN(nM, nM, n_nodes_new, N);
        PriMinus_nodes_new =NaN(nM, nM, n_nodes_new, N);

        

        
        splitWait_new = ones(n_nodes_new, 1)*resConsider*2;

        % Begin new buffers over windows with nothing-- keep all old
        % information from old buffers
        rWindow_nodes_new(:, 1:n_nodes) = rWindow_nodes;
        Fwindow_nodes_new(:, :, 1:n_nodes, :) = Fwindow_nodes;
        Hwindow_nodes_new(:, :, 1:n_nodes) = Hwindow_nodes; 
        
        
        x_nodes_new(:, 1:n_nodes, 1:kk) = x_nodes(:, :, 1:kk);
        P_nodes_new(:, :, 1:n_nodes, 1:kk) = P_nodes(:, :, :, 1:kk);
        Pminus_nodes_new(:, :, 1:n_nodes, 1:kk) = Pminus_nodes(:, :, :, 1:kk);
        w_nodes_new(1:n_nodes, 1:kk) = w_nodes(:, 1:kk);
        Phi_new(:, :, 1:n_nodes, 1:kk) = Phi(:, :, :, 1:kk);
        H_hist_new(:, :, 1:n_nodes, 1:kk) = H_hist(:, :, :, 1:kk);
        riPlus_nodes_new(:, 1:n_nodes, 1:kk) = riPlus_nodes(:, :, 1:kk);
        riMinus_nodes_new(:, 1:n_nodes, 1:kk)  = riMinus_nodes(:, :, 1:kk);
        Pri_nodes_new(:, :, 1:n_nodes, 1:kk)   = Pri_nodes(:, :, :, 1:kk);
        PriMinus_nodes_new(:, :, 1:n_nodes, 1:kk) = PriMinus_nodes(:, :, :, 1:kk);


        x_nodes_all_new(:, nodeIDs, 1:kk) = x_nodes(:, :, 1:kk);
        P_nodes_all_new(:, :, nodeIDs, 1:kk) = P_nodes(:, :, :, 1:kk);
        Pminus_nodes_all_new(:, :, nodeIDs, 1:kk) = Pminus_nodes(:, :, :, 1:kk);
        w_nodes_all_new(nodeIDs, 1:kk) = w_nodes(:, 1:kk);
        Phi_all_new(:, :, nodeIDs, 1:kk) = Phi(:, :, :, 1:kk);

        
        splitWait_new(1:n_nodes) = splitWait;



        % Must populate new nodes with relevant information, x, P, weights,
        % qkminus, Pq, Qhat
        x_nodes_splitViable = x_nodes(:, n_nodes_viable, :);
        P_nodes_splitViable = P_nodes(:, :, n_nodes_viable, :);
        Pminus_nodes_splitViable = Pminus_nodes(:, :, n_nodes_viable, :);
        w_nodes_splitViable = w_nodes(n_nodes_viable, :);
        Phi_splitViable = Phi(:, :, n_nodes_viable, :);
        riWindow_nodes_splitViable = rWindow_nodes(:, n_nodes_viable, :);

        

        x_nodes_split = x_nodes_splitViable(:, nodes_to_split, :);
        P_nodes_split = P_nodes_splitViable(:, :, nodes_to_split, :);
        Pminus_nodes_split = Pminus_nodes_splitViable(:, :, nodes_to_split, :);
        w_nodes_split = w_nodes_splitViable(nodes_to_split, :);
        Phi_split = Phi_splitViable(:, :, nodes_to_split, :);
        
        riWindow_split = riWindow_nodes_splitViable(:, nodes_to_split, :);

        H_hist = H_hist_new;
        riPlus_nodes = riPlus_nodes_new;
        riMinus_nodes = riMinus_nodes_new;
        Pri_nodes = Pri_nodes_new;
        PriMinus_nodes = PriMinus_nodes_new;



        newIdx = n_nodes+1; % Add new nodes to the back 

        idxNodestoSplit = find(nodes_to_split == true);
        for i = 1:num_nodes_to_split
            ii_nodeSplit = idxNodestoSplit(i);
            % [x_2, P_2, w_2] = accelEst_BALOO(kk, x_nodes_split(:, i, (timeStart):kk), P_nodes_split(:, :, i, timeStart:kk), w_nodes_split(i, kk-1), riWindow_split(:, i, :), yWindow, time(timeStart:kk), bufferSize, RCF1_Body2_N(:, :, timeStart:kk));
            [x_2, P_2, w_2] = accelEst_BALOO_NRHO(kk, x_nodes_split(:, i, (timeStart):kk), P_nodes_split(:, :, i, timeStart:kk), w_nodes_split(i, kk-1), riWindow_split((1:nM*buffCalcSplit), i, :), yWindow(1:nM*buffCalcSplit), time(timeStart:kk), buffCalcSplit, RCF1_Body2_N(:, :, timeStart:kk), eom, gravParam, beta_markov, accelerometer, gate, Rbody_N, RN_rot,relNavFlag, accelBias, n_elm);
            x_nodes_new(:, newIdx, kk) = x_2(:, 2);
            P_nodes_new(:, :, newIdx, kk) = P_2(:, :, 2);
            Pminus_nodes_new(:, :, newIdx, kk) = P_2(:, :, 2);
            w_nodes_new(newIdx, kk-1) = w_2(2);
            x_nodes_new(:, ii_nodeSplit, kk) = x_2(:, 1);
            P_nodes_new(:, :, ii_nodeSplit, kk) = P_2(:, :, 1);
            Pminus_nodes_new(:, :, ii_nodeSplit, kk) = P_2(:, :, 1);
            w_nodes_new(ii_nodeSplit, kk-1) = w_2(1);
            Phi_new(:, :, newIdx, kk) = Phi_split(:, :, i, kk);



            newIdx = newIdx+1;
            nodes_all = nodes_all+1;
            
            x_nodes_all_new(:, nodes_all, kk) = x_2(:, 2);
            P_nodes_all_new(:, :, nodes_all, kk) = P_2(:, :, 2);
            Pminus_nodes_all_new(:, :, nodes_all, kk) = P_2(:, :, 2);
            w_nodes_all_new(nodes_all, kk-1) = w_2(2);
            Phi_all_new(:, :, nodes_all, kk) = Phi_split(:, :, i, kk);


            nodeIDs = [nodeIDs, nodes_all];

        end


        

        w_nodes_new(:, kk-1) = w_nodes_new(:, kk-1)/(sum(w_nodes_new(:, kk-1)));

        x_nodes = x_nodes_new;
        P_nodes = P_nodes_new;
        Pminus_nodes = Pminus_nodes_new;
        w_nodes = w_nodes_new;
        Phi = Phi_new;

        

        

        rWindow_nodes = rWindow_nodes_new;
        Fwindow_nodes = Fwindow_nodes_new;
        Hwindow_nodes = Hwindow_nodes_new; 

        

        % fullEnough_new = false(1, n_nodes_new);
        % fullEnough_new(1:n_nodes) = fullEnough;
        % fullEnough = fullEnough_new;


        n_nodes = n_nodes_new;
        
        splitWait_new(nodes_to_split) = resConsider*2;
        splitWait = splitWait_new;


        

        x_nodes_all = x_nodes_all_new;
        P_nodes_all = P_nodes_all_new;
        Pminus_nodes_all = Pminus_nodes_all_new;
        w_nodes_all = w_nodes_all_new;
        Phi_all = Phi_all_new;

        


        
        
    end
end
end