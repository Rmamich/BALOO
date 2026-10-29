function [Zk, vBar] = meask_adaptivePSD_page(Ak, Bk, Mw, dy, u, bufferSize, flagIncMean, G)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prupose: Generates "measurement" Zk at time step k for adaptive process
% noise estimating psd
% 
% Inputs:
%   Ak  A matrices for time step k ((n x (m*(bufferSize + 1)) x n_nodes)
%   Bk  B matrices for time step k ((n x (n*(bufferSize)) x n_nodes)
%   dy  measurement residual from times k, k-1, ..., k-bufferSize 
%           (m*(bufferSize + 1) x n_nodes)
%   G   mapping matrix from input to state frame (
%   u   control input 
% 
% Outputs: 
%   Zk
% 
% Author: Rachel Mamich
% Date Created: 03/28/23
% Date Moified: 03/29/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n, m, n_nodes] = size(Ak);
m = m/(bufferSize +1);

if nargin < 8
    nU = size(u, 1);
    G = zeros(n*bufferSize, nU, n_nodes);
    d = 3;
    if n > 9
        G_i = [zeros(6, 3); eye(3)];
        G_i = repmat(G_i, 1, 1, n_nodes);
        for i = 1:(bufferSize)
                G((1+n*(i-1)):(i*n), (1+d*(i-1)):(i*d), :) = G_i;
        end
    end
end


Zk = pagemtimes(Ak, reshape(dy, [], 1, n_nodes)) - pagemtimes(Bk, pagemtimes(G, reshape(u, [], 1, n_nodes)));



% dyLong = reshape(dy(:, 1:bufferSize), [], 1);
dyBar = mean(reshape(dy, m, [], n_nodes), 2);




% dyLongBar = mean(dyLong, 2);
% vBar = Mw'*(Mw*Mw')*dyLong;
% wBar = mean(reshape(vBar, n, []), 2);



B_noise = zeros(n*bufferSize, 3, n_nodes);
eye_3x3xn_nodes = repmat(eye(3), 1, 1, n_nodes);
for i = 1:bufferSize
    B_noise((n*i)-2:n*i, :, :) = eye_3x3xn_nodes;
end
MwB = pagemtimes(Mw,B_noise);

vBar = reshape(pagemtimes(pagemldivide(pagemtimes(MwB, 'transpose', MwB, 'none'), eye(3)), pagemtimes(MwB, 'transpose', reshape(dy, [], 1, n_nodes), 'none')), [], n_nodes);

% fprintf('vBar Est %0.6f, %0.6f, %0.6f\n', vBar)
% fprintf('\tyBar Est %0.6f, %0.6f, %0.6f\n', dyBar)

% if flagIncMean
%     Zk = Zk + Ak(:, :, i+1)*dy(:, i+1) + wBar;
% else
%     Zk = Zk + Ak(:, :, i+1)*dy(:, i+1);
% end

end

