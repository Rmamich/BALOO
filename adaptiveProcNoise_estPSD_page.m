function [qkplus, Pkplus, Awindow, Bwindow, ZWindow, vBar] = adaptiveProcNoise_estPSD_page(qkminus, Pkminus, Xk, R, yWindow, uWindow, FWindow, Hwindow, Awindow, Bwindow, ZWindow, bufferSize, flagBufferReqMet, QW, RW, beta, dt, flagQ, flagIncMean, GWindow )
flagSafeAB = (sum(isnan(yWindow)) == 0);
qkplus = qkminus;
Pkplus = Pkminus;
[nA, ~, n_nodes, ~] = size(FWindow);
vBar = zeros(3, n_nodes);

buffReqNotMet = reshape((sum(any(Awindow(:, :, :)~=0)) == 0), [], 1);

flagBufferReqMet(buffReqNotMet) = false;

if any(flagSafeAB)
    % n_nodes_build = sum(flagSafeAB);
    Awindow_build = Awindow(:, :, flagSafeAB);
    Bwindow_build = Bwindow(:, :, flagSafeAB);
    ZWindow_build = ZWindow(:, flagSafeAB);
    FWindow_build = FWindow(:, :, flagSafeAB, :);
    HWindow_build = Hwindow(:, :, flagSafeAB);
    yWindow_build = yWindow(:, flagSafeAB);
    uWindow_build = uWindow(:, flagSafeAB);

    
    endFirst = nA*(nA+1)/2;
    nA2 = nA^2;

    [A, B, Mw] = constructAB_page(FWindow_build, HWindow_build, bufferSize);

    % Shift A and B windows back
    Awindow_build((1+nA):(nA*(bufferSize+1)), :, :) = Awindow_build(1:(nA*(bufferSize)), :, :);
    Bwindow_build((1+nA):(nA*(bufferSize)), :, :) = Bwindow_build(1:(nA*(bufferSize-1)), :, :);
    % Overwrite first slot
    Awindow_build(1:nA, :, :) = A;
    Bwindow_build(1:nA, :, :) = B;

    if nargin > 19 % not really tested
        [Zk, vBar_build] = meask_adaptivePSD_page(A, B, Mw, yWindow_build, uWindow_build, bufferSize, flagIncMean, GWindow);
        yWindowBar = mean(yWindow_build, 2);
    else
        [Zk, vBar_build] = meask_adaptivePSD_page(A, B, Mw, yWindow_build, uWindow_build, bufferSize, flagIncMean);
        yWindowBar_build = mean(yWindow_build, 2);
    end
    ZWindow_build((1+nA):(nA*(bufferSize+1)), :) = ZWindow_build(1:(nA*(bufferSize)), :);
    ZWindow_build(1:nA, :) = Zk;


    Awindow(:, :, flagSafeAB) = Awindow_build;
    Bwindow(:, :, flagSafeAB) = Bwindow_build;
    ZWindow(:, flagSafeAB) = ZWindow_build;
    % FWindow(:, :, flagSafeAB, :) = FWindow_build;
    % Hwindow(:, :, flagSafeAB) = HWindow_build;
    % yWindow(:, flagSafeAB) = yWindow_build;
    % uWindow(:, flagSafeAB) = uWindow_build;
    
    vBar(:, flagSafeAB) = vBar_build;

    

    idxFull = (flagBufferReqMet==true);
    n_nodes_update = sum(idxFull);
    Awindow_update = Awindow(:, :, idxFull);
    Bwindow_update = Bwindow(:, :, idxFull);
    ZWindow_update = ZWindow(:, idxFull);
    qkminus_update = qkminus(:, idxFull);
    Pkminus_update = Pkminus(:, :, idxFull);

    % if sum(idxFull == 2)
    % fndchadios = -0;
    % end


    if any(flagBufferReqMet)
        [H, ABlock] = formMeasMatrix_variableBufferSize_page(Awindow_update, Bwindow_update, bufferSize);
        b = zeros((bufferSize)*nA2+ endFirst, 1, n_nodes_update);
        b(1:endFirst, :, :) = kronh_page(reshape(ZWindow_update(1:nA, :), [], 1, n_nodes_update));
        for i = 2:(bufferSize + 1)
            % b((nA2*(i-2)+1+ endFirst):(nA2*(i-1) + endFirst), :, :) = pagekron(ZWindow(:, i), ZWindow(:, 1));
            b((nA2*(i-2)+1+ endFirst):(nA2*(i-1) + endFirst), :, :) = pagekron(reshape(ZWindow_update((1+nA*(i-1)):(nA*i), :), [], 1, n_nodes_update), reshape(ZWindow_update(1:nA, :), [], 1, n_nodes_update));
        end
        % TODO: Comment this back in
        b = b - pagemtimes(ABlock, vech(R)); % Remove contribution from set R
        if ~contains(flagQ, 'Inf') % Haven't really touched this
            dqds = qkminus_update;
            dqds = diag(dqds);
            H = pagemtimes(H, pagemtimes(Xk, dqds));
            [qkplus, Pkplus] = recursiveLS_psd(H, b, qkminus_update, Pkminus_update, QW, RW, beta, dt);
        else
            % The inputs and outputs of this function are in the information
            % form and will require post processing in the main script to
            % obtain the workable values
            qkminus_update = reshape(qkminus_update, [], 1, n_nodes_update);
            dqds = pagemtimes(pagemldivide(Pkminus_update, eye(3)), qkminus_update);
            dqds_mat = zeros(3, 3, n_nodes_update);
            for i = 1:n_nodes_update
                dqds_mat(:, :, i) = diag(dqds(:, :, i));
            end
            dqds = dqds_mat;
            H = pagemtimes(H, pagemtimes(Xk(:, :, idxFull), dqds));
            [qkplus_update, Pkplus_update] = recursiveLS_psd_information_page(H, b, qkminus_update, Pkminus_update, QW, RW, beta, dt);
            qkplus(:, idxFull) = qkplus_update;
            Pkplus(:, :, idxFull) = Pkplus_update;
        end

    end

end

end




