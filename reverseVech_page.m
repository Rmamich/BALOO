function Qfull = reverseVech_page(QUTVec)
    [x, ~, n_nodes] = size(QUTVec);
    n = sqrt(2*x + 1/4)-1/2;
    
    Qfull = zeros(n, n, n_nodes)*QUTVec(1);
    indexVec = 1:(n_nodes*n^2);
    indexMat = reshape(indexVec, n, n, n_nodes);
    indexLT = reshape(vech_page(indexMat), [], 1);

    Qfull(indexLT) = reshape(QUTVec, [], 1);
    
    indexMat = pagetranspose(reshape(indexVec, n, n, n_nodes));
    indexUT = vech_page(indexMat);
    
    Qfull(indexUT) = reshape(QUTVec, [], 1);
    
end