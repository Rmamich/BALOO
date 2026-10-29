n_nodes = 3;

A = rand(3,3,n_nodes);
B = rand(3,3,n_nodes);
C = zeros(3,3,n_nodes);  % preallocate the results
 
tic
for i=1:n_nodes
    C(:,:,i) = A(:,:,i) * B(:,:,i);
end
loopTime = toc;

tic
C1 = pagemtimes(A,B);
pagedTime = toc;

tic 
C(:,:,1) = A(:,:,1) * B(:,:,1);
singleMult = toc;



%%

x_nodes = [1, 2, 3; 4, 5, 6; 7, 8, 9; 1, 2, 3; 4, 5, 6; 7, 8, 9];
x_nodes(:, 1) = x_nodes(:, 1)*2;
x_nodes(:, 3) = x_nodes(:, 3)*5;

Pminus_nodes = eye(6);
Pminus_nodes(1:3, 1:3) = [1, 2, 3; 2, 4, 5; 3, 5, 6];
Pminus_nodes(4:6, 4:6) = [1, 2, 3; 2, 4, 5; 3, 5, 6];

Pminus_nodes = repmat(Pminus_nodes, 1, 1, n_nodes);
Pminus_nodes(:, :, 1) = Pminus_nodes(:, :, 1)*4;
Pminus_nodes(:, :, 3) = Pminus_nodes(:, :, 3)*3;


% w_nodes = ones(n_nodes, 5)/n_nodes;
w_nodes = [1; 5; 3];
w_nodes = w_nodes/sum(w_nodes);

x = zeros(6, 1);
P = zeros(6, 1);

for k = 1:n_nodes
    x = x + w_nodes(k, 1)*x_nodes(:, k);
end
for k = 1:n_nodes
    P = P + w_nodes(k, 1)*(Pminus_nodes(:, :, k) + (x_nodes(:, k) - x)*(x_nodes(:, k) - x)');
end


x1 = x_nodes*w_nodes(:, 1); % check
e_minus = x_nodes - x1;
e_minus = reshape(e_minus, [], 1, n_nodes);
Pminus= sum(pagemtimes(reshape(w_nodes(:, 1), 1, 1, []), Pminus_nodes + pagemtimes(e_minus, 'none', e_minus, 'transpose')), 3);


[x2, P2] = outputGMM(x_nodes, Pminus_nodes, w_nodes, n_nodes);
