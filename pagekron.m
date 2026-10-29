function [kronBA] = pagekron(B, A)
[n, m, n_nodes] = size(B);
kronBA = zeros(n^2, m^2, n_nodes);
for i = 1:n_nodes
    kronBA(:, :, i) = kron(B(:, :, i), A(:, :, i));
end
end