function L = gen_er(N, p)
% GEN_ER Generates a Laplacian for an Erdos-Renyi random network.
%   N: nodes, p: connection probability
adj = rand(N) < p;
adj = double(triu(adj, 1));
adj = adj + adj';
L = diag(sum(adj, 2)) - adj;
end
