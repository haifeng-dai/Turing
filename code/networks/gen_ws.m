function L = gen_ws(N, K_neighbor, p_rewire)
% GEN_WS Generates a Laplacian for a Watts-Strogatz small-world network.
%   N: nodes, K_neighbor: mean degree (even), p_rewire: rewiring prob.

adj = zeros(N);
for i = 1:N
    for j = 1:K_neighbor/2
        target = mod(i+j-1, N) + 1;
        adj(i, target) = 1;
        adj(target, i) = 1;
    end
end

for i = 1:N
    for j = 1:K_neighbor/2
        if rand() < p_rewire
            target_old = mod(i+j-1, N) + 1;
            adj(i, target_old) = 0;
            adj(target_old, i) = 0;

            % New random target
            candidates = find(~adj(i, :));
            candidates(candidates == i) = [];
            if ~isempty(candidates)
                new_target = candidates(randi(length(candidates)));
                adj(i, new_target) = 1;
                adj(new_target, i) = 1;
            end
        end
    end
end

L = diag(sum(adj, 2)) - adj;
end
