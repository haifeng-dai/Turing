function L = gen_ba(N, m)
% gen_ba: Barabasi-Albert model (Scale-free network)
% N: Total nodes
% m: Edges to attach from new node to existing nodes (m <= m0)

m0 = m + 1; % Initial connected nodes
adj = sparse(N, N);
% Initial clique
adj(1:m0, 1:m0) = ones(m0) - eye(m0);

degs = zeros(N, 1);
degs(1:m0) = m0 - 1;
total_deg = m0 * (m0 - 1);

for i = m0+1:N
    % Preferential attachment probabilities
    probs = degs(1:i-1) / total_deg;

    targets = zeros(1, m);
    available_nodes = 1:i-1;
    current_probs = probs;

    for k = 1:m
        cum_probs = cumsum(current_probs);
        r = rand();
        idx = find(cum_probs >= r, 1, 'first');
        picked_node = available_nodes(idx);
        targets(k) = picked_node;

        % Remove picked node and re-normalize probabilities
        current_probs(idx) = [];
        available_nodes(idx) = [];
        if ~isempty(current_probs)
            current_probs = current_probs / sum(current_probs);
        end
    end

    % Update adjacency matrix
    for t = targets
        adj(i, t) = 1;
        adj(t, i) = 1;
        degs(t) = degs(t) + 1;
    end
    degs(i) = m;
    total_deg = total_deg + 2*m;
end

L = diag(sum(adj, 2)) - adj;
end
