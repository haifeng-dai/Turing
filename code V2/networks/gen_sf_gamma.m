function L = gen_sf_gamma(N, gamma)
% gen_sf_gamma: Generate a scale-free network with a specific scaling exponent.
% N: Total nodes
% gamma: Scaling exponent (Power-law index, usually 2 < gamma < 3)

if gamma <= 2
    error('Scaling exponent gamma must be greater than 2 for a finite mean degree.');
end

% 1. Calculate weights based on the desired gamma
% According to the Static Model: alpha = 1 / (gamma - 1)
alpha = 1 / (gamma - 1);
i_list = (1:N)';
weights = i_list.^(-alpha);
weights = weights / sum(weights); % Normalize to probability distribution

% 2. Calculate target number of edges (maintaining a reasonable density)
% We aim for an average degree of approx 6-8, similar to your other networks
k_avg = 6;
M = round(N * k_avg / 2);

% 3. Generate Edges using Preferential Selection
adj = zeros(N);
edges_count = 0;

% Pre-calculate cumulative distribution for faster sampling
cum_w = cumsum(weights);

while edges_count < M
    % Select two nodes based on their weights
    r1 = rand(); r2 = rand();
    u = find(cum_w >= r1, 1, 'first');
    v = find(cum_w >= r2, 1, 'first');

    % Ensure no self-loops and no redundant edges
    if u ~= v && adj(u, v) == 0
        adj(u, v) = 1;
        adj(v, u) = 1;
        edges_count = edges_count + 1;
    end
end

% 4. Return the Laplacian matrix
L = diag(sum(adj, 2)) - adj;
end
