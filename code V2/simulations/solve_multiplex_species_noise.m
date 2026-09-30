function [t, Y, cfg] = solve_multiplex_species_noise(cfg)
% 消融专用多层求解器：噪声组分掩码 [u, v]，原 solve_multiplex 不作修改。
if ~isfield(cfg, 'noise_species_mask') || isempty(cfg.noise_species_mask)
    cfg.noise_species_mask = [1 1];
end
mask = cfg.noise_species_mask(:)';
if ~isnumeric(mask) || ~isreal(mask) || numel(mask) ~= 2 || ...
        any(~isfinite(mask)) || any(mask ~= 0 & mask ~= 1)
    error('solve_multiplex_species_noise:InvalidMask', ...
        'noise_species_mask 必须为 [u_mask, v_mask]，元素只能是 0 或 1。');
end

if ~isfield(cfg, 'dt'), cfg.dt = 0.001; end
dt = cfg.dt;
total_steps = ceil(cfg.T_END / dt);
sample_indices = round(linspace(0, total_steps, cfg.steps));
K = size(cfg.L_inter, 1);
N = cfg.N;
NK = N*K;
if length(cfg.L_intra) > K
    current_L_intra = cfg.L_intra(1:K);
elseif length(cfg.L_intra) < K
    error('网络文件仅包含 %d 层，无法满足当前的 K=%d 配置！', length(cfg.L_intra), K);
else
    current_L_intra = cfg.L_intra;
end

L_INTRA_ALL = blkdiag(current_L_intra{:});
L_INTER_BIG = kron(cfg.L_inter, speye(N));
L_EFF = cfg.alpha*L_INTRA_ALL + cfg.beta*L_INTER_BIG;
if cfg.noise ~= 0 && any(mask)
    [B_edge, edge_coeff] = buildEdgeIncidence( ...
        current_L_intra, cfg.L_inter, N, cfg.alpha, cfg.beta);
end

if isfield(cfg, 'y0') && ~isempty(cfg.y0)
    y = cfg.y0(:);
else
    u0 = 5 + cfg.init_perturb*(rand(NK,1)-0.5);
    v0 = 10 + cfg.init_perturb*(rand(NK,1)-0.5);
    y = [u0; v0];
end
u = y(1:NK);
v = y(NK+1:end);
Y = zeros(cfg.steps, 2*NK);
Y(1,:) = y';
current_sample = 2;

if isfield(cfg, 'noise_seed') && ~isempty(cfg.noise_seed)
    noise_stream = RandStream('mt19937ar', 'Seed', cfg.noise_seed);
    use_local_noise_stream = true;
else
    rng(1024, 'twister');
    use_local_noise_stream = false;
end

A_sum = 0;
A_count = 0;
A_start_step = floor(0.8*total_steps);
prev_A = -1;
detect_convergence = true;
if isfield(cfg, 'detect_convergence')
    detect_convergence = cfg.detect_convergence;
end

for n = 1:total_steps
    f = ((35+16*u-u.^2)/9-v).*u;
    g = (u-(1+0.4*v)).*v;
    Lu = -L_EFF*u;
    Lv = -L_EFF*v;
    if cfg.noise ~= 0 && any(mask)
        state_diffs = B_edge*[u,v];
        if use_local_noise_stream
            dW_edge = sqrt(dt)*randn(noise_stream, length(edge_coeff), 1);
        else
            dW_edge = sqrt(dt)*randn(length(edge_coeff), 1);
        end
        edge_fluxes = edge_coeff.*state_diffs.*dW_edge;
        noise_terms = -cfg.noise*(B_edge'*edge_fluxes);
        u = u+(f+Lu)*dt+mask(1)*noise_terms(:,1);
        v = v+(g+cfg.sigma*Lv)*dt+cfg.sigma*mask(2)*noise_terms(:,2);
    else
        u = u+(f+Lu)*dt;
        v = v+(g+cfg.sigma*Lv)*dt;
    end

    u(u < 1e-6) = 1e-6;
    v(v < 1e-6) = 1e-6;
    if current_sample <= cfg.steps && n == sample_indices(current_sample)
        Y(current_sample,:) = [u',v'];
        current_sample = current_sample+1;
    end

    if mod(n,200) == 0
        A_curr = sqrt(sum((u-5).^2+(v-10).^2)/NK);
        if ~isfinite(A_curr)
            error('数值演化发散 (NaN/Inf)：alpha=%.3f, beta=%.3f, sigma=%.3f, dt=%.4f', ...
                cfg.alpha, cfg.beta, cfg.sigma, cfg.dt);
        end
        if detect_convergence && n > 0.1*total_steps && prev_A > 0
            if abs(A_curr-prev_A) < (A_curr*1e-6+1e-8)
                if current_sample <= cfg.steps
                    Y(current_sample:end,:) = repmat([u',v'], cfg.steps-current_sample+1, 1);
                end
                break;
            end
        end
        prev_A = A_curr;
        if isfield(cfg, 'early_stop') && current_sample > 2
            stop_flag = (strcmpi(cfg.early_stop,'forward') && A_curr > 0.05) || ...
                (strcmpi(cfg.early_stop,'backward') && A_curr < 0.02);
            if stop_flag
                if current_sample <= cfg.steps
                    Y(current_sample:end,:) = repmat([u',v'], cfg.steps-current_sample+1, 1);
                end
                break;
            end
        end
    end
    if n >= A_start_step
        A_sum = A_sum+sqrt(sum((u-5).^2+(v-10).^2)/NK);
        A_count = A_count+1;
    end
end

if A_count > 0
    cfg.A_final = A_sum/A_count;
else
    cfg.A_final = sqrt(sum((u-5).^2+(v-10).^2)/NK);
end
t = linspace(0, cfg.T_END, cfg.steps);
end

function [B_edge, edge_coeff] = buildEdgeIncidence(L_intra, L_inter, N, alpha, beta)
% 每条无向边只生成一次随机增量；u/v 共享增量，掩码只关闭对应组分。
K = length(L_intra);
edge_from = cell(K+1,1);
edge_to = cell(K+1,1);
edge_coeff_cells = cell(K+1,1);
for k = 1:K
    [row,col] = find(triu(-L_intra{k},1));
    offset = (k-1)*N;
    edge_from{k} = offset+row;
    edge_to{k} = offset+col;
    edge_coeff_cells{k} = alpha*ones(length(row),1);
end
[layer_from,layer_to] = find(triu(-L_inter,1));
inter_edge_count = length(layer_from)*N;
inter_from = zeros(inter_edge_count,1);
inter_to = zeros(inter_edge_count,1);
cursor = 1;
node_ids = (1:N)';
for e = 1:length(layer_from)
    idx = cursor:(cursor+N-1);
    inter_from(idx) = (layer_from(e)-1)*N+node_ids;
    inter_to(idx) = (layer_to(e)-1)*N+node_ids;
    cursor = cursor+N;
end
edge_from{end} = inter_from;
edge_to{end} = inter_to;
edge_coeff_cells{end} = beta*ones(inter_edge_count,1);
edge_from = vertcat(edge_from{:});
edge_to = vertcat(edge_to{:});
edge_coeff = vertcat(edge_coeff_cells{:});
edge_count = length(edge_coeff);
B_edge = sparse([(1:edge_count)';(1:edge_count)'], [edge_from;edge_to], ...
    [ones(edge_count,1);-ones(edge_count,1)], edge_count, N*K);
end
