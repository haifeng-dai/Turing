function [t, Y, cfg] = solve_multiplex(cfg)
% SOLVE_MULTIPLEX 高性能向量化多层反应扩散解算器
% 核心优化：利用稀疏矩阵 Kronecker 积消除层间循环，最小化内存分配开销。

% 1. 超参数与维度设置
if ~isfield(cfg, 'dt'), cfg.dt = 0.001; end % 允许外部覆盖默认步长
dt = cfg.dt;
total_steps = ceil(cfg.T_END / dt);
% 预计算精确采样步索引 (例如第 [1..total_steps] 步中哪些步需要写回矩阵)
sample_indices = round(linspace(0, total_steps, cfg.steps));

% --- 【核心修正】按需截取层数 (Smart Slicing) ---
K = size(cfg.L_inter, 1);
N = cfg.N;
NK = N * K;

% 如果文件加载的网络层数多于需求，只截取前 K 层
if length(cfg.L_intra) > K
    current_L_intra = cfg.L_intra(1:K);
elseif length(cfg.L_intra) < K
    error('网络文件仅包含 %d 层，无法满足当前的 K=%d 配置！', length(cfg.L_intra), K);
else
    current_L_intra = cfg.L_intra;
end

% 2. 预计算全系统算子 (向量化核心)
L_INTRA_ALL = blkdiag(current_L_intra{:});
L_INTER_BIG = kron(cfg.L_inter, speye(N));
% 合成有效全系统算子
L_EFF = cfg.alpha * L_INTRA_ALL + cfg.beta * L_INTER_BIG;

% 为逐边独立噪声预构造全多层网络的关联矩阵。每条无向边仅保留一次，
% 其两端在 B 中分别取 +1 和 -1；edge_coeff 区分层内 alpha 与层间 beta。
% 该结构在整个仿真过程中复用，避免每个时间步重建随机拉普拉斯矩阵。
if cfg.noise ~= 0
    [B_edge, edge_coeff] = build_edge_incidence(current_L_intra, cfg.L_inter, N, ...
        cfg.alpha, cfg.beta);
end

% 3. 初始化状态
if isfield(cfg, 'y0') && ~isempty(cfg.y0)
    y = cfg.y0(:);
else
    % 从稳态周围引入微小扰动以激活图灵不稳定性
    u0 = 5  + cfg.init_perturb * (rand(NK, 1) - 0.5);
    v0 = 10 + cfg.init_perturb * (rand(NK, 1) - 0.5);
    y = [u0; v0];
end

% --- 【内存优化】提前提取局部向量，避免循环内反复切片 ---
u = y(1:NK);
v = y(NK+1:end);

% 4. 预分配输出矩阵
Y = zeros(cfg.steps, 2*NK);
Y(1, :) = y';
current_sample = 2;

% 5. 仿真核心逻辑
% 默认保持原噪声序列；可选 cfg.noise_seed 仅控制噪声流，不改变全局随机状态。
if isfield(cfg, 'noise_seed') && ~isempty(cfg.noise_seed)
    noise_stream = RandStream('mt19937ar', 'Seed', cfg.noise_seed);
    use_local_noise_stream = true;
else
    rng(1024, 'twister');
    use_local_noise_stream = false;
end

% --- 5.1 性能控制变量 ---
A_sum = 0; A_count = 0;
A_start_step = floor(0.8 * total_steps);
prev_A = -1; % 用于收敛探测逻辑

for n = 1:total_steps
    % --- 5.1 计算局部动力学响应 ---
    f = ((35 + 16*u - u.^2)/9 - v) .* u;
    g = (u - (1 + 0.4*v)) .* v;

    % --- 5.2 全系统扩散项计算 ---
    Lu = -L_EFF * u;
    Lv = -L_EFF * v;

    % --- 5.3 逐边独立的乘性噪声 ---
    % 同一条边在 u/v 方程中共用 dW_edge；不同无向边彼此独立。
    if cfg.noise ~= 0
        state_diffs = B_edge * [u, v];
        if use_local_noise_stream
            dW_edge = sqrt(dt) * randn(noise_stream, length(edge_coeff), 1);
        else
            dW_edge = sqrt(dt) * randn(length(edge_coeff), 1);
        end
        edge_fluxes = edge_coeff .* state_diffs .* dW_edge;
        noise_terms = -cfg.noise * (B_edge' * edge_fluxes);

        % 将 y 的更新改为 u/v 的直接累加，只在采样时写回 y。
        u = u + (f + Lu) * dt + noise_terms(:, 1);
        v = v + (g + cfg.sigma * Lv) * dt + cfg.sigma * noise_terms(:, 2);
    else
        u = u + (f + Lu) * dt;
        v = v + (g + cfg.sigma * Lv) * dt;
    end

    % 极限截断 (非负性约束向量化)
    u(u < 1e-6) = 1e-6;
    v(v < 1e-6) = 1e-6;

    % --- 5.5 结果采样与保存 (精确索引方案 B) ---
    if current_sample <= cfg.steps && n == sample_indices(current_sample)
        Y(current_sample, :) = [u', v'];
        current_sample = current_sample + 1;
    end

    % --- 5.6 增强型提前退出逻辑 (发散保护 + 收敛探测 + 早停加速) ---
    if mod(n, 200) == 0
        A_curr = sqrt(sum((u - 5).^2 + (v - 10).^2) / NK);

        % A. 发散检测 (Robustness: 拦截由于步长不稳导致的 NaN)
        if isnan(A_curr) || isinf(A_curr)
            error('数值演化发散 (NaN/Inf)：alpha=%.3f, beta=%.3f, sigma=%.3f, dt=%.4f', ...
                cfg.alpha, cfg.beta, cfg.sigma, cfg.dt);
        end

        % B. 自动稳态判定 (Convergence Analysis)
        % 如果已经运行超过 10%，且序参数不再剧烈波动，提前判定为已平衡
        if n > 0.1 * total_steps && prev_A > 0
            if abs(A_curr - prev_A) < (A_curr * 1e-6 + 1e-8)
                if current_sample <= cfg.steps
                    Y(current_sample:end, :) = repmat([u', v'], cfg.steps - current_sample + 1, 1);
                end
                break;
            end
        end
        prev_A = A_curr;

        % C. 早停加速 (与 find_thresholds.m 的归一化序参数阈值对齐)
        if isfield(cfg, 'early_stop') && current_sample > 2
            stop_flag = false;
            if strcmpi(cfg.early_stop, 'forward') && A_curr > 0.05
                stop_flag = true;
            elseif strcmpi(cfg.early_stop, 'backward') && A_curr < 0.02
                stop_flag = true;
            end
            if stop_flag
                if current_sample <= cfg.steps
                    Y(current_sample:end, :) = repmat([u', v'], cfg.steps - current_sample + 1, 1);
                end
                break;
            end
        end
    end

    % --- 5.7 时间窗口平均计算 (Time-Averaging) ---
    if n >= A_start_step
        A_sum = A_sum + sqrt(sum((u - 5).^2 + (v - 10).^2) / NK);
        A_count = A_count + 1;
    end
end

% 计算出口平均序参数 (基于最后 10% 时间窗口)
if A_count > 0
    cfg.A_final = A_sum / A_count;
else
    cfg.A_final = sqrt(sum((u - 5).^2 + (v - 10).^2) / NK);
end
t = linspace(0, cfg.T_END, cfg.steps);
end

function [B_edge, edge_coeff] = build_edge_incidence(L_intra, L_inter, N, alpha, beta)
% BUILD_EDGE_INCIDENCE 由拉普拉斯矩阵的上三角非零元构造一次无向边关联矩阵。
K = length(L_intra);
edge_from = cell(K + 1, 1);
edge_to = cell(K + 1, 1);
edge_coeff_cells = cell(K + 1, 1);

for k = 1:K
    [row, col] = find(triu(-L_intra{k}, 1));
    offset = (k - 1) * N;
    edge_from{k} = offset + row;
    edge_to{k} = offset + col;
    edge_coeff_cells{k} = alpha * ones(length(row), 1);
end

% 层间边连接相同 node ID 的两个层；L_inter 的上三角给出每条无向层间边。
[layer_from, layer_to] = find(triu(-L_inter, 1));
inter_edge_count = length(layer_from) * N;
inter_from = zeros(inter_edge_count, 1);
inter_to = zeros(inter_edge_count, 1);
cursor = 1;
node_ids = (1:N)';
for e = 1:length(layer_from)
    idx = cursor:(cursor + N - 1);
    inter_from(idx) = (layer_from(e) - 1) * N + node_ids;
    inter_to(idx) = (layer_to(e) - 1) * N + node_ids;
    cursor = cursor + N;
end
edge_from{end} = inter_from;
edge_to{end} = inter_to;
edge_coeff_cells{end} = beta * ones(inter_edge_count, 1);

edge_from = vertcat(edge_from{:});
edge_to = vertcat(edge_to{:});
edge_coeff = vertcat(edge_coeff_cells{:});
edge_count = length(edge_coeff);
B_edge = sparse([(1:edge_count)'; (1:edge_count)'], [edge_from; edge_to], ...
    [ones(edge_count, 1); -ones(edge_count, 1)], edge_count, N * K);
end
