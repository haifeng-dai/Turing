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
% --- 5.0 预生成全过程的公共布朗增量 (强制锁定随机种子以确保噪声一致性) ---
rng(1024, 'twister');
dW_all = sqrt(dt) * randn(total_steps, 1);

% --- 5.1 性能控制变量 ---
A_sum = 0; A_count = 0;
A_start_step = floor(0.8 * total_steps);
prev_A = -1; % 用于收敛探测逻辑

for n = 1:total_steps
    % --- 5.1 提取当前状态并领取随机增量 ---
    dW = dW_all(n);

    % --- 5.2 计算局部动力学响应 (公式级硬编码以获得巅峰性能) ---
    f = ((35 + 16*u - u.^2)/9 - v) .* u;
    g = (u - (1 + 0.4*v)) .* v;

    % --- 5.3 全系统扩散项计算 ---
    Lu = -L_EFF * u;
    Lv = -L_EFF * v;

    % --- 5.4 状态原地更新 (优化：彻底消除循环内下标索引开销) ---
    % 将 y 的更新改为 u/v 的直接累加，只在采样时写回 y
    u = u + (f + Lu) * dt + (cfg.noise * Lu) * dW;
    v = v + (g + cfg.sigma * Lv) * dt + (cfg.noise * cfg.sigma * Lv) * dW;

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
        A_curr = sqrt(sum((u - 5).^2 + (v - 10).^2));

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

        % C. 早停加速 (与 find_thresholds.m 判定准则 0.05 强对齐)
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
        A_sum = A_sum + sqrt(sum((u - 5).^2 + (v - 10).^2)); % 基于图片公式总结
        A_count = A_count + 1;
    end
end

% 计算出口平均序参数 (基于最后 10% 时间窗口)
if A_count > 0
    cfg.A_final = A_sum / A_count;
else
    cfg.A_final = sqrt(sum((u - 5).^2 + (v - 10).^2));
end
t = linspace(0, cfg.T_END, cfg.steps);
end
