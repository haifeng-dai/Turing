%% sigma_alpha_beta.m: 分析层内扩散(alpha)对滞回阈值(sigma)的影响，并对比不同层间耦合强度(beta)
clear; clc; close all;
addpath('simulations', 'networks');

RUN_SIMULATION = false;  % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)

%% 1. 实验参数配置
TOPO_TYPE = 'ER';
DYNA.N = 200;
DYNA.K = 5;
P_VAL_FIXED = 0.03;      % 固定连接概率
DYNA.noise = 0.01;       % 固定噪声强度
DYNA.T_END = 500;        % 扫描点平衡时间
DYNA.steps = 2;          % 采样步数
DYNA.init_perturb = 0.1;

% 搜索范围设置 (Diffusion Ratio sigma)
DYNA.sigma_min = 0;
DYNA.sigma_max = 80;

% 扫描变量定义
ALPHA_LIST = 0.01:0.01:0.15;  % X轴：层内扩散强度
RATIO_LIST = [0.1, 1, 10]; % 这里的 ratio 是倍数参数 (beta = alpha * ratio)

% 数据保存路径
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end
data_path = fullfile(res_dir, sprintf('scan_alpha_ratio_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, P_VAL_FIXED));

%% 2. 核心并行计算区 [64核优化版]
num_a = length(ALPHA_LIST);
num_r = length(RATIO_LIST);

if RUN_SIMULATION
    % === 并行池初始化：充分利用所有64个核 ===
    if isempty(gcp('nocreate'))
        parpool('Threads', 64);
        fprintf('[POOL] 已初始化 64 线程并行池\n');
    else
        pool = gcp;
        fprintf('[POOL] 使用现有并行池 (%d 个worker)\n', pool.NumWorkers);
    end

    [R_idx, A_idx] = meshgrid(1:num_r, 1:num_a);
    R_idx = R_idx(:);
    A_idx = A_idx(:);
    total_tasks = length(R_idx);

    Sigma_F_Results = zeros(total_tasks, 1);
    Sigma_B_Results = zeros(total_tasks, 1);

    % --- A. 进度监控器初始化 ---
    q = parallel.pool.DataQueue;
    hWait = waitbar(0, '正在初始化...', 'Name', '并行任务进度');
    cleanupObj = onCleanup(@() delete(hWait));
    tStart_sim = tic;
    afterEach(q, @(~) nUpdateProgress(total_tasks, hWait, tStart_sim));

    % --- B. 构造层间拉普拉斯算子 ---
    adj_inter = (ones(DYNA.K) - eye(DYNA.K));
    DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

    % --- C. 加载全局种子 (保证一致性) ---
    seed_file = fullfile('results', 'evolution_er_N200_K5_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat');
    if ~exist(seed_file, 'file'), error('找不到全局种子文件，请先运行 pattern_evolution.m！'); end
    tmp_seed = load(seed_file);
    y_universal_seed = tmp_seed.Y(end, :)';

    % --- D. 加载固定拓扑网络 ---
    net_path = fullfile('results', 'topology', 'ER', sprintf('N%d_p%.3f.mat', DYNA.N, P_VAL_FIXED));
    tmp_net = load(net_path);
    L_intra_lib = tmp_net.nets;

    % 1. 粗略扫描 (Rough Scan) - 用于预测引导
    fprintf('\n[SIM] 执行粗略扫描以提取预测线 (Rough Scan)... \n');
    n_rough_a = 8;
    idx_rough_a = round(linspace(1, num_a, n_rough_a));
    alpha_rough = ALPHA_LIST(idx_rough_a);

    [R_rough, A_rough_idxs] = meshgrid(1:num_r, 1:n_rough_a);
    R_rough = R_rough(:);
    A_rough_idxs = A_rough_idxs(:);
    n_rough_tasks = length(R_rough);

    sf_rough_vec = zeros(n_rough_tasks, 1);
    sb_rough_vec = zeros(n_rough_tasks, 1);

    parfor task_id = 1:n_rough_tasks
        c = DYNA;
        cur_alpha = alpha_rough(A_rough_idxs(task_id));
        c.alpha = cur_alpha;
        c.beta  = cur_alpha * RATIO_LIST(R_rough(task_id));
        c.L_intra = L_intra_lib;
        c.y_seed  = y_universal_seed;
        res = find_thresholds(TOPO_TYPE, P_VAL_FIXED, c);
        sf_rough_vec(task_id) = res(1);
        sb_rough_vec(task_id) = res(2);
    end

    sf_rough = reshape(sf_rough_vec, [n_rough_a, num_r]);
    sb_rough = reshape(sb_rough_vec, [n_rough_a, num_r]);

    % 2. 预测插值
    sigma_pred_f = zeros(num_a, num_r);
    sigma_pred_b = zeros(num_a, num_r);
    for r_i = 1:num_r
        sigma_pred_f(:, r_i) = interp1(alpha_rough, sf_rough(:, r_i), ALPHA_LIST, 'pchip', DYNA.sigma_max);
        sigma_pred_b(:, r_i) = interp1(alpha_rough, sb_rough(:, r_i), ALPHA_LIST, 'pchip', DYNA.sigma_max);
    end

    fprintf('\n[SIM] 开始局域精细并行扫描 (Guided Fine Scan): %d 组任务... \n', total_tasks);
    simStart = tic;

    parfor t = 1:total_tasks
        dyna_local = DYNA;
        cur_alpha = ALPHA_LIST(A_idx(t));
        dyna_local.alpha = cur_alpha;
        dyna_local.beta  = cur_alpha * RATIO_LIST(R_idx(t));

        dyna_local.L_intra = L_intra_lib;
        dyna_local.y_seed  = y_universal_seed;

        % 利用预测值缩小搜索区间
        pred_center = max(sigma_pred_f(A_idx(t), R_idx(t)), sigma_pred_b(A_idx(t), R_idx(t)));
        if pred_center > DYNA.sigma_max - 5
            dyna_local.sigma_min = 0;
            dyna_local.sigma_max = DYNA.sigma_max;
        else
            dyna_local.sigma_min = max(0, pred_center - 10);
            dyna_local.sigma_max = min(DYNA.sigma_max, pred_center + 10);
        end

        res = find_thresholds(TOPO_TYPE, P_VAL_FIXED, dyna_local);
        Sigma_F_Results(t) = res(1);
        Sigma_B_Results(t) = res(2);
        send(q, t);
    end

    simTime = toc(simStart);
    fprintf('[DONE] 仿真完成！总耗时: %.2f 秒。\n', simTime);

    SF_Matrix = reshape(Sigma_F_Results, [num_a, num_r]);
    SB_Matrix = reshape(Sigma_B_Results, [num_a, num_r]);
    save(data_path, 'SF_Matrix', 'SB_Matrix', 'ALPHA_LIST', 'RATIO_LIST', 'P_VAL_FIXED');
    fprintf('[SAVE] 结果已保存到: %s\n', data_path);
end

%% 3. 数据可视化
fprintf('\n[VIS] 准备数据可视化...\n');
if exist(data_path, 'file')
    load(data_path);
    fprintf('[LOAD] 已加载结果文件: %s\n', data_path);
else
    fprintf('[WARN] 结果文件不存在！请先运行仿真 (RUN_SIMULATION=true)。\n');
    return;
end

h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.45, 0.5]);
ax = axes('Position', [0.15, 0.15, 0.75, 0.75]); hold on;

colors = jet(num_r);
for j = 1:num_r
    % 正向扫描 (实线)
    plot(ALPHA_LIST, SF_Matrix(:, j), '-', 'Color', colors(j,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('\\eta = %.1f (\\sigma_f)', RATIO_LIST(j)));
    % 反向扫描 (虚线)
    plot(ALPHA_LIST, SB_Matrix(:, j), '--', 'Color', colors(j,:), 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
end

xlabel('Intra Diffusion \alpha', 'FontSize', 14);
ylabel('Critical Diffusion Ratio \sigma', 'FontSize', 14);
title(sprintf('Impact of Alpha on Hysteresis for Selected Ratios \\eta=\\beta/\\alpha (ER, p=%.3f)', P_VAL_FIXED), 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 10, 'NumColumns', 1);
grid off; box on;
set(ax, 'FontSize', 14);

plots_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
out_img = fullfile(plots_dir, sprintf('scan_alpha_ratio_%s_p%.3f.png', lower(TOPO_TYPE), P_VAL_FIXED));
exportgraphics(h, out_img, 'Resolution', 300);
fprintf('[DONE] 绘图已生成: %s\n', out_img);

% ========================================================
function nUpdateProgress(total, hWait, tStart)
persistent count
if isempty(count), count = 0; end
count = count + 1;
progress = count / total;
waitbar(progress, hWait, sprintf('处理进度: %d/%d (%.0f%%) | 耗时: %.1f 秒', ...
    count, total, progress*100, toc(tStart)));
end
