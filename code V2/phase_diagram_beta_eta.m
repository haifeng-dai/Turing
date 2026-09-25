%% phase_diagram_beta_eta.m: 层间耦合(beta)与噪声(eta)二维相图分析 (方案一: 热力云图 + 临界相边界)
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'simulations'), fullfile(script_dir, 'networks'));

RUN_SIMULATION = false;  % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)

%% 1. 实验参数配置
TOPO_TYPE = 'ER';
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
P_VAL_FIXED = 0.01;      % 固定连接概率
DYNA.T_END = 500;
DYNA.steps = 2;
DYNA.init_perturb = 0.1;

DYNA.sigma_min = 0;
DYNA.sigma_max = 80;

% --- 扫描变量定义 (兼顾分辨率与 64 核算力) ---
% beta 覆盖 0.0025 到 0.5 (即 beta/alpha 约为 0.05 到 10)
BETA_FACTORS = [0.05, 0.1, 0.2, 0.5, 1, 2, 4, 6, 8, 10];
BETA_LIST    = BETA_FACTORS * DYNA.alpha;   % 层间耦合强度
ETA_LIST     = linspace(0, 1.5, 31);        % 噪声强度 (步长 0.05)

% 数据保存路径
res_dir = fullfile(script_dir, 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end
data_path = fullfile(res_dir, sprintf('phase_diagram_beta_eta_%s_N%d_K%d_a%.3f_p%.3f.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, P_VAL_FIXED));

num_e = length(ETA_LIST);
num_b = length(BETA_LIST);

%% 2. 核心并行计算区
if RUN_SIMULATION
    if isempty(gcp('nocreate'))
        parpool('Threads', 64);
        fprintf('[POOL] 已初始化 64 线程并行池\n');
    else
        pool = gcp;
        fprintf('[POOL] 使用现有并行池 (%d 个worker)\n', pool.NumWorkers);
    end

    [B_idx, E_idx] = meshgrid(1:num_b, 1:num_e);
    B_idx = B_idx(:);
    E_idx = E_idx(:);
    total_tasks = length(B_idx);

    Sigma_F_Results = zeros(total_tasks, 1);
    Sigma_B_Results = zeros(total_tasks, 1);

    q = parallel.pool.DataQueue;
    hWait = waitbar(0, '正在初始化...', 'Name', '并行任务进度');
    cleanupObj = onCleanup(@() delete(hWait));
    tStart_sim = tic;
    afterEach(q, @(~) nUpdateProgress(total_tasks, hWait, tStart_sim));

    % A. 构造层间算子
    adj_inter = (ones(DYNA.K) - eye(DYNA.K));
    DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

    % B. 加载全局种子
    seed_file = fullfile(res_dir, sprintf('evolution_er_N%d_K%d_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', ...
        DYNA.N, DYNA.K));
    if ~exist(seed_file, 'file')
        dyna_seed = DYNA;
        dyna_seed.sigma = 100;
        dyna_seed.noise = 0;
        pattern_evolution(dyna_seed, TOPO_TYPE, 0.030, false);
    end
    tmp_seed = load(seed_file);
    y_universal_seed = tmp_seed.Y(end, :)';

    % C. 加载拓扑网络
    net_path = fullfile(res_dir, 'topology', 'ER', sprintf('N%d_p%.3f.mat', DYNA.N, P_VAL_FIXED));
    tmp_net = load(net_path);
    L_intra_lib = tmp_net.nets;

    % 1. 粗扫描预测引导 (Rough Scan)
    fprintf('\n[SIM] 执行粗略扫描以提取预测线 (Rough Scan)... \n');
    n_rough_b = min(5, num_b);
    idx_rough_b = round(linspace(1, num_b, n_rough_b));
    beta_rough = BETA_LIST(idx_rough_b);

    [B_rough, E_rough] = meshgrid(1:n_rough_b, 1:num_e);
    B_rough = B_rough(:);
    E_rough = E_rough(:);
    n_rough_tasks = length(E_rough);

    sf_rough_vec = zeros(n_rough_tasks, 1);
    sb_rough_vec = zeros(n_rough_tasks, 1);

    parfor task_id = 1:n_rough_tasks
        c = DYNA;
        c.noise = ETA_LIST(E_rough(task_id));
        c.beta  = beta_rough(B_rough(task_id));
        c.L_intra = L_intra_lib;
        c.y_seed  = y_universal_seed;
        res = find_thresholds(TOPO_TYPE, P_VAL_FIXED, c);
        sf_rough_vec(task_id) = res(1);
        sb_rough_vec(task_id) = res(2);
    end

    sf_rough = reshape(sf_rough_vec, [num_e, n_rough_b]);
    sb_rough = reshape(sb_rough_vec, [num_e, n_rough_b]);

    % 2. 插值预测引导
    sigma_pred_f = zeros(num_e, num_b);
    sigma_pred_b = zeros(num_e, num_b);
    for e_i = 1:num_e
        sigma_pred_f(e_i, :) = interp1(beta_rough, sf_rough(e_i, :), BETA_LIST, 'pchip', DYNA.sigma_max);
        sigma_pred_b(e_i, :) = interp1(beta_rough, sb_rough(e_i, :), BETA_LIST, 'pchip', DYNA.sigma_max);
    end

    % 3. 精细扫描
    fprintf('\n[SIM] 开始局域精细并行扫描 (Guided Fine Scan): %d 组任务... \n', total_tasks);
    simStart = tic;

    parfor t = 1:total_tasks
        dyna_local = DYNA;
        dyna_local.beta  = BETA_LIST(B_idx(t));
        dyna_local.noise = ETA_LIST(E_idx(t));
        dyna_local.L_intra = L_intra_lib;
        dyna_local.y_seed  = y_universal_seed;

        pred_center = max(sigma_pred_f(E_idx(t), B_idx(t)), sigma_pred_b(E_idx(t), B_idx(t)));
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

    SF_Matrix = reshape(Sigma_F_Results, [num_e, num_b]);
    SB_Matrix = reshape(Sigma_B_Results, [num_e, num_b]);
    save(data_path, 'SF_Matrix', 'SB_Matrix', 'ETA_LIST', 'BETA_LIST', 'P_VAL_FIXED', 'BETA_FACTORS');
    fprintf('[SAVE] 仿真数据已保存至: %s\n', data_path);
end

%% ============================================================
%% 3. 【方案一】二维相图可视化 (热力云图 + 临界相边界等高线)
%% ============================================================
fprintf('\n[VIS] 正在生成二维相图...\n');
if ~exist(data_path, 'file')
    error('未找到数据文件: %s，请将 RUN_SIMULATION 设为 true 运行仿真。', data_path);
end
data = load(data_path);
ETA_LIST  = data.ETA_LIST;
BETA_LIST = data.BETA_LIST;
SF_Matrix = data.SF_Matrix;
SB_Matrix = data.SB_Matrix;

% 计算滞后区间宽度 Delta sigma
Delta_Sigma = SF_Matrix - SB_Matrix;
Delta_Sigma(Delta_Sigma < 0) = 0; % 数值截断

% 网格平滑插值 (提升相图渲染质感)
[B_grid, E_grid] = meshgrid(BETA_LIST, ETA_LIST);
b_fine = linspace(min(BETA_LIST), max(BETA_LIST), 200);
e_fine = linspace(min(ETA_LIST), max(ETA_LIST), 200);
[B_fine, E_fine] = meshgrid(b_fine, e_fine);
Delta_fine = interp2(B_grid, E_grid, Delta_Sigma, B_fine, E_fine, 'spline');
Delta_fine(Delta_fine < 0) = 0;

% 创建图窗
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.48, 0.46]);
ax = axes('Position', [0.15, 0.17, 0.68, 0.74]);
hold(ax, 'on');

% 1. 绘制热力云图 (pcolor 平滑过渡)
p = pcolor(ax, B_fine, E_fine, Delta_fine);
set(p, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax, parula);
cb = colorbar(ax, 'Position', [0.85, 0.17, 0.03, 0.74]);
ylabel(cb, 'Hysteresis Width $\Delta\sigma = \sigma_c^f - \sigma_c^b$', ...
    'FontSize', 14, 'Interpreter', 'latex');
set(cb, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

% 2. 叠加等高线 (相变分界线: Delta_sigma 临界消亡线)
tol_boundary = 0.10; % 滞后消失判断准则
[C, hc] = contour(ax, B_fine, E_fine, Delta_fine, [tol_boundary, tol_boundary], ...
    'w--', 'LineWidth', 2.5, 'DisplayName', 'Phase Boundary ($\Delta\sigma \approx 0$)');

% 3. 区域文字标注 (双稳态 vs 单稳态)
text(ax, mean(BETA_LIST)*0.6, 0.2, 'Bistable Region', ...
    'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold', ...
    'Interpreter', 'latex', 'Color', 'w');
text(ax, mean(BETA_LIST)*0.8, 1.2, 'Monostable Region', ...
    'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold', ...
    'Interpreter', 'latex', 'Color', [0.95 0.95 0.95]);

% 4. 坐标轴与排版精修
xlabel(ax, 'Inter-layer Coupling Strength $\beta$', 'FontSize', 15, 'Interpreter', 'latex');
ylabel(ax, 'Noise Intensity $\eta$', 'FontSize', 15, 'Interpreter', 'latex');
xlim(ax, [min(BETA_LIST), max(BETA_LIST)]);
ylim(ax, [min(ETA_LIST), max(ETA_LIST)]);
box(ax, 'on');
set(ax, 'FontSize', 14, 'Layer', 'top', 'TickLabelInterpreter', 'latex');

% 导出高清图
test_plots_dir = fullfile(script_dir, 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end
out_png = fullfile(test_plots_dir, 'phase_diagram_beta_eta.png');
exportgraphics(h, out_png, 'Resolution', 300);
fprintf('[DONE] 二维相图已保存至: %s\n', out_png);

% ========================================================
function nUpdateProgress(total, hWait, tStart)
persistent count
if isempty(count), count = 0; end
count = count + 1;
progress = count / total;
waitbar(progress, hWait, sprintf('处理进度: %d/%d (%.0f%%) | 耗时: %.1f 秒', ...
    count, total, progress*100, toc(tStart)));
end
