%% A_sigma_connectivity_ws.m: 不同重连概率(p)下的WS网络滞后效应对比
clear; clc; close all;
addpath('simulations', 'networks');


% ==========================================
% 0. 实验参数配置
% ==========================================
TOPO_TYPE = 'WS';
% K_VALS    = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20];  % WS 网络的平均度 K (degree)
P_REWIRE  = 0.01;
K_VALS    = [6, 8, 20];  % WS 网络的平均度 K (degree)

DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta  = 0.1 * DYNA.alpha;
DYNA.noise = 0.1;
DYNA.T_END = 200;
DYNA.steps = 2;
DYNA.init_perturb = 0.1;

% 扫描范围
DYNA.sigma_min  = 10;
DYNA.sigma_max  = 30;
DYNA.sigma_npts = 76;   % 对应 0.1 的步长

% 层间结构 (全连通)
adj_inter = ones(DYNA.K) - eye(DYNA.K);
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;
num_k = length(K_VALS);

% ==========================================
% 1. 【核心计算区】(缺失时自动生成种子与扫描)
% ==========================================
results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

FORCE_RERUN = false; % true: 强制重算并覆盖旧数据；false: 优先读取已有缓存

% 检查结果是否缺失或强制重算
missing_result = false(1, num_k);
for i = 1:num_k
    mat_name = sprintf('hysteresis_%s_N%d_K%d_pr%.2f_mK%d_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, P_REWIRE, K_VALS(i), DYNA.alpha, DYNA.beta, DYNA.noise);
    missing_result(i) = FORCE_RERUN || ~exist(fullfile(results_dir, mat_name), 'file');
end

if any(missing_result)
    % 检查并自愈标准种子
    seed_name = sprintf('evolution_er_N%d_K%d_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', ...
        DYNA.N, DYNA.K);
    seed_file = fullfile(results_dir, seed_name);
    if ~exist(seed_file, 'file')
        fprintf('[INFO] 未检测到种子文件，正在生成参考斑图种子...\n');
        dyna_seed = DYNA;
        dyna_seed.sigma = 100;
        dyna_seed.noise = 0;
        pattern_evolution(dyna_seed, 'ER', 0.030, false);
    end

    fprintf('\n[SIM] 正在尝试启动 %d 组缺失的仿真任务 (WS Network)... \n', sum(missing_result));
    simStart = tic;
    for i = find(missing_result)
        dyna_temp = DYNA;
        dyna_temp.ws_k = K_VALS(i);
        dyna_temp.ws_p = P_REWIRE;
        sweep_param_ws(TOPO_TYPE, K_VALS(i), P_REWIRE, dyna_temp);
    end
    fprintf('[SIM] 扫描仿真已完成。总运行时间: %.2f 秒。\n', toc(simStart));
end

% ==========================================
% 2. 【数据加载区】
% ==========================================
Results = cell(num_k, 1);
for i = 1:num_k
    mat_name = sprintf('hysteresis_%s_N%d_K%d_pr%.2f_mK%d_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, P_REWIRE, K_VALS(i), DYNA.alpha, DYNA.beta, DYNA.noise);
    data_path = fullfile(results_dir, mat_name);
    Results{i} = load(data_path);
    fprintf('[LOAD] 数据已由磁盘载入内存: %s\n', mat_name);
end

% ==========================================
% 3. 【绘图展示区】(对比 A vs Sigma)
% ==========================================
% 检查是否存在绘图数据
if isempty(Results) || isempty(Results{1})
    error('关键缺失：未找到任何可绘制的结果，请先运行 Section 1。');
end

h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.45]);
ax = axes('Position', [0.18, 0.18, 0.75, 0.75]); hold on;
colors = lines(num_k);

% 绘制所有正向扫描 (Forward)
for i = 1:num_k
    if ~isempty(Results{i})
        res = Results{i};
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $k=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:4:length(res.sigma_range));
    end
end

% 绘制所有反向扫描 (Backward)
for i = 1:num_k
    if ~isempty(Results{i})
        res = Results{i};
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $k=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:4:length(res.sigma_range));
    end
end

% title('Hysteresis Loop across WS Network Rewiring Probability (p)', 'FontSize', 12);
xlabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex'); ylabel('$A(\sigma)$', 'FontSize', 18, 'Interpreter', 'latex');
xlim([10, 30]);
% ylim([0, 140]);
legend('Location', 'NorthWest', 'FontSize', 18, 'Interpreter', 'latex', 'NumColumns', 2);
grid on; set(ax, 'Box', 'on', 'FontSize', 18, 'TickLabelInterpreter', 'latex');

% 图像导出
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'V2', 'manuscirpt', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
% fname = fullfile(plots_dir, 'fig3b.eps');
% exportgraphics(h, fname);
test_plots_dir = fullfile(fileparts(mfilename('fullpath')), 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end
test_fname = fullfile(test_plots_dir, 'fig3b.png');
exportgraphics(h, test_fname, 'Resolution', 300);
fprintf('[DONE] 测试图片已保存: %s\n', test_fname);
