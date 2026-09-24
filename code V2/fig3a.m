%% A_sigma_connectivity_er.m: 不同连通概率(p)下的滞后效应对比 (注释切换版)
clear; clc; close all;
addpath('simulations', 'networks');

RUN_SIMULATION = false;  % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)

% ==========================================
% 0. 实验参数配置
% ==========================================
TOPO_TYPE = 'ER';
% P_VALS    = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10];
P_VALS    = [0.03, 0.04, 0.10];

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
num_p = length(P_VALS);

% ==========================================
% 1. 【核心计算区】(由 RUN_SIMULATION 控制)
% ==========================================
if RUN_SIMULATION
    fprintf('\n[SIM] 正在尝试启动 %d 并行仿真任务... \n', num_p);
    simStart = tic;
    for i = 1:num_p
        % 只要这里的 sweep_param 执行过一次，磁盘 results/ 目录下就会有数据
        sweep_param(TOPO_TYPE, P_VALS(i), DYNA);
    end
    fprintf('[SIM] 扫描仿真已完成。总运行时间: %.2f 秒。\n', toc(simStart));
end

% ==========================================
% 2. 【数据加载区】(保证绘图区始终有数据源)
% ==========================================
Results = cell(num_p, 1);
for i = 1:num_p
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, P_VALS(i), DYNA.alpha, DYNA.beta, DYNA.noise);
    data_path = fullfile('results', mat_name);

    if exist(data_path, 'file')
        Results{i} = load(data_path);
        fprintf('[LOAD] 数据已由磁盘载入内存: %s\n', mat_name);
    else
        warning('由于历史数据不存在(%s)，绘图区可能无法正常工作。', mat_name);
    end
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
colors = lines(num_p);

% 绘制所有正向扫描 (Forward)
for i = 1:num_p
    res = Results{i};
    plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $p=%.2f$', P_VALS(i)), ...
        'MarkerIndices', 1:4:length(res.sigma_range));
end

% 绘制所有反向扫描 (Backward)
for i = 1:num_p
    res = Results{i};
    plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $p=%.2f$', P_VALS(i)), ...
        'MarkerIndices', 1:4:length(res.sigma_range));
end

xlabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex'); ylabel('$A(\sigma)$', 'FontSize', 18, 'Interpreter', 'latex');
xlim([10, 25]);
% ylim([0, 140]);
legend('Location', 'NorthWest', 'FontSize', 18, 'NumColumns', 2, 'Interpreter', 'latex');
grid on; set(ax, 'Box', 'on', 'FontSize', 18, 'TickLabelInterpreter', 'latex');

% 图像导出
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'V2', 'manuscirpt', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
% fname = fullfile(plots_dir, 'fig3a.eps');
% exportgraphics(h, fname);
test_plots_dir = fullfile(fileparts(mfilename('fullpath')), 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end
test_fname = fullfile(test_plots_dir, 'fig3a.png');
exportgraphics(h, test_fname, 'Resolution', 300);
fprintf('[DONE] 测试图片已保存: %s\n', test_fname);
