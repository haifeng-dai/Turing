%% A_sigma_er_N.m: 不同网络规模(N)下的滞后效应对比
clear; clc; close all;
addpath('simulations', 'networks');

RUN_SIMULATION = false; % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)

% ==========================================
% 0. 实验参数配置
% ==========================================
TOPO_TYPE = 'ER';
TOPO_PARAM = 0.03;      % 固定层内连接概率 (p)

N_VALS = [100, 500, 1000];
num_n = length(N_VALS);

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
DYNA.sigma_npts = 201;

% 层间结构 (全连通)
adj_inter = ones(DYNA.K) - eye(DYNA.K);
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

% 为每种网络规模生成配置
DYNA_list = cell(num_n, 1);
for i = 1:num_n
    dyna_local = DYNA;
    dyna_local.N = N_VALS(i);
    DYNA_list{i} = dyna_local;
end

% ==========================================
% 1. 【核心计算区】(由 RUN_SIMULATION 控制)
% ==========================================
if RUN_SIMULATION
    fprintf('\n[SIM] 正在尝试启动 %d 并行仿真任务... \n', num_n);
    simStart = tic;
    for i = 1:num_n
        dyna_local = DYNA_list{i};
        fprintf('[SIM] 现正扫描: N=%d ...\n', dyna_local.N);
        sweep_param(TOPO_TYPE, TOPO_PARAM, dyna_local);
    end
    fprintf('[SIM] 扫描仿真已完成。总运行时间: %.2f 秒。\n', toc(simStart));
end

% ==========================================
% 2. 【数据加载区】(保证绘图区始终有数据源)
% ==========================================
Results = cell(num_n, 1);
for i = 1:num_n
    dyna_local = DYNA_list{i};
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), dyna_local.N, dyna_local.K, TOPO_PARAM, dyna_local.alpha, dyna_local.beta, dyna_local.noise);
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
colors = lines(num_n);

% 绘制所有正向扫描 (Forward)
for i = 1:num_n
    if isempty(Results{i}), continue; end
    res = Results{i};
    plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $N=%d$', N_VALS(i)), ...
        'MarkerIndices', 1:4:length(res.sigma_range));
end

% 绘制所有反向扫描 (Backward)
for i = 1:num_n
    if isempty(Results{i}), continue; end
    res = Results{i};
    plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $N=%d$', N_VALS(i)), ...
        'MarkerIndices', 1:4:length(res.sigma_range));
end

xlabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex'); ylabel('$A(\sigma)$', 'FontSize', 18, 'Interpreter', 'latex');
grid on; set(ax, 'Box', 'on', 'FontSize', 18, 'TickLabelInterpreter', 'latex');
legend('Location', 'NorthWest', 'FontSize', 18, 'Interpreter', 'latex', 'NumColumns', 2);

% 图像导出
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
fname = fullfile(plots_dir, 'fig4a.eps');
exportgraphics(h, fname);
fprintf('[DONE] 绘图已更新: %s\n', fname);

