%% A_sigma_er_K.m: 不同网络层数(K)下的ER网络滞后效应对比
clear; clc; close all;
addpath('simulations', 'networks');

RUN_SIMULATION = false; % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)

% ==========================================
% 0. 实验参数配置
% ==========================================
TOPO_TYPE = 'ER';
TOPO_PARAM = 0.03;      % 固定层内连接概率 (p)

K_VALS = 2:4:10;        % 对比层数 (2, 6, 10)
num_k  = length(K_VALS);

DYNA.N = 200;           % 固定网络规模
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

% 为每种层数生成配置
DYNA_list = cell(num_k, 1);
for i = 1:num_k
    dyna_local = DYNA;
    dyna_local.K = K_VALS(i);
    % 构造层间结构 (全连通)
    adj_inter = ones(dyna_local.K) - eye(dyna_local.K);
    dyna_local.L_inter = diag(sum(adj_inter, 2)) - adj_inter;
    DYNA_list{i} = dyna_local;
end

% ==========================================
% 1. 【核心计算区】(由 RUN_SIMULATION 控制)
% ==========================================
if RUN_SIMULATION
    fprintf('\n[SIM] 正在尝试启动 %d 并行仿真任务 (ER Network Layers)... \n', num_k);
    simStart = tic;
    for i = 1:num_k
        dyna_local = DYNA_list{i};
        fprintf('[SIM] 现正扫描: K=%d ...\n', dyna_local.K);
        sweep_param(TOPO_TYPE, TOPO_PARAM, dyna_local);
    end
    fprintf('[SIM] 扫描仿真已完成。总运行时间: %.2f 秒。\n', toc(simStart));
end

% ==========================================
% 2. 【数据加载区】(保证绘图区始终有数据源)
% ==========================================
Results = cell(num_k, 1);
for i = 1:num_k
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
if isempty(Results) || (num_k > 0 && isempty(Results{1}))
    fprintf('[WARN] 未找到完整结果数据，跳过绘图步骤。\n');
else
    h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.45]);
    ax = axes('Position', [0.18, 0.18, 0.75, 0.75]); hold on;
    colors = lines(num_k);

    % 绘制所有正向扫描 (Forward)
    for i = 1:num_k
        if isempty(Results{i}), continue; end
        res = Results{i};
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $K=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:10:length(res.sigma_range));
    end

    % 绘制所有反向扫描 (Backward)
    for i = 1:num_k
        if isempty(Results{i}), continue; end
        res = Results{i};
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $K=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:10:length(res.sigma_range));
    end

    xlabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex'); ylabel('$A(\sigma)$', 'FontSize', 18, 'Interpreter', 'latex');
    grid on; set(ax, 'Box', 'on', 'FontSize', 18, 'TickLabelInterpreter', 'latex');
    legend('Location', 'NorthWest', 'FontSize', 18, 'Interpreter', 'latex', 'NumColumns', 2);
    % title(sprintf('ER Network Hysteresis Loop vs Layers (N=%d, p=%.3f)', DYNA.N, TOPO_PARAM));

    % 图像导出
    plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'figures');
    if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
    fname = fullfile(plots_dir, 'fig5a.eps');
    exportgraphics(h, fname);
    fprintf('[DONE] 绘图已更新: %s\n', fname);
end
