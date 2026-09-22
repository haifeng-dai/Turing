%% A_sigma_connectivity_ws.m: 不同重连概率(p)下的WS网络滞后效应对比
clear; clc; close all;
addpath('simulations', 'networks');

RUN_SIMULATION = false;  % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)

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
% 1. 【核心计算区】(由 RUN_SIMULATION 控制)
% ==========================================
if RUN_SIMULATION
    fprintf('\n[SIM] 正在尝试启动 %d 并行仿真任务 (WS Network)... \n', num_k);
    simStart = tic;
    for i = 1:num_k
        % WS 网络参数：K 为平均度，P_REWIRE 为固定的重连概率
        % 创建临时配置以传递两个参数
        dyna_temp = DYNA;
        dyna_temp.ws_k = K_VALS(i);          % 平均度
        dyna_temp.ws_p = P_REWIRE;           % 重连概率
        sweep_param_ws(TOPO_TYPE, K_VALS(i), P_REWIRE, dyna_temp);
    end
    fprintf('[SIM] 扫描仿真已完成。总运行时间: %.2f 秒。\n', toc(simStart));
end

% ==========================================
% 2. 【数据加载区】(保证绘图区始终有数据源)
% ==========================================
Results = cell(num_k, 1);
for i = 1:num_k
    % WS网络文件名包含平均度K和重连概率p_rewire
    mat_name = sprintf('hysteresis_%s_N%d_K%d_pr%.2f_mK%d_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, P_REWIRE, K_VALS(i), DYNA.alpha, DYNA.beta, DYNA.noise);
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
xlim([10, 30]); ylim([0, 140]);
legend('Location', 'NorthWest', 'FontSize', 18, 'Interpreter', 'latex', 'NumColumns', 2);
grid on; set(ax, 'Box', 'on', 'FontSize', 18, 'TickLabelInterpreter', 'latex');

% 图像导出
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'V2', 'manuscirpt', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
fname = fullfile(plots_dir, 'fig3b.eps');
exportgraphics(h, fname);
fprintf('[DONE] 绘图已更新: %s\n', fname);
