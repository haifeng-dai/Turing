%% fig3c_permuted.m: BA 多层网络随机置换节点编号后的滞后效应对比
clear; clc; close all;
addpath('simulations', 'networks');


% ==========================================
% 0. 实验参数配置
% ==========================================
TOPO_TYPE = 'BA';
% M_VALS    = [2, 3, 4, 5, 6, 7, 8, 9, 10]; % BA 网络的优先连接参数 (每个新节点增加的边数)
M_VALS    = [3, 4, 10];
PERMUTATION_SEED = 20260927;
PERMUTATION_TAG = sprintf('nodeperm_s%d', PERMUTATION_SEED);

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
num_m = length(M_VALS);

% ==========================================
% 1. 【核心计算区】(缺失时自动生成种子与扫描)
% ==========================================
results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

FORCE_RERUN = false; % true: 强制重算并覆盖旧数据；false: 优先读取已有缓存

% 检查结果是否缺失或强制重算
missing_result = false(1, num_m);
for i = 1:num_m
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_%s_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, M_VALS(i), DYNA.alpha, DYNA.beta, ...
        DYNA.noise, PERMUTATION_TAG);
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

    fprintf('\n[SIM] 正在尝试启动 %d 组缺失的仿真任务 (BA Network)... \n', sum(missing_result));
    simStart = tic;
    for i = find(missing_result)
        dyna_local = DYNA;
        net_file = fullfile(results_dir, 'topology', 'BA', ...
            sprintf('N%d_m%d.mat', DYNA.N, M_VALS(i)));
        net_data = load(net_file, 'nets');
        if ~isfield(net_data, 'nets') || numel(net_data.nets) < DYNA.K
            error('BA 网络文件中的层数不足：%s', net_file);
        end

        % 各层独立打乱节点编号；L(p,p) 与邻接矩阵重标号等价。
        rng(PERMUTATION_SEED + M_VALS(i), 'twister');
        dyna_local.L_intra = net_data.nets(1:DYNA.K);
        for layer_idx = 1:DYNA.K
            node_permutation = randperm(DYNA.N);
            dyna_local.L_intra{layer_idx} = ...
                dyna_local.L_intra{layer_idx}(node_permutation, node_permutation);
        end
        dyna_local.cache_tag = PERMUTATION_TAG;
        sweep_param(TOPO_TYPE, M_VALS(i), dyna_local);
    end
    fprintf('[SIM] 扫描仿真已完成。总运行时间: %.2f 秒。\n', toc(simStart));
end

% ==========================================
% 2. 【数据加载区】
% ==========================================
Results = cell(num_m, 1);
for i = 1:num_m
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_%s_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, M_VALS(i), DYNA.alpha, DYNA.beta, ...
        DYNA.noise, PERMUTATION_TAG);
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
colors = lines(num_m);

% 绘制所有正向扫描 (Forward)
for i = 1:num_m
    res = Results{i};
    plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $m=%d$', M_VALS(i)));
end

% 绘制所有反向扫描 (Backward)
for i = 1:num_m
    res = Results{i};
    plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $m=%d$', M_VALS(i)));
end

% title('Hysteresis Loop across BA Network Preferential Attachment (m)', 'FontSize', 12);
xlabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex'); ylabel('$A(\sigma)$', 'FontSize', 18, 'Interpreter', 'latex');
xlim([10, 18]);
% ylim([0, 80]);
legend('Location', 'NorthWest', 'FontSize', 18, 'Interpreter', 'latex', 'NumColumns', 2);
grid on; set(ax, 'Box', 'on', 'FontSize', 18, 'TickLabelInterpreter', 'latex');

% 图像导出
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'V2', 'manuscirpt', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
% fname = fullfile(plots_dir, 'fig3c_permuted.eps');
% exportgraphics(h, fname);
test_plots_dir = fullfile(fileparts(mfilename('fullpath')), 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end
test_fname = fullfile(test_plots_dir, 'fig3c_permuted.png');
exportgraphics(h, test_fname, 'Resolution', 300);
fprintf('[DONE] 测试图片已保存: %s\n', test_fname);
