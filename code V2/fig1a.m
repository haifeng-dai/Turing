%% A_sigma_noise.m: 序参数(A) 随扩散比(sigma)的变化曲线，对比不同噪声强度(eta)下的滞后效应
clear; clc; close all;
addpath('simulations', 'networks');

% ==========================================
% 1. 基本参数设置 (针对噪声影响 分析)
% ==========================================
TOPO_TYPE  = 'ER';
TOPO_PARAM = 0.03;

DYNA.N = 200;           % 节点数
DYNA.K = 5;             % 层数
DYNA.alpha = 0.05;      % 层内扩散强度
DYNA.beta  = 0.1 * DYNA.alpha;  % 层间耦合强度
DYNA.T_END = 500;       % 扫描点仿真时长
DYNA.steps = 2;         % 输出记录点数 (对 A 统计影响较小)
DYNA.init_perturb = 0.1;

% --- 待对比的噪声强度列表 (eta) ---
ETA_LIST = [0, 0.9, 1];

% --- 扫描范围设置 (Sigma Range) ---
DYNA.sigma_min  = 10;
DYNA.sigma_max  = 24;
DYNA.sigma_npts = 71;   % 对应 0.1 的步长

% 基于全连通配置层间拉普拉斯矩阵
adj_inter = ones(DYNA.K) - eye(DYNA.K);
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;
results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');

% ==========================================
% 2. 【核心计算区】(批量运行滞后扫描)
% ==========================================
% 只在扫描结果缺失时生成标准种子并运行对应噪声的扫描
missing_result = false(size(ETA_LIST));
for i = 1:length(ETA_LIST)
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, ETA_LIST(i));
    missing_result(i) = ~exist(fullfile(results_dir, mat_name), 'file');
end

if any(missing_result)
    seed_name = sprintf('evolution_er_N%d_K%d_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', ...
        DYNA.N, DYNA.K);
    seed_file = fullfile(results_dir, seed_name);
    if ~exist(seed_file, 'file')
        dyna_seed = DYNA;
        dyna_seed.sigma = 100;
        dyna_seed.noise = 0;
        pattern_evolution(dyna_seed, TOPO_TYPE, TOPO_PARAM, false);
    end

    fprintf('\n[INFO] 正在运行 %d 组缺失的滞后扫描 (Parallel)...\n', sum(missing_result));
    simStart = tic;
    for i = find(missing_result)
        dyna_local = DYNA;
        dyna_local.noise = ETA_LIST(i);
        sweep_param(TOPO_TYPE, TOPO_PARAM, dyna_local);
    end
    fprintf('[INFO] 缺失扫描已完成。总耗时: %.2f 秒。\n', toc(simStart));
end

% ==========================================
% 3. 【绘图展示区】(合一对比 A vs Sigma)
% ==========================================
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.45]);
ax = axes('Position', [0.18, 0.18, 0.75, 0.75]);
hold on;
colors = lines(length(ETA_LIST));

% 绘制所有正向扫描 (Forward)
for i = 1:length(ETA_LIST)
    cur_eta = ETA_LIST(i);
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, cur_eta);
    data_file = fullfile(results_dir, mat_name);
    if exist(data_file, 'file')
        res = load(data_file);
        plot(res.sigma_range, res.A_fwd , '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), ...
            'MarkerIndices', 1:4:length(res.sigma_range), ...
            'DisplayName', sprintf('Fwd, $\\eta = %.1f$', cur_eta));
    end
end

% 绘制所有反向扫描 (Backward)
for i = 1:length(ETA_LIST)
    cur_eta = ETA_LIST(i);
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, cur_eta);
    data_file = fullfile(results_dir, mat_name);
    if exist(data_file, 'file')
        res = load(data_file);
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, ...
            'MarkerIndices', 3:4:length(res.sigma_range), ...
            'DisplayName', sprintf('Bwd, $\\eta = %.1f$', cur_eta));
    end
end

% 属性精修
xlabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('$A(\sigma)$', 'FontSize', 18, 'Interpreter', 'latex');
xlim([10, 24]);
% ylim([0, 1]);
legend('Location', 'NorthWest', 'FontSize', 16, 'NumColumns', 2, 'Interpreter', 'latex');
grid on;
set(ax, 'FontSize', 18, 'Box', 'on', 'TickLabelInterpreter', 'latex');

% 添加面板标签 (a) 在左上角外部
annotation(h, 'textbox', [0.026, 0.86, 0.08, 0.08], ...
    'String', '(a)', 'FontSize', 18, ...
    'BackgroundColor', 'none', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% 导出图像
% plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'V2', 'manuscirpt', 'figures');
% if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
% out_img = fullfile(plots_dir, 'fig1a.eps');
% exportgraphics(h, out_img);
% fprintf('[DONE] 绘图已更新: %s\n', out_img);

% 测试环境 PNG 导出
plots_dir = fullfile(fileparts(mfilename('fullpath')), 'fig');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
out_img = fullfile(plots_dir, 'fig1a.png');
exportgraphics(h, out_img, 'Resolution', 300);
fprintf('[DONE] 绘图已保存: %s\n', out_img);
