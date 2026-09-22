%% A_sigma_single.m: 序参数(A) 随扩散比(sigma)的变化曲线 (单曲线基准演示)
clear; clc; close all;
addpath('simulations', 'networks');

RUN_SIMULATION = true;  % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)

% ==========================================
% 基本参数设置 (针对单一拓扑)
% ==========================================
TOPO_TYPE  = 'ER';
TOPO_PARAM = 0.03;

DYNA.N = 200;
DYNA.K = 5;             % 待测试层数
DYNA.alpha = 0.01;      % 层内扩散强度
DYNA.beta  = 0.1 * DYNA.alpha;  % 层间耦合强度
DYNA.noise = 0.01;       % 噪声强度
DYNA.T_END = 500;       % 每个点的平衡时长
DYNA.init_perturb = 0.1;

% --- 扫描范围设置 ---
DYNA.sigma_min  = 20;
DYNA.sigma_max  = 50;
DYNA.sigma_npts = 50;

% --- 层间结构 (默认全极通) ---
adj_inter = ones(DYNA.K) - eye(DYNA.K);
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

% --- 1. 执行扫描仿真 ---
if RUN_SIMULATION
    simStart = tic;
    sweep_param(TOPO_TYPE, TOPO_PARAM, DYNA);
    simTime = toc(simStart);
    fprintf('[DONE] 仿真任务已完成！总耗时: %.2f 秒。\n', simTime);
end

% --- 2. 加载结果并绘图 ---
mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, DYNA.noise);
data_path = fullfile(fileparts(mfilename('fullpath')), 'results', mat_name);
load(data_path); % 加载 cfg, sigma_range, A_fwd, A_bwd

% 结果可视化 (滞后环对比)
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.45]);
ax = axes('Position', [0.18, 0.18, 0.75, 0.75]); hold on;
plot(sigma_range, A_fwd, '-o', 'Color', [0 0.447 0.741], 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0 0.447 0.741], 'DisplayName', 'Forward (从均匀态开始)');
plot(sigma_range, A_bwd, '--^', 'Color', [0.85 0.325 0.098], 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Backward (从斑图态开始)');

xlabel('Diffusion Ratio \sigma', 'FontSize', 14); ylabel('Order Parameter A(\sigma)', 'FontSize', 14);
% ylim([0, 100]);
% title(sprintf('Hysteresis Loop Baseline (p=%.3f, \\eta=%.2f)', TOPO_PARAM, DYNA.noise), 'FontSize', 12);
legend('Location', 'NorthWest', 'FontSize', 14);
grid on; box on; set(ax, 'FontSize', 14);

% 保存结果图
plots_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
out_img = fullfile(plots_dir, sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f.png', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, DYNA.noise));
exportgraphics(gcf, out_img, 'Resolution', 300);
fprintf('[DONE] 绘图已更新: %s\n', out_img);
