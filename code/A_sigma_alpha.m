%% A_sigma_alpha.m: 序参数(A) 随扩散比(sigma)的变化曲线，对比不同层内扩散强度(alpha)下的滞后效应
clear; clc; close all;
addpath('simulations', 'networks');

RUN_SIMULATION = true;  % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)

% ==========================================
% 1. 基本参数设置 (针对层内扩散分析)
% ==========================================
TOPO_TYPE  = 'ER';
TOPO_PARAM = 0.03;

DYNA.N = 200;           % 节点数
DYNA.K = 5;             % 层数
DYNA.noise = 0.01;      % 固定噪声强度
DYNA.T_END = 500;       % 扫描点仿真时长
DYNA.steps = 2;         % 输出记录点数 (对 A 统计影响较小)
DYNA.init_perturb = 0.1;

% --- 待对比的层内扩散强度列表 (alpha) ---
ALPHA_LIST = [0.01, 0.05, 0.10];
% beta 固定为 alpha 的 0.1 倍
BETA_LIST = 0.1 * ALPHA_LIST;

% --- 扫描范围设置 (Sigma Range) ---
DYNA.sigma_min  = 10;
DYNA.sigma_max  = 60;
DYNA.sigma_npts = 201;   % 对应 0.2 的步长

% 基于全连通配置层间拉普拉斯矩阵
adj_inter = ones(DYNA.K) - eye(DYNA.K);
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

% ==========================================
% 2. 【核心计算区】(批量运行滞后扫描)
% ==========================================
if RUN_SIMULATION
    fprintf('\n[INFO] 正在启动 %d 批量滞后扫描仿真 (Parallel)...\n', length(ALPHA_LIST));
    simStart = tic;
    for i = 1:length(ALPHA_LIST)
        dyna_local = DYNA;
        dyna_local.alpha = ALPHA_LIST(i);
        dyna_local.beta = BETA_LIST(i);  % beta = 0.1 * alpha

        % 调用更新后的仿真接口
        sweep_param(TOPO_TYPE, TOPO_PARAM, dyna_local);
    end
    fprintf('[INFO] 批量扫描完成。总耗时: %.2f 秒。\n', toc(simStart));
end

% ==========================================
% 3. 【绘图展示区】(合一对比 A vs Sigma)
% ==========================================
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.45]);
ax = axes('Position', [0.18, 0.18, 0.75, 0.75]);
hold on;
colors = lines(length(ALPHA_LIST));

for i = 1:length(ALPHA_LIST)
    cur_alpha = ALPHA_LIST(i);
    cur_beta = BETA_LIST(i);

    % 根据约定规则动态生成文件名并加载结果
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, cur_alpha, cur_beta, DYNA.noise);
    data_file = fullfile('results', mat_name);

    if exist(data_file, 'file')
        res = load(data_file);

        % 绘制正向扫描 (Forward)
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), ...
            'MarkerIndices', 1:4:length(res.sigma_range), ...
            'DisplayName', sprintf('Forward (\\alpha = %.3f)', cur_alpha));
        % 绘制反向扫描 (Backward)
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, ...
            'MarkerIndices', 3:4:length(res.sigma_range), ...
            'DisplayName', sprintf('Backward (\\alpha = %.3f)', cur_alpha));
    else
        warning('未找到结果数据: %s，请确认计算区仿真已执行。', mat_name);
    end
end

% 属性精修
xlabel('\sigma (Diffusion Ratio)', 'FontSize', 14);
ylabel('Order Parameter A(\sigma)', 'FontSize', 14);
xlim([DYNA.sigma_min, DYNA.sigma_max]); ylim([0, 150]);
legend('Location', 'NorthWest', 'FontSize', 10, 'NumColumns', 2);
xlim([10, 40]);
grid off; box on;
set(ax, 'FontSize', 14, 'Box', 'on');

% 导出图像
plots_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
out_img = fullfile(plots_dir, sprintf('hysteresis_comp_alpha_%s_N%d_K%d_p%.3f_b_0.1a_n%.0f.png', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.noise));
exportgraphics(h, out_img, 'Resolution', 300);
fprintf('[DONE] 绘图已生成: %s\n', out_img);
