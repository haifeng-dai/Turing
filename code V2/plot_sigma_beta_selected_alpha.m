%% plot_sigma_beta_selected_alpha.m
% 读取 scan_beta_alpha_*.mat 结果并绘制特定 alpha 下的 sigma_beta 图
clear; clc; close all;

%% 1. 配置
% 指定要读取的数据文件 (根据 sigma_beta_alpha.m 生成的文件名)
TOPO_TYPE = 'er';
N = 200;
K = 5;
P_VAL_FIXED = 0.030;
data_file = sprintf('results/scan_beta_alpha_%s_N%d_K%d_p%.3f.mat', TOPO_TYPE, N, K, P_VAL_FIXED);

% 指定要绘制的 alpha 值 (如果数据中不存在精确匹配，则寻找最近的值)
TARGET_ALPHAS = [0.05, 0.10, 0.15];

%% 2. 加载数据
if ~exist(data_file, 'file')
    error('找不到结果文件: %s\n请先运行 sigma_beta_alpha.m 生成数据。', data_file);
end

load(data_file); % 预期变量: SF_Matrix, SB_Matrix, ALPHA_LIST, BETA_LIST

%% 3. 绘图
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.5, 0.6]);
ax = axes; hold on;

% 颜色映射
colors = lines(length(TARGET_ALPHAS));

legend_entries = {};
for i = 1:length(TARGET_ALPHAS)
    target = TARGET_ALPHAS(i);
    [val, idx] = min(abs(ALPHA_LIST - target));

    % 检查是否为有效匹配 (容差设为 1e-4)
    if val > 1e-4
        warning('数据源中不存在精确匹配的 alpha=%.3f, 跳过该绘图。可选范围: [%.3f, %.3f]', ...
            target, min(ALPHA_LIST), max(ALPHA_LIST));
        continue;
    end

    actual_alpha = ALPHA_LIST(idx);

    % 绘制正向扫描阈值
    p_fwd = plot(BETA_LIST, SF_Matrix(idx, :), '-o', 'Color', colors(i,:), ...
        'LineWidth', 2, 'MarkerSize', 6);
    % 绘制反向扫描阈值
    p_bwd = plot(BETA_LIST, SB_Matrix(idx, :), '--s', 'Color', colors(i,:), ...
        'LineWidth', 1.5, 'MarkerSize', 6);

    % 只为正向扫描添加图例项
    legend_entries = [legend_entries, {sprintf('\\alpha = %.3f (Fwd)', actual_alpha), ...
        sprintf('\\alpha = %.3f (Bwd)', actual_alpha)}];
end

% 配置坐标轴
xlabel('Inter-layer Coupling Strength \beta', 'FontSize', 14);
ylabel('Critical Diffusion Ratio \sigma', 'FontSize', 14);
title(sprintf('Hysteresis Thresholds for Selected \\alpha (ER, p=%.3f)', P_VAL_FIXED), 'FontSize', 14);
legend(legend_entries, 'Location', 'best', 'FontSize', 12);
grid on;
set(ax, 'FontSize', 14);

% 保存图片
plots_dir = 'plots';
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
out_img = fullfile(plots_dir, sprintf('selected_alpha_sigma_beta_%s.png', TOPO_TYPE));
exportgraphics(h, out_img, 'Resolution', 300);
fprintf('[DONE] 绘图已保存至: %s\n', out_img);
