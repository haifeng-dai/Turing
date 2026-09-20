%% plot_sigma_beta_alpha_selected_alpha.m: 绘制指定alpha值下的beta/alpha扫描曲线
clear; clc; close all;
addpath('simulations', 'networks');

% ==========================================
% 0. 参数配置
% ==========================================
TOPO_TYPE = 'ER';
P_VAL_FIXED = 0.03;
SELECTED_ALPHAS = [0.03, 0.05, 0.10, 0.15];  % 指定要绘制的 alpha 值列表

% ==========================================
% 1. 数据加载
% ==========================================
data_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
data_file = fullfile(data_dir, sprintf('scan_beta_alpha_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), 200, 5, P_VAL_FIXED));

if ~exist(data_file, 'file')
    error('数据文件不存在: %s\n请先运行 sigma_beta_alpha.m 生成扫描数据！', data_file);
end

load(data_file);  % 加载 SF_Matrix, SB_Matrix, ALPHA_LIST, RATIO_LIST
fprintf('[LOAD] 已加载数据: %s\n', data_file);

% ==========================================
% 2. 数据提取与处理
% ==========================================
% 找到最接近每个 SELECTED_ALPHAS 的索引
num_alpha = length(SELECTED_ALPHAS);
alpha_indices = zeros(num_alpha, 1);
actual_alphas = zeros(num_alpha, 1);
sigma_f_all = zeros(num_alpha, length(RATIO_LIST));
sigma_b_all = zeros(num_alpha, length(RATIO_LIST));

for i = 1:num_alpha
    [~, alpha_indices(i)] = min(abs(ALPHA_LIST - SELECTED_ALPHAS(i)));
    actual_alphas(i) = ALPHA_LIST(alpha_indices(i));
    sigma_f_all(i, :) = SF_Matrix(alpha_indices(i), :);
    sigma_b_all(i, :) = SB_Matrix(alpha_indices(i), :);
    fprintf('[DATA] 选择 alpha=%.4f (最接近值: %.4f)\n', SELECTED_ALPHAS(i), actual_alphas(i));
end

fprintf('[DATA] 扫描 %d 条曲线 (ratio: %.1f - %.1f)\n', length(RATIO_LIST), RATIO_LIST(1), RATIO_LIST(end));

% ==========================================
% 3. 绘图设置
% ==========================================
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.15, 0.15, 0.7, 0.7]);
ax = axes('Position', [0.12, 0.12, 0.78, 0.78]); hold on;

% 为不同 alpha 分配颜色
colors = jet(num_alpha);

% 循环绘制多条曲线
for i = 1:num_alpha
    % 绘制正向扫描曲线（实线）
    plot(RATIO_LIST, sigma_f_all(i, :), '-', 'Color', colors(i,:), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('\\alpha=%.4f (\\sigma_f)', actual_alphas(i)));

    % 绘制反向扫描曲线（虚线）
    plot(RATIO_LIST, sigma_b_all(i, :), '--', 'Color', colors(i,:), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('\\alpha=%.4f (\\sigma_b)', actual_alphas(i)));
end

% ==========================================
% 4. 图形属性
% ==========================================
xlabel('Coupling Ratio \\eta = \\beta/\\alpha', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Critical Diffusion Ratio \sigma', 'FontSize', 14, 'FontWeight', 'bold');

% 生成标题（显示 alpha 范围）
alpha_range_str = sprintf('%.4f-%.4f', min(actual_alphas), max(actual_alphas));
title(sprintf('Hysteresis vs Coupling Ratio (ER, p=%.3f, \\alpha\\in[%s])', ...
    P_VAL_FIXED, alpha_range_str), 'FontSize', 15, 'FontWeight', 'bold');

% 坐标轴属性
set(ax, 'FontSize', 13, 'Box', 'on');
% xlim([min(RATIO_LIST(RATIO_LIST>0))*0.8, max(RATIO_LIST)*1.2]);
% ylim([min(min(min(sigma_f_all)), min(min(sigma_b_all)))*0.9, ...
%       max(max(max(sigma_f_all)), max(max(sigma_b_all)))*1.1]);

% 网格与图例
grid on; grid minor;
legend('Location', 'best', 'FontSize', 11, 'Box', 'on', 'LineWidth', 1.5, 'NumColumns', 2);

% ==========================================
% 5. 统计滞后信息
% ==========================================
for i = 1:num_alpha
    hysteresis_width = sigma_f_all(i, :) - sigma_b_all(i, :);
    [max_hyst, max_hyst_idx] = max(hysteresis_width);
    fprintf('[HYST] alpha=%.4f: 最大滞后宽度=%.4f (在 ratio=%.1f 处)\n', ...
        actual_alphas(i), max_hyst, RATIO_LIST(max_hyst_idx));
end

% ==========================================
% 6. 图像导出
% ==========================================
plots_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

% 生成文件名（包含 alpha 值信息）
alpha_str = sprintf('a%.4f_to_%.4f', min(actual_alphas), max(actual_alphas));
out_img = fullfile(plots_dir, sprintf('plot_sigma_beta_alpha_%s_%s_p%.3f.png', ...
    alpha_str, lower(TOPO_TYPE), P_VAL_FIXED));
exportgraphics(h, out_img, 'Resolution', 300);
fprintf('[DONE] 图像已生成: %s\n', out_img);
