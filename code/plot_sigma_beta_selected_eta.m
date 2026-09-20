%% plot_sigma_beta_selected_eta.m
% 读取 scan_beta_noise_*.mat 结果并绘制特定 eta 下的 sigma_beta 图
clear; clc; close all;

%% 1. 参数配置
TOPO_TYPE   = 'er';
N           = 200;
K           = 5;
alpha       = 0.050;
P_VAL_FIXED = 0.030;

% 指定想要绘制的 eta (噪声强度) 值列表
% 注意：代码会自动寻找数据中最接近的值，若偏差过大则跳过
TARGET_ETAS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];

%% 2. 构建数据路径并加载
res_dir   = fullfile(fileparts(mfilename('fullpath')), 'results');
data_file = fullfile(res_dir, sprintf('scan_beta_noise_%s_N%d_K%d_a%.3f_p%.3f.mat', ...
    lower(TOPO_TYPE), N, K, alpha, P_VAL_FIXED));

if ~exist(data_file, 'file')
    error('找不到结果文件: %s\n请先运行 sigma_beta_noise_er.m 生成扫描数据。', data_file);
end

load(data_file); % 预期加载: SF_Matrix, SB_Matrix, ETA_LIST, BETA_LIST, P_VAL_FIXED

%% 3. 绘图：横轴为 Beta，纵轴为 Sigma
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.5, 0.6]);
ax = axes('Position', [0.15, 0.15, 0.75, 0.75]); hold on;

% 颜色方案
colors = lines(length(TARGET_ETAS));

for i = 1:length(TARGET_ETAS)
    target = TARGET_ETAS(i);

    % 在 ETA_LIST 中找到最接近 target 的实际索引
    [val, idx] = min(abs(ETA_LIST - target));

    % 检查是否为有效匹配 (容差 1e-4)
    if val > 1e-4
        warning('数据源中不存在精确匹配的 eta=%.3f, 已跳过。可选范围: [%.2f, %.2f]', ...
            target, min(ETA_LIST), max(ETA_LIST));
        continue;
    end

    actual_eta = ETA_LIST(idx);

    % 提取该 eta 对应的正向和反向扫描结果 (SF_Matrix 维度是 [num_e, num_b])
    sf_curve = SF_Matrix(idx, :);
    sb_curve = SB_Matrix(idx, :);

    % 1) 绘制正向扫描 (实线 + 圆圈)
    plot(BETA_LIST, sf_curve, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), ...
        'DisplayName', sprintf('\\eta = %.1f (\\sigma_f)', actual_eta));

    % 2) 绘制反向扫描 (虚线 + 三角)
    plot(BETA_LIST, sb_curve, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 6, ..., ...
        'DisplayName', sprintf('\\eta = %.1f (\\sigma_b)', actual_eta));
end

% 坐标轴与图件标注
xlabel('Inter-layer Coupling Strength \beta', 'FontSize', 14);
ylabel('Critical Diffusion Ratio \sigma', 'FontSize', 14);
title(sprintf('Hysteresis Thresholds for Selected \\eta (ER, p=%.3f)', P_VAL_FIXED), 'FontSize', 14);
xlim([0, 1]);
legend('Location', 'best', 'FontSize', 12);
grid on; box on;
set(ax, 'FontSize', 14);

% 4. 保存图表
plots_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
out_img = fullfile(plots_dir, sprintf('selected_eta_sigma_beta_%s.png', lower(TOPO_TYPE)));
exportgraphics(h, out_img, 'Resolution', 300);

fprintf('[DONE] 绘图已保存至: %s\n', out_img);
