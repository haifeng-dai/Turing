%% plot_sigma_alpha_phase_selected_ratio.m
% 为多个指定 Ratio (beta/alpha) 值分别绘制 (alpha, sigma) 相图: Uniform / Bistable / Pattern
% 参考 plot_sigma_beta_phase_selected_alpha.m 的逻辑，改为以 alpha 为横轴
clear; clc; close all;

%% 1. 配置
TOPO_TYPE = 'ER';
N = 200;
K = 5;
P_VAL_FIXED = 0.030;
SELECTED_RATIOS = [0.1, 1.0, 10.0];  % 要绘制的多个 Ratio (eta) 值
RATIO_TOL = 1e-4;       % 匹配容差

% ========== 为每个 Ratio 配置标注文字和箭头坐标 ==========
% 采用相对坐标 (0-1)，易于微调
% config(i).ratio = ...
% ----------------- 第一个子图 (d): Ratio = 0.1 -----------------
ratio_config(1).ratio = 0.1;
% Forward 临界线标示 (蓝色)
ratio_config(1).forward.txt_x_norm = 0.15;         % "\sigma_c^f" 符号横向位置 (0-1)
ratio_config(1).forward.txt_y_norm = 0.87;         % "\sigma_c^f" 符号纵向位置 (0-1)
ratio_config(1).forward.arrow_start_x_norm = 0.15; % 蓝色箭头起点 (X)
ratio_config(1).forward.arrow_start_y_norm = 0.85; % 蓝色箭头起点 (Y)
ratio_config(1).forward.arrow_end_x_norm = 0.05;   % 蓝色箭头终点 (X)
ratio_config(1).forward.arrow_end_y_norm = 0.74;   % 蓝色箭头终点 (Y)
% Backward 临界线标示 (橙色)
ratio_config(1).backward.txt_x_norm = 0.08;        % "\sigma_c^b" 符号横向位置 (0-1)
ratio_config(1).backward.txt_y_norm = 0.10;        % "\sigma_c^b" 符号纵向位置 (0-1)
ratio_config(1).backward.arrow_start_x_norm = 0.11; % 橙色箭头起点 (X)
ratio_config(1).backward.arrow_start_y_norm = 0.13; % 橙色箭头起点 (Y)
ratio_config(1).backward.arrow_end_x_norm = 0.21;   % 橙色箭头终点 (X)
ratio_config(1).backward.arrow_end_y_norm = 0.23;   % 橙色箭头终点 (Y)
% 区域文字位置 (Uniform, Bistable, Pattern)
ratio_config(1).region_x_norms = [0.7, 0.4, 0.6];   % 三个词各自的 X 相对位置 (0-1)
ratio_config(1).region_y_norms = [0.08, 0.2, 0.72]; % 三个词各自的 Y 相对位置 (0-1)

% ----------------- 第二个子图 (e): Ratio = 1.0 -----------------
ratio_config(2).ratio = 1.0;
% Forward 临界线标示 (蓝色)
ratio_config(2).forward.txt_x_norm = 0.15;         % "\sigma_c^f" 符号横向位置
ratio_config(2).forward.txt_y_norm = 0.87;         % "\sigma_c^f" 符号纵向位置
ratio_config(2).forward.arrow_start_x_norm = 0.15; % 蓝色箭头起点 (X)
ratio_config(2).forward.arrow_start_y_norm = 0.85; % 蓝色箭头起点 (Y)
ratio_config(2).forward.arrow_end_x_norm = 0.05;   % 蓝色箭头终点 (X)
ratio_config(2).forward.arrow_end_y_norm = 0.74;   % 蓝色箭头终点 (Y)
% Backward 临界线标示 (橙色)
ratio_config(2).backward.txt_x_norm = 0.02;        % "\sigma_c^b" 符号横向位置
ratio_config(2).backward.txt_y_norm = 0.10;        % "\sigma_c^b" 符号纵向位置
ratio_config(2).backward.arrow_start_x_norm = 0.07; % 橙色箭头起点 (X)
ratio_config(2).backward.arrow_start_y_norm = 0.15; % 橙色箭头起点 (Y)
ratio_config(2).backward.arrow_end_x_norm = 0.17;   % 橙色箭头终点 (X)
ratio_config(2).backward.arrow_end_y_norm = 0.27;   % 橙色箭头终点 (Y)
% 区域文字位置
ratio_config(2).region_x_norms = [0.7, 0.3, 0.6];   % 三个词各自的 X 相对位置
ratio_config(2).region_y_norms = [0.08, 0.25, 0.72]; % 三个词各自的 Y 相对位置

% ----------------- 第三个子图 (f): Ratio = 10.0 -----------------
ratio_config(3).ratio = 10.0;
% Forward 临界线标示 (蓝色)
ratio_config(3).forward.txt_x_norm = 0.17;         % "\sigma_c^f" 符号横向位置
ratio_config(3).forward.txt_y_norm = 0.49;         % "\sigma_c^f" 符号纵向位置
ratio_config(3).forward.arrow_start_x_norm = 0.2; % 蓝色箭头起点 (X)
ratio_config(3).forward.arrow_start_y_norm = 0.43; % 蓝色箭头起点 (Y)
ratio_config(3).forward.arrow_end_x_norm = 0.20;   % 蓝色箭头终点 (X)
ratio_config(3).forward.arrow_end_y_norm = 0.3;   % 蓝色箭头终点 (Y)
% Backward 临界线标示 (橙色)
ratio_config(3).backward.txt_x_norm = 0.35;        % "\sigma_c^b" 符号横向位置
ratio_config(3).backward.txt_y_norm = 0.1;        % "\sigma_c^b" 符号纵向位置
ratio_config(3).backward.arrow_start_x_norm = 0.35; % 橙色箭头起点 (X)
ratio_config(3).backward.arrow_start_y_norm = 0.12; % 橙色箭头起点 (Y)
ratio_config(3).backward.arrow_end_x_norm = 0.25;   % 橙色箭头终点 (X)
ratio_config(3).backward.arrow_end_y_norm = 0.2;   % 橙色箭头终点 (Y)
% 区域文字位置
ratio_config(3).region_x_norms = [0.8, 0.65, 0.8];   % 三个词各自的 X 相对位置
ratio_config(3).region_y_norms = [0.08, 0.5, 0.9]; % 三个词各自的 Y 相对位置

res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
data_path = fullfile(res_dir, sprintf('scan_alpha_ratio_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), N, K, P_VAL_FIXED));

%% 2. 加载数据
if ~exist(data_path, 'file')
    error('未找到数据文件: %s\n请先运行 sigma_alpha_beta.m 生成数据。', data_path);
end
data = load(data_path);

%% 3. 为每个 Ratio 值循环绘图
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

for r_idx = 1:length(SELECTED_RATIOS)
    target_ratio = SELECTED_RATIOS(r_idx);

    % 寻找最接近的 Ratio 索引
    [r_diff, idx_r] = min(abs(data.RATIO_LIST - target_ratio));
    if r_diff > RATIO_TOL
        warning('未找到精确 Ratio=%.2f，使用最近值 %.2f', target_ratio, data.RATIO_LIST(idx_r));
    end
    actual_ratio = data.RATIO_LIST(idx_r);

    % 获取配置
    config = ratio_config(r_idx);

    x_data = data.ALPHA_LIST(:)';  % 横轴：Alpha
    sf_vec = data.SF_Matrix(:, idx_r)'; % 对应列
    sb_vec = data.SB_Matrix(:, idx_r)';

    % 绘图准备 (保持与 noise_impact_thresholds.m 一致)
    h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.45]);
    ax = axes('Position', [0.18, 0.18, 0.75, 0.75]); hold on;

    max_val = max(sf_vec) * 1.15;
    min_val = 0;
    if r_idx == 3
        min_val = 10;
    end
    fill_max_vec = ones(size(x_data)) * max_val;
    fill_min_vec = ones(size(x_data)) * min_val;

    % 填色区域
    fill([x_data, fliplr(x_data)], [sb_vec, fill_min_vec], ...
        [0.90 0.95 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');  % Uniform
    fill([x_data, fliplr(x_data)], [sf_vec, fliplr(sb_vec)], ...
        [0.90 0.90 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');  % Bistable
    fill([x_data, fliplr(x_data)], [fill_max_vec, fliplr(sf_vec)], ...
        [1.00 0.95 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');  % Pattern

    % 临界边界线
    plot(x_data, sf_vec, '--', 'Color', [0.00 0.447 0.741], 'LineWidth', 2.5, 'DisplayName', 'Forward Critical \sigma_c^f');
    plot(x_data, sb_vec, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 2.5, 'DisplayName', 'Backward Critical \sigma_c^b');

    xlabel('$\alpha$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex');
    box on; grid off;
    xlim([0.01, 0.1]);  % 限制横坐标范围
    ylim([min_val, max_val]);
    set(ax, 'FontSize', 18, 'Layer', 'top', 'TickLabelInterpreter', 'latex');

    % --- 预先计算坐标映射 (供区域文字和箭头使用) ---
    ax_pos = get(ax, 'Position');
    xl = xlim(ax); yl = ylim(ax);
    norm2data_x = @(xn) xl(1) + xn * (xl(2) - xl(1));
    norm2data_y = @(yn) yl(1) + yn * (yl(2) - yl(1));
    data2fig_x = @(xd) ax_pos(1) + (xd - xl(1)) / (xl(2) - xl(1)) * ax_pos(3);
    data2fig_y = @(yd) ax_pos(2) + (yd - yl(1)) / (yl(2) - yl(1)) * ax_pos(4);
    clamp01 = @(v) max(0, min(1, v));

    % 区域标签 - 使用解耦的归一化坐标映射
    text(norm2data_x(config.region_x_norms(1)), min_val + config.region_y_norms(1)*(max_val-min_val), 'Uniform', ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
    text(norm2data_x(config.region_x_norms(2)), min_val + config.region_y_norms(2)*(max_val-min_val), 'Bistable', ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
    text(norm2data_x(config.region_x_norms(3)), min_val + config.region_y_norms(3)*(max_val-min_val), 'Pattern', ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

    % Forward 标注
    af_sx = clamp01(data2fig_x(norm2data_x(config.forward.arrow_start_x_norm)));
    af_sy = clamp01(data2fig_y(norm2data_y(config.forward.arrow_start_y_norm)));
    af_ex = clamp01(data2fig_x(norm2data_x(config.forward.arrow_end_x_norm)));
    af_ey = clamp01(data2fig_y(norm2data_y(config.forward.arrow_end_y_norm)));
    annotation(h, 'arrow', [af_sx, af_ex], [af_sy, af_ey], 'Color', [0.00 0.447 0.741], 'LineWidth', 2, 'HeadStyle', 'vback2', 'HeadLength', 8);
    text(norm2data_x(config.forward.txt_x_norm), norm2data_y(config.forward.txt_y_norm), '$\sigma_c^f$', 'Color', [0.00 0.447 0.741], 'FontSize', 18, 'Interpreter', 'latex');

    % Backward 标注
    ab_sx = clamp01(data2fig_x(norm2data_x(config.backward.arrow_start_x_norm)));
    ab_sy = clamp01(data2fig_y(norm2data_y(config.backward.arrow_start_y_norm)));
    ab_ex = clamp01(data2fig_x(norm2data_x(config.backward.arrow_end_x_norm)));
    ab_ey = clamp01(data2fig_y(norm2data_y(config.backward.arrow_end_y_norm)));
    annotation(h, 'arrow', [ab_sx, ab_ex], [ab_sy, ab_ey], 'Color', [0.85 0.325 0.098], 'LineWidth', 2, 'HeadStyle', 'vback2', 'HeadLength', 8);
    text(norm2data_x(config.backward.txt_x_norm), norm2data_y(config.backward.txt_y_norm), '$\sigma_c^b$', 'Color', [0.85 0.325 0.098], 'FontSize', 18, 'Interpreter', 'latex');

    % 保存图片
    fig_suffix = char('d' + r_idx - 1);
    fname = fullfile(plots_dir, sprintf('fig2%s.eps', fig_suffix));
    exportgraphics(h, fname);
    fprintf('[DONE] Ratio=%.2f 的相图已导出: %s\n', actual_ratio, fname);
end
