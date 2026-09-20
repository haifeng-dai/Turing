%% plot_sigma_beta_phase_selected_alpha.m
% 为多个指定 alpha 值分别绘制 (beta, sigma) 相图: Uniform / Bistable / Pattern
% 每个 alpha 对应一个完全独立的图形，文字和箭头坐标完全解耦
clear; clc; close all;

%% 1. 配置
TOPO_TYPE = 'ER';
N = 200;
K = 5;
P_VAL_FIXED = 0.030;
SELECTED_ALPHAS = [0.03, 0.05, 0.10];  % 要绘制的多个 alpha 值
ALPHA_TOL = 1e-4;      % alpha 匹配容差

% ========== 为每个 alpha 配置标注文字和箭头坐标 ==========
% 重要: 所有坐标使用相对值 (0-1)，其中 0 代表轴的起点，1 代表轴的终点
% 结构: alpha_config(i).forward / alpha_config(i).backward
% 每个包含: txt_x_norm, txt_y_norm (文字位置, 0-1),
%          arrow_start_x_norm, arrow_start_y_norm, arrow_end_x_norm, arrow_end_y_norm (0-1)

% ----------------- 第一个子图 (a): alpha = 0.03 -----------------
alpha_config(1).alpha = 0.03;
% Forward 临界线标示 (蓝色)
alpha_config(1).forward.txt_x_norm = 0.25;         % 符号 "\sigma_c^f" 的横向位置 (0-1)
alpha_config(1).forward.txt_y_norm = 0.65;         % 符号 "\sigma_c^f" 的纵向位置 (0-1)
alpha_config(1).forward.arrow_start_x_norm = 0.26; % 蓝色箭头起点 (X)
alpha_config(1).forward.arrow_start_y_norm = 0.6;  % 蓝色箭头起点 (Y)
alpha_config(1).forward.arrow_end_x_norm = 0.23;   % 蓝色箭头指向 (X)
alpha_config(1).forward.arrow_end_y_norm = 0.5;   % 蓝色箭头指向 (Y)
% Backward 临界线标示 (橙色)
alpha_config(1).backward.txt_x_norm = 0.06;        % 符号 "\sigma_c^b" 的横向位置
alpha_config(1).backward.txt_y_norm = 0.08;        % 符号 "\sigma_c^b" 的纵向位置
alpha_config(1).backward.arrow_start_x_norm = 0.1;% 橙色箭头起点 (X)
alpha_config(1).backward.arrow_start_y_norm = 0.13;% 橙色箭头起点 (Y)
alpha_config(1).backward.arrow_end_x_norm = 0.15;  % 橙色箭头指向 (X)
alpha_config(1).backward.arrow_end_y_norm = 0.23;  % 橙色箭头指向 (Y)
% 区域文字位置 (Uniform, Bistable, Pattern)
alpha_config(1).region_x_norms = [0.7, 0.35, 0.7]; % 三个词各自的 X 相对位置
alpha_config(1).region_y_norms = [0.10, 0.3, 0.72]; % 三个词各自的 Y 相对位置

% ----------------- 第二个子图 (b): alpha = 0.05 -----------------
alpha_config(2).alpha = 0.05;
% Forward 临界线标示 (蓝色)
alpha_config(2).forward.txt_x_norm = 0.20;         % 符号 "\sigma_c^f" 的横向位置 (0-1)
alpha_config(2).forward.txt_y_norm = 0.62;         % 符号 "\sigma_c^f" 的纵向位置 (0-1)
alpha_config(2).forward.arrow_start_x_norm = 0.2; % 蓝色箭头起点 (X)
alpha_config(2).forward.arrow_start_y_norm = 0.6; % 蓝色箭头起点 (Y)
alpha_config(2).forward.arrow_end_x_norm = 0.1;   % 蓝色箭头指向 (X)
alpha_config(2).forward.arrow_end_y_norm = 0.48;   % 蓝色箭头指向 (Y)
% Backward 临界线标示 (橙色)
alpha_config(2).backward.txt_x_norm = 0.05;        % 符号 "\sigma_c^b" 的横向位置 (0-1)
alpha_config(2).backward.txt_y_norm = 0.33;        % 符号 "\sigma_c^b" 的纵向位置 (0-1)
alpha_config(2).backward.arrow_start_x_norm = 0.08;% 橙色箭头起点 (X)
alpha_config(2).backward.arrow_start_y_norm = 0.27;% 橙色箭头起点 (Y)
alpha_config(2).backward.arrow_end_x_norm = 0.15;  % 橙色箭头指向 (X)
alpha_config(2).backward.arrow_end_y_norm = 0.14;  % 橙色箭头指向 (Y)
% 区域文字位置 (Uniform, Bistable, Pattern)
alpha_config(2).region_x_norms = [0.7, 0.35, 0.6]; % 三个词各自的 X 相对位置
alpha_config(2).region_y_norms = [0.10, 0.2, 0.70]; % 三个词各自的 Y 相对位置

% ----------------- 第三个子图 (c): alpha = 0.10 -----------------
alpha_config(3).alpha = 0.10;
% Forward 临界线标示 (蓝色)
alpha_config(3).forward.txt_x_norm = 0.05;         % 符号 "\sigma_c^f" 的横向位置 (0-1)
alpha_config(3).forward.txt_y_norm = 0.65;         % 符号 "\sigma_c^f" 的纵向位置 (0-1)
alpha_config(3).forward.arrow_start_x_norm = 0.08; % 蓝色箭头起点 (X)
alpha_config(3).forward.arrow_start_y_norm = 0.6; % 蓝色箭头起点 (Y)
alpha_config(3).forward.arrow_end_x_norm = 0.12;   % 蓝色箭头指向 (X)
alpha_config(3).forward.arrow_end_y_norm = 0.45;   % 蓝色箭头指向 (Y)
% Backward 临界线标示 (橙色)
alpha_config(3).backward.txt_x_norm = 0.08;        % 符号 "\sigma_c^b" 的横向位置
alpha_config(3).backward.txt_y_norm = 0.05;        % 符号 "\sigma_c^b" 的纵向位置
alpha_config(3).backward.arrow_start_x_norm = 0.13;% 橙色箭头起点 (X)
alpha_config(3).backward.arrow_start_y_norm = 0.08;% 橙色箭头起点 (Y)
alpha_config(3).backward.arrow_end_x_norm = 0.2;  % 橙色箭头指向 (X)
alpha_config(3).backward.arrow_end_y_norm = 0.2;  % 橙色箭头指向 (Y)
% 区域文字位置 (Uniform, Bistable, Pattern)
alpha_config(3).region_x_norms = [0.55, 0.7, 0.55]; % 三个词各自的 X 相对位置
alpha_config(3).region_y_norms = [0.1, 0.5, 0.92]; % 三个词各自的 Y 相对位置

res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
data_path = fullfile(res_dir, sprintf('scan_beta_alpha_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), N, K, P_VAL_FIXED));

%% 2. 加载数据
if ~exist(data_path, 'file')
    error('未找到数据文件: %s\n请先运行 sigma_beta_alpha.m 生成数据。', data_path);
end

data = load(data_path);

required_fields = {'SF_Matrix', 'SB_Matrix', 'ALPHA_LIST'};
for i = 1:numel(required_fields)
    if ~isfield(data, required_fields{i})
        error('数据缺少字段: %s', required_fields{i});
    end
end

% 兼容两种横轴字段: RATIO_LIST (beta/alpha) 或 BETA_LIST (beta)
if isfield(data, 'RATIO_LIST')
    is_ratio_axis = true;
else
    is_ratio_axis = false;
end

%% 3. 为每个 alpha 值循环绘图
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

for alpha_idx = 1:length(SELECTED_ALPHAS)
    target_alpha = SELECTED_ALPHAS(alpha_idx);

    % 从配置中获取该 alpha 对应的标注参数
    config = alpha_config(alpha_idx);

    % 确定横轴内容
    if is_ratio_axis
        ratio_vec = data.RATIO_LIST(:)';
        x_data = ratio_vec;  % 使用比值作为横轴
        x_label = '$\beta/\alpha$';
    else
        x_data = data.BETA_LIST(:)';
        x_label = '$\beta$';
    end

    % 寻找与目标 alpha 最接近的一行
    [alpha_diff, idx_a] = min(abs(data.ALPHA_LIST - target_alpha));
    if alpha_diff > ALPHA_TOL
        warning(['未找到精确 alpha=%.4f (容差 %.1e)，改用最近值 alpha=%.4f。\n' ...
            '若需精确值，请调整 sigma_beta_alpha.m 的 ALPHA_LIST。'], ...
            target_alpha, ALPHA_TOL, data.ALPHA_LIST(idx_a));
    end
    actual_alpha = data.ALPHA_LIST(idx_a);

    sf_vec = data.SF_Matrix(idx_a, :);
    sb_vec = data.SB_Matrix(idx_a, :);

    % 避免填充时出现负宽度: 确保 sf >= sb
    swap_mask = sf_vec < sb_vec;
    if any(swap_mask)
        tmp = sf_vec(swap_mask);
        sf_vec(swap_mask) = sb_vec(swap_mask);
        sb_vec(swap_mask) = tmp;
    end

    %% 为本次 alpha 创建独立图形
    h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.46, 0.52]);
    ax = axes('Position', [0.14, 0.20, 0.80, 0.74]); hold on;

    max_val = max(sf_vec) * 1.12;
    min_val = max(0, min(sb_vec) * 0.88);
    if alpha_idx == 1
        min_val = 8;
    end
    fill_max_vec = ones(size(x_data)) * max_val;
    fill_min_vec = ones(size(x_data)) * min_val;

    % 下层: Uniform
    fill([x_data, fliplr(x_data)], [sb_vec, fill_min_vec], ...
        [0.90 0.95 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');

    % 中层: Bistable
    fill([x_data, fliplr(x_data)], [sf_vec, fliplr(sb_vec)], ...
        [0.90 0.90 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');

    % 上层: Pattern
    fill([x_data, fliplr(x_data)], [fill_max_vec, fliplr(sf_vec)], ...
        [1.00 0.95 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');

    % 临界边界
    plot(x_data, sf_vec, '--', 'Color', [0.00 0.447 0.741], 'LineWidth', 2.5, ...
        'DisplayName', 'Forward Critical \sigma_c^f');
    plot(x_data, sb_vec, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 2.5, ...
        'DisplayName', 'Backward Critical \sigma_c^b');

    % 区域文字
    xlabel(x_label, 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex');
    box on; grid off;
    xlim([0, 10]);  % 设置横轴范围为 0-10
    ylim([min_val, max_val]);
    set(ax, 'FontSize', 18, 'Layer', 'top', 'TickLabelInterpreter', 'latex');

    % --- 预先计算坐标映射 (供区域文字和箭头使用) ---
    ax_pos = get(ax, 'Position');
    xl = xlim(ax);
    yl = ylim(ax);
    norm2data_x = @(xn) xl(1) + xn * (xl(2) - xl(1));
    norm2data_y = @(yn) yl(1) + yn * (yl(2) - yl(1));
    data2fig_x = @(xd) ax_pos(1) + (xd - xl(1)) / (xl(2) - xl(1)) * ax_pos(3);
    data2fig_y = @(yd) ax_pos(2) + (yd - yl(1)) / (yl(2) - yl(1)) * ax_pos(4);
    clamp01 = @(v) max(0, min(1, v));

    % 区域文字 - 使用解耦的归一化坐标映射
    t1 = text(norm2data_x(config.region_x_norms(1)), min_val + config.region_y_norms(1)*(max_val-min_val), 'Uniform', ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
    t2 = text(norm2data_x(config.region_x_norms(2)), min_val + config.region_y_norms(2)*(max_val-min_val), 'Bistable', ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
    t3 = text(norm2data_x(config.region_x_norms(3)), min_val + config.region_y_norms(3)*(max_val-min_val), 'Pattern', ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

    % ========== Forward (蓝色) - 完全独立配置 ==========
    % 文字的数据坐标
    txt_f_data_x = norm2data_x(config.forward.txt_x_norm);
    txt_f_data_y = norm2data_y(config.forward.txt_y_norm);

    % 箭头的 figure 归一化坐标
    arrow_f_start_fig_x = clamp01(data2fig_x(norm2data_x(config.forward.arrow_start_x_norm)));
    arrow_f_start_fig_y = clamp01(data2fig_y(norm2data_y(config.forward.arrow_start_y_norm)));
    arrow_f_end_fig_x = clamp01(data2fig_x(norm2data_x(config.forward.arrow_end_x_norm)));
    arrow_f_end_fig_y = clamp01(data2fig_y(norm2data_y(config.forward.arrow_end_y_norm)));

    annotation(h, 'arrow', [arrow_f_start_fig_x, arrow_f_end_fig_x], [arrow_f_start_fig_y, arrow_f_end_fig_y], ...
        'Color', [0.00 0.447 0.741], 'LineWidth', 2, 'HeadStyle', 'vback2', 'HeadLength', 8);
    tf = text(txt_f_data_x, txt_f_data_y, '$\sigma_c^f$', ...
        'HorizontalAlignment', 'left', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.00 0.447 0.741]);

    % ========== Backward (橙色) - 完全独立配置 ==========
    % 文字的数据坐标
    txt_b_data_x = norm2data_x(config.backward.txt_x_norm);
    txt_b_data_y = norm2data_y(config.backward.txt_y_norm);

    % 箭头的 figure 归一化坐标
    arrow_b_start_fig_x = clamp01(data2fig_x(norm2data_x(config.backward.arrow_start_x_norm)));
    arrow_b_start_fig_y = clamp01(data2fig_y(norm2data_y(config.backward.arrow_start_y_norm)));
    arrow_b_end_fig_x = clamp01(data2fig_x(norm2data_x(config.backward.arrow_end_x_norm)));
    arrow_b_end_fig_y = clamp01(data2fig_y(norm2data_y(config.backward.arrow_end_y_norm)));

    annotation(h, 'arrow', [arrow_b_start_fig_x, arrow_b_end_fig_x], [arrow_b_start_fig_y, arrow_b_end_fig_y], ...
        'Color', [0.85 0.325 0.098], 'LineWidth', 2, 'HeadStyle', 'vback2', 'HeadLength', 8);
    tb = text(txt_b_data_x, txt_b_data_y, '$\sigma_c^b$', ...
        'HorizontalAlignment', 'left', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.85 0.325 0.098]);

    uistack([t1, t2, t3, tf, tb], 'top');

    % 导出图片
    fig_suffix = char('a' + alpha_idx - 1);
    out_img = fullfile(plots_dir, sprintf('fig2%s.eps', fig_suffix));
    exportgraphics(h, out_img);
    fprintf('[DONE] 第 %d 个相图已导出 (alpha=%.4f): %s\n', alpha_idx, actual_alpha, out_img);

end
fprintf('\n所有 %d 个相图绘制完成！\n', length(SELECTED_ALPHAS));
