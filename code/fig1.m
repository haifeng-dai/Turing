%% fig1abcd.m: 综合绘图脚本 - 在统一图窗中展示fig1a、fig1b、fig1c、fig1d
% 布局: 2x2 网格，其中 (c) 和 (d) 内部包含多个子图
clear; clc; close all;
addpath('simulations', 'networks');

%% ============================================================
%% 初始化参数
%% ============================================================
TOPO_TYPE  = 'ER';
TOPO_PARAM = 0.03;
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta  = 0.1 * DYNA.alpha;

% Figure 1a 参数
ETA_LIST = [0, 0.3, 0.5];
DYNA.sigma_min  = 0;
DYNA.sigma_max  = 24;
DYNA.sigma_npts = 121;
adj_inter = ones(DYNA.K) - eye(DYNA.K);
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

% Figure 1cd 参数
DYNA_CD.sigma = 16.0;
scan_mode = 'bwd';
noise_list = [0.00, 0.30];

%% ============================================================
%% 创建主图窗和布局 (2x2 网格)
%% ============================================================
hf = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.9]);
tl_main = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

%% ============================================================
%% (a) 上左 (1,1): Figure 1a - A vs Sigma 曲线
%% ============================================================
ax_a = nexttile(tl_main, 1);
hold on;
colors = lines(length(ETA_LIST));

% 加载并绘制正向扫描 (Forward)
for i = 1:length(ETA_LIST)
    cur_eta = ETA_LIST(i);
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, cur_eta);
    data_file = fullfile('results', mat_name);
    if exist(data_file, 'file')
        res = load(data_file);
        plot(res.sigma_range, res.A_fwd , '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 5, 'MarkerFaceColor', colors(i,:), ...
            'MarkerIndices', 1:4:length(res.sigma_range), ...
            'DisplayName', sprintf('Fwd, $\\eta = %.1f$', cur_eta));
    end
end

% 加载并绘制反向扫描 (Backward)
for i = 1:length(ETA_LIST)
    cur_eta = ETA_LIST(i);
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, cur_eta);
    data_file = fullfile('results', mat_name);
    if exist(data_file, 'file')
        res = load(data_file);
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 5, ...
            'MarkerIndices', 3:4:length(res.sigma_range), ...
            'DisplayName', sprintf('Bwd, $\\eta = %.1f$', cur_eta));
    end
end

xlabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$A(\sigma)$', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 24]); ylim([0, 148]);
legend('Location', 'NorthWest', 'FontSize', 18, 'NumColumns', 1, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
grid on;
set(ax_a, 'FontSize', 18, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
% 添加面板标签 (a)
annotation(hf, 'textbox', [0.005, 0.93, 0.05, 0.05], 'String', '(a)', ...
    'FontSize', 18, 'FontWeight', 'bold', 'BackgroundColor', 'none', 'EdgeColor', 'none', 'FontName', 'Times New Roman');

%% ============================================================
%% (b) 上右 (1,2): Figure 1b - 三相图
%% ============================================================
ax_b = nexttile(tl_main, 2);
hold on;

% 加载数据
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
data_path = fullfile(res_dir, sprintf('scan_eta_p_%s_N%d_K%d_a%.3f_b%.3f.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta));

if exist(data_path, 'file')
    data = load(data_path);
    idx_p = find(abs(data.P_LIST - TOPO_PARAM) < 1e-6, 1);

    if ~isempty(idx_p)
        eta_range = data.ETA_LIST;
        sf_vec = data.SF_Matrix(idx_p, :);
        sb_vec = data.SB_Matrix(idx_p, :);

        % 动态确定 Y 轴边界
        max_val = max(sf_vec) * 1.15;
        min_val = min(sb_vec) * 0.85;
        fill_max_vec = ones(size(eta_range)) * max_val;
        fill_min_vec = ones(size(eta_range)) * min_val;

        % 区域填充
        fill([eta_range, fliplr(eta_range)], [sb_vec, fill_min_vec], ...
            [0.9 0.95 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off', 'Parent', ax_b);
        fill([eta_range, fliplr(eta_range)], [sf_vec, fliplr(sb_vec)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.8, 'HandleVisibility', 'off', 'Parent', ax_b);
        fill([eta_range, fliplr(eta_range)], [fill_max_vec, fliplr(sf_vec)], ...
            [1.0 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off', 'Parent', ax_b);

        % 边界曲线
        plot(ax_b, eta_range, sf_vec, '--', 'Color', [0 0.447 0.741], 'LineWidth', 2.5, ...
            'DisplayName', 'Forward $\sigma_c^f$');
        plot(ax_b, eta_range, sb_vec, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 2.5, ...
            'DisplayName', 'Backward $\sigma_c^b$');

        % 区域文字标签
        text(ax_b, 0.3, 11, 'Stable', 'HorizontalAlignment', 'center', 'FontSize', 18, ...
            'Interpreter', 'latex', 'Color', [0.2 0.2 0.2], 'FontName', 'Times New Roman');
        text(ax_b, 0.13, 18, 'Hysteresis', 'HorizontalAlignment', 'center', 'FontSize', 18, ...
            'Interpreter', 'latex', 'Color', [0.2 0.2 0.2], 'FontName', 'Times New Roman');
        text(ax_b, 0.9, 15, 'Turing', 'HorizontalAlignment', 'center', 'FontSize', 18, ...
            'Interpreter', 'latex', 'Color', [0.2 0.2 0.2], 'FontName', 'Times New Roman');

        xlabel(ax_b, '$\eta$', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
        ylabel(ax_b, '$\sigma$', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
        xlim(ax_b, [0, 1]); ylim(ax_b, [9, max_val]);
        legend(ax_b, 'Location', 'northeast', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    end
end

grid(ax_b, 'off');
box(ax_b, 'on');
set(ax_b, 'FontSize', 18, 'Layer', 'top', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
% 添加面板标签 (b)
annotation(hf, 'textbox', [0.485, 0.93, 0.05, 0.05], 'String', '(b)', ...
    'FontSize', 18, 'FontWeight', 'bold', 'BackgroundColor', 'none', 'EdgeColor', 'none', 'FontName', 'Times New Roman');

%% ============================================================
%% (c) 下左 (2,1): 轨迹对比 (eta=0.0 vs eta=0.3)
%% ============================================================
U_data = cell(2, 1);
t_data = cell(2, 1);

% 加载轨迹数据
results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
for i = 1:2
    cur_noise = noise_list(i);
    mat_name = sprintf('evolution_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_s%.1f_n%.2f_%s_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, DYNA_CD.sigma, cur_noise, scan_mode);
    file_path = fullfile(results_dir, mat_name);

    if exist(file_path, 'file')
        data = load(file_path, 't', 'Y');
        t_data{i} = data.t;
        U_data{i} = data.Y(:, 1:DYNA.N);
    else
        warning('未找到数据文件: %s', mat_name);
    end
end

% 自动探测特定节点 (突变轨线) - 同 fig1cd.m
idx_switch = [];
if ~isempty(U_data{2})
    u2 = U_data{2};
    u2_end = u2(end, :);
    u2_early = u2(round(size(u2, 1)*0.2), :);
    idx_switch = find(u2_early > 4 & u2_end < 2, 1);
    if isempty(idx_switch)
        [~, idx_switch] = max(max(u2) - min(u2));
    end
end

% 绘制轨迹
c_u = [0.0 0.447 0.741];
% 嵌套布局
tl_c = tiledlayout(tl_main, 2, 1, 'Padding', 'none', 'TileSpacing', 'compact');
tl_c.Layout.Tile = 3;
for i = 1:2
    ax_c_sub = nexttile(tl_c);
    hold(ax_c_sub, 'on');

    if ~isempty(U_data{i})
        t_cur = t_data{i};
        u_cur = U_data{i};
        u_end = u_cur(end, :);
        [~, idx_up] = max(u_end);
        [~, idx_down] = min(u_end);

        plot(ax_c_sub, t_cur, u_cur(:, 1:5:end), 'Color', [0.85 0.85 0.85]);
        plot(ax_c_sub, t_cur, u_cur(:, idx_up), 'Color', c_u, 'LineWidth', 1.5);
        plot(ax_c_sub, t_cur, u_cur(:, idx_down), 'Color', c_u*0.6, 'LineWidth', 1.5);
        if i == 2 && ~isempty(idx_switch)
            plot(ax_c_sub, t_cur, u_cur(:, idx_switch), 'Color', [0.9 0.4 0], 'LineWidth', 2.0);
        end
    end

    ylabel(ax_c_sub, '$u(t)$', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    xlim(ax_c_sub, [0, 500]);
    grid(ax_c_sub, 'on');
    set(ax_c_sub, 'FontSize', 18, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');

    if ~isempty(U_data{i})
        text(ax_c_sub, 0.95, 0.8, sprintf('$\\eta = %.1f$', noise_list(i)), ...
            'Units', 'normalized', 'FontSize', 18, 'Interpreter', 'latex', ...
            'BackgroundColor', 'w', 'EdgeColor', 'w', 'Margin', 2, 'HorizontalAlignment', 'right', 'FontName', 'Times New Roman');
    end

    if i == 2
        xlabel(ax_c_sub, '$t$', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    end
end

% 添加面板标签 (c)
annotation(hf, 'textbox', [0.005, 0.45, 0.05, 0.05], 'String', '(c)', ...
    'FontSize', 18, 'FontWeight', 'bold', 'BackgroundColor', 'none', 'EdgeColor', 'none', 'FontName', 'Times New Roman');

%% ============================================================
%% (d) 第3-4行右侧: 末态节点分布对比
%% ============================================================

% 绘制散点分布 (带有共享侧边色条)
if ~isempty(U_data{1}) && ~isempty(U_data{2})
    u1_end = U_data{1}(end, :);
    u2_end = U_data{2}(end, :);
    c_min = min([min(u1_end), min(u2_end)]);
    c_max = max([max(u1_end), max(u2_end)]);

    % 嵌套布局: 32列高分辨率网格，散点图占29列，留出3/32右边距
    tl_d = tiledlayout(tl_main, 2, 32, 'Padding', 'none', 'TileSpacing', 'compact');
    tl_d.Layout.Tile = 4;
    ax_d_subs = gobjects(2, 1);
    for i = 1:2
        u_end = U_data{i}(end, :);
        % 每行起始格子基于32列计算，占据前29列
        tile_idx = (i-1)*32 + 1;
        ax_d_subs(i) = nexttile(tl_d, tile_idx, [1 29]);
        hold(ax_d_subs(i), 'on');

        scatter(ax_d_subs(i), 1:DYNA.N, u_end, 20, u_end, 'filled');
        colormap(ax_d_subs(i), 'jet');

        ylabel(ax_d_subs(i), '$u(x)$', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
        set(ax_d_subs(i), 'FontSize', 18, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
        clim(ax_d_subs(i), [c_min c_max]);
        ylim(ax_d_subs(i), [min(u_end), max(u_end)]);
        grid(ax_d_subs(i), 'on');

        text(ax_d_subs(i), 0.95, 0.3, sprintf('$\\eta = %.1f$', noise_list(i)), ...
            'Units', 'normalized', 'FontSize', 18, 'Interpreter', 'latex', ...
            'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 2, 'HorizontalAlignment', 'right', 'FontName', 'Times New Roman');

        if i == 2
            xlabel(ax_d_subs(i), 'Node ID', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
        end
    end

    % --- 终极手动坐标控制方式 ---
    drawnow;

    % 在这里下达色条指令
    cb = colorbar(ax_d_subs(2));
    set(cb, 'FontSize', 18, 'FontName', 'Times New Roman');
    cb.Label.FontSize = 18;

    % 手动设置全局标准化坐标 [横轴x, 纵轴y, 宽度w, 高度h]
    % 调整 x 的值（第一个参数），让它对齐 (b) 图右边缘
    cb.Position = [0.923, 0.08, 0.012, 0.39];
end

% 添加面板标签 (d)
annotation(hf, 'textbox', [0.485, 0.45, 0.05, 0.05], 'String', '(d)', ...
    'FontSize', 18, 'FontWeight', 'bold', 'BackgroundColor', 'none', 'EdgeColor', 'none', 'FontName', 'Times New Roman');

% %% ============================================================
% %% 导出最终图像
% %% ============================================================
out_fig = fullfile(plots_dir, 'fig1.eps');
exportgraphics(hf, out_fig, 'ContentType', 'vector');
fprintf('[SUCCESS] 综合图表已导出: %s\n', out_fig);
