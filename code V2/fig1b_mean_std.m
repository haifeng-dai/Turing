%% fig1b_mean_std.m: fig1b 的正向/反向阈值多种子均值与标准差
clear; clc; close all;

% 与 fig1b.m 保持相同的模型参数、网络和扫描范围
TOPO_TYPE = 'ER';
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta = 0.1 * DYNA.alpha;
DYNA.noise = 0.01;
DYNA.T_END = 500;
DYNA.steps = 2;
DYNA.init_perturb = 0.1;
DYNA.sigma_min = 0;
DYNA.sigma_max = 35;

P_VAL = 0.03;
ETA_LIST = linspace(0, 2, 201);
SEED_LIST = 1:5;
FORCE_RERUN = false;
if numel(SEED_LIST) < 2 || numel(unique(SEED_LIST)) ~= numel(SEED_LIST) || ...
        any(SEED_LIST < 0 | SEED_LIST ~= floor(SEED_LIST))
    error('SEED_LIST 至少需要两个互不重复的非负整数种子。');
end

script_dir = fileparts(mfilename('fullpath'));
res_dir = fullfile(script_dir, 'results');
num_eta = numel(ETA_LIST);
SF_REPEATS = zeros(numel(SEED_LIST), num_eta);
SB_REPEATS = zeros(numel(SEED_LIST), num_eta);

% 每个种子对应一条正向和一条反向阈值曲线；已有文件直接复用。
for seed_idx = 1:numel(SEED_LIST)
    seed_value = SEED_LIST(seed_idx);
    data_path = fullfile(res_dir, sprintf( ...
        'scan_eta_p%.3f_%s_N%d_K%d_a%.3f_b%.3f_seed%.0f.mat', ...
        P_VAL, lower(TOPO_TYPE), DYNA.N, DYNA.K, ...
        DYNA.alpha, DYNA.beta, seed_value));
    if FORCE_RERUN || ~exist(data_path, 'file')
        fprintf('[SIM] 计算 p=%.3f、seed=%d 的阈值扫描。\n', P_VAL, seed_value);
        dyna_seed = DYNA;
        dyna_seed.noise_seed = seed_value;
        sigma_eta_connectivity_er(dyna_seed, TOPO_TYPE, P_VAL, ETA_LIST);
    else
        fprintf('[CACHE] 读取 seed=%d：%s\n', seed_value, data_path);
    end

    data = load(data_path, 'p_val', 'ETA_LIST', 'sf_vec', 'sb_vec');
    if ~isfield(data, 'p_val') || abs(data.p_val-P_VAL) > 1e-10 || ...
            ~isfield(data, 'ETA_LIST') || ...
            ~isequal(data.ETA_LIST(:), ETA_LIST(:)) || ...
            ~isfield(data, 'sf_vec') || numel(data.sf_vec) ~= num_eta || ...
            ~isfield(data, 'sb_vec') || numel(data.sb_vec) ~= num_eta || ...
            any(~isfinite(data.sf_vec(:))) || any(~isfinite(data.sb_vec(:)))
        error('seed=%d 的阈值文件与当前扫描参数不符：%s', seed_value, data_path);
    end
    SF_REPEATS(seed_idx, :) = data.sf_vec(:)';
    SB_REPEATS(seed_idx, :) = data.sb_vec(:)';
end

eta_range = ETA_LIST;
sf_vec = mean(SF_REPEATS, 1);
sb_vec = mean(SB_REPEATS, 1);
sf_var = var(SF_REPEATS, 0, 1);
sb_var = var(SB_REPEATS, 0, 1);
sf_std = sqrt(sf_var);
sb_std = sqrt(sb_var);

stats_path = fullfile(res_dir, sprintf( ...
    'scan_eta_p%.3f_%s_N%d_K%d_a%.3f_b%.3f_threshold_stats.mat', ...
    P_VAL, lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta));
save(stats_path, 'P_VAL', 'ETA_LIST', 'SEED_LIST', 'SF_REPEATS', 'SB_REPEATS', ...
    'sf_vec', 'sb_vec', 'sf_var', 'sb_var', 'sf_std', 'sb_std');

% 以下布局、颜色、区域文字和箭头坐标与 fig1b.m 一致。
label_forward.txt_x = 0.47;
label_forward.txt_y = 30;
label_forward.arrow_start_x = 0.47;
label_forward.arrow_start_y = 29.5;
label_forward.arrow_end_x = 0.40;
label_forward.arrow_end_y = 25.7;

label_backward.txt_x = 0.80;
label_backward.txt_y = 8;
label_backward.arrow_start_x = 0.8;
label_backward.arrow_start_y = 9.0;
label_backward.arrow_end_x = 1;
label_backward.arrow_end_y = 12;

h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.43]);
ax = axes('Position', [0.18, 0.18, 0.75, 0.75]);
hold(ax, 'on');

max_val = max(sf_vec) * 1.15;
min_val = min(sb_vec) * 0.85;
fill_max_vec = ones(size(eta_range)) * max_val;
fill_min_vec = ones(size(eta_range)) * min_val;

% 三个区域由两条均值曲线分界，与原图的定义相同。
fill(ax, [eta_range, fliplr(eta_range)], [sb_vec, fill_min_vec], ...
    [0.9 0.95 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');
fill(ax, [eta_range, fliplr(eta_range)], [sf_vec, fliplr(sb_vec)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.8, 'HandleVisibility', 'off');
fill(ax, [eta_range, fliplr(eta_range)], [fill_max_vec, fliplr(sf_vec)], ...
    [1.0 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');

% 均值 ± 1 个样本标准差；标准差与纵轴 sigma 的单位相同。
forward_color = [0 0.447 0.741];
backward_color = [0.85 0.325 0.098];
fill(ax, [eta_range, fliplr(eta_range)], ...
    [sf_vec + sf_std, fliplr(sf_vec - sf_std)], forward_color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.18, 'HandleVisibility', 'off');
fill(ax, [eta_range, fliplr(eta_range)], ...
    [sb_vec + sb_std, fliplr(sb_vec - sb_std)], backward_color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.18, 'HandleVisibility', 'off');

plot(ax, eta_range, sf_vec, '--', 'Color', forward_color, 'LineWidth', 2.5, ...
    'DisplayName', 'Forward Critical \sigma_c^f');
plot(ax, eta_range, sb_vec, '--', 'Color', backward_color, 'LineWidth', 2.5, ...
    'DisplayName', 'Backward Critical \sigma_c^b');

t1 = text(0.3, 12, 'Uniform', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
t2 = text(0.23, 24, 'Bistable', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
t3 = text(1.5, 26, 'Pattern', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

xlabel('$\eta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex');
grid off; box on;
xlim([0, 2]); ylim([min_val, 32]);
set(ax, 'FontSize', 18, 'Layer', 'top', 'TickLabelInterpreter', 'latex');

ax_pos = get(ax, 'Position');
xl = xlim(ax); yl = ylim(ax);
d2f_x = @(xd) ax_pos(1) + (xd - xl(1))/(xl(2)-xl(1)) * ax_pos(3);
d2f_y = @(yd) ax_pos(2) + (yd - yl(1))/(yl(2)-yl(1)) * ax_pos(4);

annotation(h, 'arrow', ...
    [d2f_x(label_forward.arrow_start_x), d2f_x(label_forward.arrow_end_x)], ...
    [d2f_y(label_forward.arrow_start_y), d2f_y(label_forward.arrow_end_y)], ...
    'Color', forward_color, 'LineWidth', 2, 'HeadStyle', 'vback2', 'HeadLength', 8);
tf = text(label_forward.txt_x, label_forward.txt_y, '$\sigma_c^f$', ...
    'HorizontalAlignment', 'left', 'FontSize', 18, ...
    'Interpreter', 'latex', 'Color', forward_color);

annotation(h, 'arrow', ...
    [d2f_x(label_backward.arrow_start_x), d2f_x(label_backward.arrow_end_x)], ...
    [d2f_y(label_backward.arrow_start_y), d2f_y(label_backward.arrow_end_y)], ...
    'Color', backward_color, 'LineWidth', 2, 'HeadStyle', 'vback2', 'HeadLength', 8);
tb = text(label_backward.txt_x, label_backward.txt_y, '$\sigma_c^b$', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, ...
    'Interpreter', 'latex', 'Color', backward_color);

uistack([t1, t2, t3, tf, tb], 'top');
annotation(h, 'textbox', [0.05, 0.86, 0.08, 0.08], ...
    'String', '(b)', 'FontSize', 18, ...
    'BackgroundColor', 'none', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

test_plots_dir = fullfile(script_dir, 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end
test_out_img = fullfile(test_plots_dir, 'fig1b_mean_std.png');
exportgraphics(h, test_out_img, 'Resolution', 300);
fprintf('[DONE] 多种子均值与标准差相图已保存：%s\n', test_out_img);
fprintf('[SAVE] 正向/反向阈值统计数据已保存：%s\n', stats_path);
