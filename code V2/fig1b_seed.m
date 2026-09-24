%% fig1b_seed.m: 基于多个噪声随机种子的临界阈值均值与方差
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'simulations'), fullfile(script_dir, 'networks'));

% 实验配置 (与 fig1b 保持一致)
TOPO_TYPE = 'ER';
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta = 0.005;
DYNA.noise = 0.01;
DYNA.T_END = 500;
DYNA.steps = 2;
DYNA.init_perturb = 0.1;
DYNA.sigma_min = 0;
DYNA.sigma_max = 25;

P_LIST = 0.01;
ETA_LIST = linspace(0, 3, 301);
SEED_LIST = 1:5; % 可按需增减重复实验种子
if numel(SEED_LIST) < 2 || numel(unique(SEED_LIST)) ~= numel(SEED_LIST)
    error('SEED_LIST 至少需要包含两个互不重复的随机种子。');
end

res_dir = fullfile(script_dir, 'results');
SF_REPEATS = zeros(numel(SEED_LIST), numel(ETA_LIST));
SB_REPEATS = zeros(numel(SEED_LIST), numel(ETA_LIST));

% 调用方逐个检查种子结果缓存；仅缺失时运行该次扫描。
for seed_idx = 1:numel(SEED_LIST)
    cur_seed = SEED_LIST(seed_idx);
    data_name = sprintf('scan_eta_p_%s_N%d_K%d_a%.3f_b%.3f_seed%.0f.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta, cur_seed);
    data_path = fullfile(res_dir, data_name);

    if ~exist(data_path, 'file')
        seed_file = fullfile(res_dir, sprintf( ...
            'evolution_er_N%d_K%d_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', ...
            DYNA.N, DYNA.K));
        if ~exist(seed_file, 'file')
            rng(0, 'twister');
            dyna_seed = DYNA;
            dyna_seed.sigma = 100;
            dyna_seed.noise = 0;
            pattern_evolution(dyna_seed, TOPO_TYPE, 0.03, false);
        end

        dyna_run = DYNA;
        dyna_run.noise_seed = cur_seed;
        sigma_eta_connectivity_er(dyna_run, TOPO_TYPE, P_LIST, ETA_LIST);
    end

    data = load(data_path);
    idx_p = find(abs(data.P_LIST - 0.01) < 1e-6, 1);
    if isempty(idx_p)
        error('种子 %.0f 的数据中未找到 p=0.01。', cur_seed);
    end
    if ~isequal(data.P_LIST, P_LIST)
        error('种子 %.0f 的数据 P_LIST 与当前配置不一致。', cur_seed);
    end
    if ~isequal(data.ETA_LIST, ETA_LIST)
        error('种子 %.0f 的数据 ETA_LIST 与当前配置不一致。', cur_seed);
    end
    SF_REPEATS(seed_idx, :) = data.SF_Matrix(idx_p, :);
    SB_REPEATS(seed_idx, :) = data.SB_Matrix(idx_p, :);
end

DELTA_REPEATS = SF_REPEATS - SB_REPEATS;
DELTA_MEAN = mean(DELTA_REPEATS, 1);
DELTA_VAR = var(DELTA_REPEATS, 0, 1);
DELTA_STD = sqrt(DELTA_VAR);

% find_thresholds.m 的二分搜索精度 epsilon=0.05；
% |Delta sigma| 在此范围内视作数值上无法区分于零。
DELTA_ZERO_TOL_LIST = [0.025, 0.05, 0.10];
P_ZERO_BY_TOL = zeros(numel(DELTA_ZERO_TOL_LIST), numel(ETA_LIST));
for tol_idx = 1:numel(DELTA_ZERO_TOL_LIST)
    P_ZERO_BY_TOL(tol_idx, :) = mean( ...
        abs(DELTA_REPEATS) <= DELTA_ZERO_TOL_LIST(tol_idx), 1);
end
zero_tol_idx = find(abs(DELTA_ZERO_TOL_LIST - 0.05) < eps, 1);
DELTA_ZERO_TOL = DELTA_ZERO_TOL_LIST(zero_tol_idx);
P_ZERO = P_ZERO_BY_TOL(zero_tol_idx, :);

stats_path = fullfile(res_dir, sprintf( ...
    'scan_eta_p_%s_N%d_K%d_a%.3f_b%.3f_seed_stats.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta));
save(stats_path, 'SEED_LIST', 'P_LIST', 'ETA_LIST', ...
    'SF_REPEATS', 'SB_REPEATS', 'DELTA_REPEATS', 'DELTA_MEAN', ...
    'DELTA_VAR', 'DELTA_STD', 'DELTA_ZERO_TOL_LIST', 'P_ZERO_BY_TOL', ...
    'DELTA_ZERO_TOL', 'P_ZERO');

% 双 y 轴：左侧为平均滞回宽度及 ±1 标准差，右侧为近零概率。
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.48, 0.45]);
ax = axes('Position', [0.16, 0.17, 0.75, 0.75]);
hold(ax, 'on');

yyaxis(ax, 'left');
delta_lower = DELTA_MEAN - DELTA_STD;
delta_upper = DELTA_MEAN + DELTA_STD;
fill(ax, [ETA_LIST, fliplr(ETA_LIST)], ...
    [delta_upper, fliplr(delta_lower)], [0 0.447 0.741], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.18, 'HandleVisibility', 'off');
h_mean = plot(ax, ETA_LIST, DELTA_MEAN, '-', ...
    'Color', [0 0.447 0.741], 'LineWidth', 2);
plot(ax, [ETA_LIST(1), ETA_LIST(end)], [0, 0], 'k:', ...
    'LineWidth', 1.2, 'HandleVisibility', 'off');
ylabel(ax, '$\langle\Delta\sigma\rangle$', 'FontSize', 16, 'Interpreter', 'latex');

delta_min = min(delta_lower);
delta_max = max(delta_upper);
delta_pad = max(0.08 * (delta_max - delta_min), DELTA_ZERO_TOL);
ylim(ax, [delta_min - delta_pad, delta_max + delta_pad]);

yyaxis(ax, 'right');
h_zero = plot(ax, ETA_LIST, P_ZERO, '-o', ...
    'Color', [0.85 0.325 0.098], 'LineWidth', 1.8, ...
    'MarkerSize', 4, 'MarkerIndices', 1:15:numel(ETA_LIST));
ylabel(ax, sprintf('$P(|\\Delta\\sigma|\\leq %.2f)$', DELTA_ZERO_TOL), ...
    'FontSize', 16, 'Interpreter', 'latex');
ylim(ax, [0, 1]);
yticks(ax, 0:0.2:1);

xlabel(ax, '$\eta$', 'FontSize', 16, 'Interpreter', 'latex');
xlim(ax, [0, 3]);
grid(ax, 'on');
box(ax, 'on');
set(ax, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
legend(ax, [h_mean, h_zero], ...
    {'Mean $\Delta\sigma$', sprintf('$P(|\\Delta\\sigma|\\leq %.2f)$', DELTA_ZERO_TOL)}, ...
    'Location', 'best', 'Interpreter', 'latex');

annotation(h, 'textbox', [0.04, 0.89, 0.06, 0.07], 'String', '(b)', ...
    'FontSize', 16, 'BackgroundColor', 'none', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

test_plots_dir = fullfile(script_dir, 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end
test_out_img = fullfile(test_plots_dir, 'fig1b_seed.png');
exportgraphics(h, test_out_img, 'Resolution', 300);
fprintf('[DONE] 滞回宽度均值与近零概率图已保存: %s\n统计数据: %s\n', ...
    test_out_img, stats_path);
