%% fig1b_seed.m: 基于多个噪声随机种子的临界阈值均值与方差 (多 p 合一绘图)
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'simulations'), fullfile(script_dir, 'networks'));

% 实验配置 (与 fig1b 保持一致)
TOPO_TYPE = 'ER';
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta = 0.005;
DYNA.T_END = 500;
DYNA.steps = 2;
DYNA.init_perturb = 0.1;
DYNA.sigma_min = 0;
DYNA.sigma_max = 35;

% 扫描配置：外层遍历 P_LIST，内层遍历 SEED_LIST
P_LIST = [0.02, 0.04, 0.06]; % 配置需要对比的多个 p
ETA_LIST = linspace(0, 2, 201);
SEED_LIST = 1:5;             % 重复实验种子列表
if numel(SEED_LIST) < 2 || numel(unique(SEED_LIST)) ~= numel(SEED_LIST)
    error('SEED_LIST 至少需要包含两个互不重复的随机种子。');
end

FORCE_RERUN = true;         % 是否强制重算: true = 强制重新计算, false = 优先读取已有缓存
res_dir = fullfile(script_dir, 'results');
test_plots_dir = fullfile(script_dir, 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end

num_p = length(P_LIST);
num_e = length(ETA_LIST);

% 预分配跨 p 的统计结果存储
ALL_DELTA_MEAN = zeros(num_p, num_e);
ALL_DELTA_STD  = zeros(num_p, num_e);
ALL_P_ZERO     = zeros(num_p, num_e);
DELTA_ZERO_TOL = 0.05;

%% ============================================================
%% 1. 数据计算与加载区 (外层遍历每个 p，独立持久化)
%% ============================================================
for p_idx = 1:num_p
    P_VAL = P_LIST(p_idx);
    fprintf('\n============================================================\n');
    fprintf('[DATA] 正在处理 p = %.3f (%d/%d)...\n', P_VAL, p_idx, num_p);
    fprintf('============================================================\n');

    SF_REPEATS = zeros(numel(SEED_LIST), num_e);
    SB_REPEATS = zeros(numel(SEED_LIST), num_e);

    for seed_idx = 1:numel(SEED_LIST)
        cur_seed = SEED_LIST(seed_idx);
        data_name = sprintf('scan_eta_p%.3f_%s_N%d_K%d_a%.3f_b%.3f_seed%.0f.mat', ...
            P_VAL, lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta, cur_seed);
        data_path = fullfile(res_dir, data_name);

        if FORCE_RERUN || ~exist(data_path, 'file')
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
            sigma_eta_connectivity_er(dyna_run, TOPO_TYPE, P_VAL, ETA_LIST);
        end

        data = load(data_path);
        if abs(data.p_val - P_VAL) > 1e-6
            error('种子 %.0f 的数据 p_val=%.3f 与当前配置 P_VAL=%.3f 不一致。', cur_seed, data.p_val, P_VAL);
        end
        if ~isequal(data.ETA_LIST, ETA_LIST)
            error('种子 %.0f 的数据 ETA_LIST 与当前配置不一致。', cur_seed);
        end
        SF_REPEATS(seed_idx, :) = data.sf_vec;
        SB_REPEATS(seed_idx, :) = data.sb_vec;
    end

    % 统计当前 p 的滞回宽度均值、方差
    DELTA_REPEATS = SF_REPEATS - SB_REPEATS;
    DELTA_MEAN = mean(DELTA_REPEATS, 1);
    DELTA_STD  = sqrt(var(DELTA_REPEATS, 0, 1));
    P_ZERO     = mean(abs(DELTA_REPEATS) <= DELTA_ZERO_TOL, 1);

    % 暂存到全局矩阵供统一绘图
    ALL_DELTA_MEAN(p_idx, :) = DELTA_MEAN;
    ALL_DELTA_STD(p_idx, :)  = DELTA_STD;
    ALL_P_ZERO(p_idx, :)     = P_ZERO;

    % 保存单个 p 的统计结果
    stats_path = fullfile(res_dir, sprintf( ...
        'scan_eta_p%.3f_%s_N%d_K%d_a%.3f_b%.3f_seed_stats.mat', ...
        P_VAL, lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta));
    save(stats_path, 'SEED_LIST', 'P_VAL', 'ETA_LIST', ...
        'SF_REPEATS', 'SB_REPEATS', 'DELTA_REPEATS', 'DELTA_MEAN', ...
        'DELTA_STD', 'DELTA_ZERO_TOL', 'P_ZERO');
end

%% ============================================================
%% 2. 统一绘图区 (在同一图窗中对比展示所有 p)
%% ============================================================
fprintf('\n[PLOT] 正在生成多 p 合一对比图...\n');
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.52, 0.48]);
ax = axes('Position', [0.15, 0.17, 0.72, 0.75]);
hold(ax, 'on');

% 颜色与标记配置
colors = lines(num_p);
markers = {'o', 's', '^', 'v', 'd', 'p', 'h'};
h_mean_list = gobjects(num_p, 1);
h_zero_list = gobjects(num_p, 1);

% --- 左 Y 轴：平均滞回宽度及其误差带 ---
yyaxis(ax, 'left');
for p_idx = 1:num_p
    c = colors(p_idx, :);
    d_mean = ALL_DELTA_MEAN(p_idx, :);
    d_std  = ALL_DELTA_STD(p_idx, :);

    % 绘制半透明误差带 (阴影)
    fill(ax, [ETA_LIST, fliplr(ETA_LIST)], ...
        [d_mean + d_std, fliplr(d_mean - d_std)], c, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');

    % 绘制均值主实线
    h_mean_list(p_idx) = plot(ax, ETA_LIST, d_mean, '-', ...
        'Color', c, 'LineWidth', 2.2, ...
        'DisplayName', sprintf('$\\langle\\Delta\\sigma\\rangle$ ($p=%.2f$)', P_LIST(p_idx)));
end
plot(ax, [ETA_LIST(1), ETA_LIST(end)], [0, 0], 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
ylabel(ax, '$\langle\Delta\sigma\rangle$', 'FontSize', 16, 'Interpreter', 'latex');

% --- 右 Y 轴：近零概率 ---
yyaxis(ax, 'right');
for p_idx = 1:num_p
    c = colors(p_idx, :);
    m = markers{mod(p_idx - 1, length(markers)) + 1};
    p_zero = ALL_P_ZERO(p_idx, :);

    % 绘制带标记的虚线
    h_zero_list(p_idx) = plot(ax, ETA_LIST, p_zero, '--', ...
        'Color', c, 'LineWidth', 1.8, 'Marker', m, 'MarkerSize', 5, ...
        'MarkerIndices', (1 + (p_idx-1)*3):15:num_e, ...
        'DisplayName', sprintf('$P(|\\Delta\\sigma|\\leq %.2f)$ ($p=%.2f$)', DELTA_ZERO_TOL, P_LIST(p_idx)));
end
ylabel(ax, sprintf('$P(|\\Delta\\sigma|\\leq %.2f)$', DELTA_ZERO_TOL), ...
    'FontSize', 16, 'Interpreter', 'latex');
ylim(ax, [0, 1.05]);
yticks(ax, 0:0.2:1);

% 坐标轴公共修饰
xlabel(ax, '$\eta$', 'FontSize', 16, 'Interpreter', 'latex');
xlim(ax, [0, ETA_LIST(end)]);
grid(ax, 'on');
box(ax, 'on');
set(ax, 'FontSize', 14, 'TickLabelInterpreter', 'latex');

% 图例设置 (包含均值曲线与概率曲线)
legend(ax, [h_mean_list; h_zero_list], 'Location', 'best', ...
    'Interpreter', 'latex', 'NumColumns', 2, 'FontSize', 11);

% 面板标签
annotation(h, 'textbox', [0.03, 0.90, 0.06, 0.07], 'String', '(b)', ...
    'FontSize', 16, 'BackgroundColor', 'none', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% 保存合一图片
test_out_img = fullfile(test_plots_dir, 'fig1b_seed_combined.png');
exportgraphics(h, test_out_img, 'Resolution', 300);
fprintf('[DONE] 多 p 合一对比图已保存: %s\n', test_out_img);
