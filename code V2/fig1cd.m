%% compare_noise_evolution_u.m: 专属定制脚本，对比不同噪声(eta=0 vs eta=0.3)下L1层的u演化与稳态分布
% 提取自 pattern_evolution.m 逻辑，仅绘制第一层 u 分量的轨迹和末态分布快照图
clear; clc; close all;
addpath('simulations', 'networks');

% 1. 设置实验参数 (与数据文件匹配)
TOPO_TYPE  = 'ER';
TOPO_PARAM = 0.03;
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta  = 0.005;
DYNA.sigma = 16.0;
scan_mode = 'bwd'; % 控制初始种子方向

noise_list = [0.00, 0.30];
U_data = cell(2, 1);
t_data = cell(2, 1);

results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
c_u = [0.0 0.447 0.741]; % 深蓝颜色

%% 2. 加载两次仿真的数据 (0噪声 与 0.3噪声)
fprintf('[INFO] 正在尝试加载数据...\n');
for i = 1:2
    cur_noise = noise_list(i);
    mat_name = sprintf('evolution_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_s%.1f_n%.2f_%s_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, DYNA.sigma, cur_noise, scan_mode);
    file_path = fullfile(results_dir, mat_name);
    fprintf("\n%s\n", mat_name)

    if ~exist(file_path, 'file')
        error('未找到数据文件: %s\n请先在 pattern_evolution.m 中分别将 DYNA.noise 设为 0.00 和 0.30，运行(RUN_SIMULATION=true)以生成数据！', mat_name);
    end

    data = load(file_path, 't', 'Y');
    t_data{i} = data.t;
    U_data{i} = data.Y(:, 1:DYNA.N); % 我们只取第1层 (1:N) 的 u
end

% 提取绘图核心线段
t1 = t_data{1}; u1 = U_data{1}; u1_end = u1(end, :);
[~, idx_up1] = max(u1_end); [~, idx_down1] = min(u1_end);

t2 = t_data{2}; u2 = U_data{2}; u2_end = u2(end, :);
[~, idx_up2] = max(u2_end); [~, idx_down2] = min(u2_end);

% 自动探测从上分支跳变到下分支的特定节点 (突变轨线)
% 判据：在仿真前半段处于高位(>4)，但在演化结束时处于低位(<2)
u2_early = u2(round(length(t2)*0.2), :); % 取演化 20% 处的采样
idx_switch = find(u2_early > 4 & u2_end < 2, 1);
if isempty(idx_switch)
    [~, idx_switch] = max(max(u2) - min(u2)); % 如果没找到严格跳变的，就找波动最大的那个
end

plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'V2', 'manuscirpt', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

%% 3. 画轨迹对比 (仅 L1 的 u)
fprintf('[PLOT] 正在生成轨迹对比图...\n');
f1 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.2, 0.35, 0.45]);
t1_lay = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% 上图: eta = 0.00
nexttile(t1_lay); hold on; set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');
plot(t1, u1(:, 1:5:end), 'Color', [0.85 0.85 0.85]);
plot(t1, u1(:, idx_up1), 'Color', c_u, 'LineWidth', 1.5);
plot(t1, u1(:, idx_down1), 'Color', c_u*0.6, 'LineWidth', 1.5);
ylabel('u(t)', 'FontSize', 18, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylim_val = ylim;
text(0.95, 0.8, '$\eta = 0.0$', 'Units', 'normalized', 'FontSize', 18, ...
    'Interpreter', 'latex', 'BackgroundColor', 'w', 'EdgeColor', 'w', 'Margin', 2, 'HorizontalAlignment', 'right');
xlim([0, 500]); grid on; box on;

% 下图: eta = 0.30
nexttile(t1_lay); hold on; set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');
plot(t2, u2(:, 1:5:end), 'Color', [0.85 0.85 0.85]);
plot(t2, u2(:, idx_up2), 'Color', c_u, 'LineWidth', 1.5);
plot(t2, u2(:, idx_down2), 'Color', c_u*0.6, 'LineWidth', 1.5);
% 绘制并突出跳变波形的轨线
if ~isempty(idx_switch)
    plot(t2, u2(:, idx_switch), 'Color', [0.9 0.4 0], 'LineWidth', 2.5); % 保持之前配色，但加粗突出
end
ylabel('$u(t)$', 'FontSize', 18, 'Interpreter', 'latex'); xlabel('$t$', 'FontSize', 18, 'Interpreter', 'latex');
ylim(ylim_val);
text(0.95, 0.8, '$\eta = 0.3$', 'Units', 'normalized', 'FontSize', 18, ...
    'Interpreter', 'latex', 'BackgroundColor', 'w', 'EdgeColor', 'w', 'Margin', 2, 'HorizontalAlignment', 'right');
xlim([0, 500]); grid on; box on;

out_traj = fullfile(plots_dir, 'fig1c.eps');
exportgraphics(f1, out_traj);

%% 4. 画末态节点分布快照对比 (仅 L1 的 u)
fprintf('[PLOT] 正在生成快照对比图...\n');
f2 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.55, 0.2, 0.35, 0.45]);
t2_lay = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% 预计算统一的色带范围
c_min = min([min(u1_end), min(u2_end)]);
c_max = max([max(u1_end), max(u2_end)]);

% 上图: eta = 0.00
ax1 = nexttile(t2_lay); hold on; set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');
scatter(1:DYNA.N, u1_end, 20, u1_end, 'filled');
colormap(gca, 'jet');
ylabel('$u(x)$', 'FontSize', 18, 'Interpreter', 'latex');
text(0.95, 0.3, '$\eta = 0.0$', 'Units', 'normalized', 'FontSize', 18, ...
    'Interpreter', 'latex', 'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 2, 'HorizontalAlignment', 'right');
ylim([min(u1_end), max(u1_end)]);
grid on; box on; clim(ax1, [c_min c_max]);

% 下图: eta = 0.30
ax2 = nexttile(t2_lay); hold on; set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');
scatter(1:DYNA.N, u2_end, 20, u2_end, 'filled');
colormap(gca, 'jet');
ylabel('$u(x)$', 'FontSize', 18, 'Interpreter', 'latex'); xlabel('Node ID', 'FontSize', 18, 'Interpreter', 'latex');
text(0.95, 0.3, '$\eta = 0.3$', 'Units', 'normalized', 'FontSize', 18, ...
    'Interpreter', 'latex', 'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 2, 'HorizontalAlignment', 'right');
ylim([min(u2_end), max(u2_end)]);
grid on; box on; clim(ax2, [c_min c_max]);

% 挂载统一共享色条
cb = colorbar(ax2);
cb.Layout.Tile = 'east';
% cb.Label.String = 'u state';
cb.Label.FontSize = 18;
cb.Label.Interpreter = 'latex';

out_snap = fullfile(plots_dir, 'fig1d.eps');
exportgraphics(f2, out_snap);

fprintf('\n[SUCCESS] 噪声对比图集 (L1的u分量) 已完成！\n');
fprintf('  => %s\n', out_traj);
fprintf('  => %s\n', out_snap);
