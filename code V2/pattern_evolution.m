%% pattern_evolution.m: 图样演化保存版 (中英混合对齐)
clear; clc; close all;
addpath('simulations', 'networks');

RUN_SIMULATION = true;  % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)
USE_SEED       = true;  % 控制初值模式 (true: 从斑图种子开始/反向模拟, false: 从随机噪声开始/正向模拟)

%% 1. 实验参数配置 (核心动力学参数)
TOPO_TYPE  = 'ER';
TOPO_PARAM = 0.03;

DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta  = 0.1 * DYNA.alpha;
DYNA.noise = 0.80;
DYNA.T_END = 500;
DYNA.steps = 50;
DYNA.init_perturb = 0.1;
DYNA.sigma = 16.0;

%% 2. 文件路径定义 & 仿真执行区
% 自动映射分支扫描方向
if USE_SEED, scan_mode = 'bwd'; else, scan_mode = 'fwd'; end

mat_name = sprintf('evolution_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_s%.1f_n%.2f_%s_results.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, DYNA.sigma, DYNA.noise, scan_mode);
out_path = fullfile(fileparts(mfilename('fullpath')), 'results', mat_name);

if RUN_SIMULATION
    % --- A. 动态加载层内拓扑库 (适配扁平化映射) ---
    topo_base = fullfile(fileparts(mfilename('fullpath')), 'results', 'topology');
    type_u = upper(TOPO_TYPE);
    switch type_u
        case 'ER', net_file = fullfile(topo_base, 'ER', sprintf('N%d_p%.3f.mat', DYNA.N, TOPO_PARAM));
        case 'WS', net_file = fullfile(topo_base, 'WS', sprintf('N%d_pr%.2f.mat', DYNA.N, TOPO_PARAM));
        case 'BA', net_file = fullfile(topo_base, 'BA', sprintf('N%d_m%d.mat', DYNA.N, TOPO_PARAM));
        case 'SF', net_file = fullfile(topo_base, 'SF', sprintf('N%d_g%.1f.mat', DYNA.N, TOPO_PARAM));
        otherwise, error('未知拓扑类型: %s', TOPO_TYPE);
    end

    if ~exist(net_file, 'file'), error('找不到网络文件: %s', net_file); end
    data = load(net_file, 'nets');
    DYNA.L_intra = data.nets; % 此处正式填入 L_intra 字段

    % --- B. 构造层间拉普拉斯算子 L_inter ---
    adj_inter = ones(DYNA.K) - eye(DYNA.K);
    DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

    if USE_SEED
        fprintf('[PRE] 正在加载全局强斑图种子 (sigma=100) 作为演化初值 (反向模式)...\n');
        % 使用与 sweep_param.m 相同的标准种子路径逻辑
        seed_name = sprintf('evolution_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_s100.0_n0.00_fwd_results.mat', ...
            lower(TOPO_TYPE), DYNA.N, DYNA.K, 0.03, 0.05, 0.005);
        seed_file = fullfile(fileparts(mfilename('fullpath')), 'results', seed_name);

        if ~exist(seed_file, 'file')
            error('找不到标准种子文件: %s\n请先设置 USE_SEED=false, sigma=100.0 运行以生成种子。', seed_file);
        end
        tmp_seed = load(seed_file);
        DYNA.y0 = tmp_seed.Y(end, :)';
    else
        fprintf('[PRE] 从随机噪声开始模拟 (正向模式)...\n');
        DYNA.y0 = []; % 触发 solve_multiplex 内部的随机初值生成
    end

    %% 4. 执行单点演化仿真 (主任务)
    fprintf('[SIM] 正在执行单点演化仿真 (sigma=%.2f)...\n', DYNA.sigma);
    [t, Y, cfg] = solve_multiplex(DYNA); % 核心解算

    % 根据第一层节点的度进行排序 (用于后续轨迹分析)
    degs = -diag(DYNA.L_intra{1});
    [~, sort_idx] = sort(degs, 'descend');

    % --- 结果持久化存档 ---
    save(out_path, 't', 'Y', 'cfg', 'sort_idx');
    fprintf('[DONE] 仿真分析已完成并存入: %s\n', mat_name);
end

%% 5. 结果处理与数据还原
% 加载刚才存入的结果 (保持索引完全同步)
load(out_path);

% 提取维度信息
N = cfg.N; K = size(cfg.L_inter, 1); NK = N * K;
U_all = Y(:, 1:NK);
V_all = Y(:, NK+1:end);

% 计算包含 u 和 v 的综合序参数演化 (L2-norm)
A_t = sqrt(sum((U_all - 5).^2 + (V_all - 10).^2, 2) / (N * K));

% 颜色定义
c_u = [0.0 0.447 0.741]; % 深蓝
c_v = [0.85 0.325 0.098]; % 深橙

plots_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

%% 6. 绘图: 各层节点演化轨迹 (前 5 层示意)
% f1 = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.45, 0.8]);
% K_show = min(K, 5);
% for i = 1:K_show
%     % 绘制 u 分量演化
%     subplot(K_show, 2, 2*i-1); hold on;
%     layer_u = U_all(:, (i-1)*N+1 : i*N);
%     u_end = layer_u(end, :);
%     [~, idx_up] = max(u_end);   % 找上分支代表
%     [~, idx_down] = min(u_end); % 找下分支代表

%     plot(t, layer_u(:, 1:5:end), 'Color', [0.85 0.85 0.85]); % 背景线
%     plot(t, layer_u(:, idx_up), 'Color', c_u, 'LineWidth', 1.5);    % 上分支高亮
%     plot(t, layer_u(:, idx_down), 'Color', c_u*0.6, 'LineWidth', 1.5); % 下分支高亮(颜色稍深)
%     ylabel(['L', num2str(i), ' (u)'], 'FontSize', 10);
%     grid on; set(gca, 'FontSize', 9);
%     % if i == 1, title('Evolution of u'); end

%     % 绘制 v 分量演化
%     subplot(K_show, 2, 2*i); hold on;
%     layer_v = V_all(:, (i-1)*N+1 : i*N);
%     v_end = layer_v(end, :);
%     [~, idx_up_v] = max(v_end);
%     [~, idx_down_v] = min(v_end);

%     plot(t, layer_v(:, 1:5:end), 'Color', [0.85 0.85 0.85]); % 背景线
%     plot(t, layer_v(:, idx_up_v), 'Color', c_v, 'LineWidth', 1.5);    % 上分支高亮
%     plot(t, layer_v(:, idx_down_v), 'Color', c_v*0.6, 'LineWidth', 1.5); % 下分支高亮
%     ylabel(['L', num2str(i), ' (v)'], 'FontSize', 10);
%     grid on; set(gca, 'FontSize', 9);
%     % if i == 1, title('Evolution of v'); end
% end
% xlabel('Time t');
% 生成动态文件名
fname_base = sprintf('evolution_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_s%.1f_n%.2f_%s', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, DYNA.sigma, DYNA.noise, scan_mode);

out_traj = fullfile(plots_dir, [fname_base, '_trajectories.png']);
out_order = fullfile(plots_dir, [fname_base, '_order.png']);
out_snapshot = fullfile(plots_dir, [fname_base, '_snapshot.png']);

% exportgraphics(f1, out_traj, 'Resolution', 300);

%% 7. 绘图: 全系统序参数 A(t) 的时间响应
f2 = figure('Units', 'normalized', 'Position', [0.56, 0.55, 0.3, 0.3]);
plot(t, A_t, 'Color', c_u, 'LineWidth', 2);
xlabel('Time t'); ylabel('A(t)');
% title('Evolution of Order Parameter A');
grid on;
exportgraphics(f2, out_order, 'Resolution', 300);

return;
%% 8. 绘图: 各层稳态空间分布快照 (Snapshot)
f3 = figure('Units', 'normalized', 'Position', [0.56, 0.1, 0.35, 0.4]);
for i = 1:K_show
    % u 的节点分布
    subplot(K_show, 2, 2*i-1);
    layer_u_end = U_all(end, (i-1)*N+1 : i*N);

    % 计算分支分布
    n_up = sum(layer_u_end > 5); n_down = N - n_up;

    scatter(1:N, layer_u_end, 15, layer_u_end, 'filled');
    colormap(gca, 'parula');
    cb = colorbar;
    cb.TickLabels = cellstr(num2str(cb.Ticks', '%.1f'));
    ytickformat('%.1f');
    ylabel(['L', num2str(i), ' (u)']); grid on;
    % title(sprintf('L%d u (Up:%d, Down:%d)', i, n_up, n_down), 'FontSize', 10);

    % v 的节点分布
    subplot(K_show, 2, 2*i);
    layer_v_end = V_all(end, (i-1)*N+1 : i*N);

    % 计算分支分布
    n_up_v = sum(layer_v_end > 10); n_down_v = N - n_up_v;

    scatter(1:N, layer_v_end, 15, layer_v_end, 'filled');
    colormap(gca, 'spring');
    cb = colorbar;
    cb.TickLabels = cellstr(num2str(cb.Ticks', '%.1f'));
    ytickformat('%.1f');
    ylabel(['L', num2str(i), ' (v)']); grid on;
    % title(sprintf('L%d v (Up:%d, Down:%d)', i, n_up_v, n_down_v), 'FontSize', 10);
end
xlabel('Node ID (Original)');
exportgraphics(f3, out_snapshot, 'Resolution', 300);

fprintf('\n[SUCCESS] 全部演化分析图表已生成：\n');
fprintf('  1. 分层轨迹图: %s\n', out_traj);
fprintf('  2. 序参数演化: %s\n', out_order);
fprintf('  3. 末态分布图: %s\n', out_snapshot);

%% 9. 额外绘图: 每一层的深入观察 (Focused View per Layer)
fprintf('\n[INFO] 正在生成各单层焦点分析图 (L1 - L%d)...\n', K);

for k = 1:K
    % --- A. 焦点轨线图 (u上 v下) ---
    f4 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1+0.02*k, 0.25, 0.35, 0.45]);

    % u 分量
    subplot(2, 1, 1); hold on;
    set(gca, 'FontSize', 14);
    layer_u = U_all(:, (k-1)*N+1 : k*N);
    u_end = layer_u(end, :); [~, idx_up] = max(u_end); [~, idx_down] = min(u_end);
    plot(t, layer_u(:, 1:5:end), 'Color', [0.85 0.85 0.85]);
    plot(t, layer_u(:, idx_up), 'Color', c_u, 'LineWidth', 1.5);
    plot(t, layer_u(:, idx_down), 'Color', c_u*0.6, 'LineWidth', 1.5);
    ylabel(['L', num2str(k), ' (u)'], 'FontSize', 14);
    % title(sprintf('Focused Trajectories: Layer %d (\\sigma=%.1f)', k, DYNA.sigma), 'FontSize', 14);
    grid on; box on;

    % v 分量
    subplot(2, 1, 2); hold on;
    set(gca, 'FontSize', 14);
    layer_v = V_all(:, (k-1)*N+1 : k*N);
    v_end = layer_v(end, :); [~, idx_up_v] = max(v_end); [~, idx_down_v] = min(v_end);
    plot(t, layer_v(:, 1:5:end), 'Color', [0.85 0.85 0.85]);
    plot(t, layer_v(:, idx_up_v), 'Color', c_v, 'LineWidth', 1.5);
    plot(t, layer_v(:, idx_down_v), 'Color', c_v*0.6, 'LineWidth', 1.5);
    ylabel(['L', num2str(k), ' (v)'], 'FontSize', 14); xlabel('Time t', 'FontSize', 14);
    grid on; box on;

    out_focus_traj = fullfile(plots_dir, [fname_base, sprintf('_L%d_focus_traj.png', k)]);
    exportgraphics(f4, out_focus_traj, 'Resolution', 300);

    % --- B. 焦点快照图 (u上 v下) ---
    f5 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.48+0.02*k, 0.25, 0.35, 0.45]);

    % u 空间分布
    subplot(2, 1, 1);
    set(gca, 'FontSize', 14);
    layer_u_end = U_all(end, (k-1)*N+1 : k*N);
    n_up = sum(layer_u_end > 5); n_down = N - n_up;
    scatter(1:N, layer_u_end, 20, layer_u_end, 'filled'); colormap(gca, 'parula'); colorbar;
    ylabel('u state', 'FontSize', 14); xlabel('Node ID', 'FontSize', 14); grid on;
    % title(sprintf('L%d u (Up:%d, Down:%d)', k, n_up, n_down), 'FontSize', 14);

    % v 空间分布
    subplot(2, 1, 2);
    set(gca, 'FontSize', 14);
    layer_v_end = V_all(end, (k-1)*N+1 : k*N);
    n_up_v = sum(layer_v_end > 10); n_down_v = N - n_up_v;
    scatter(1:N, layer_v_end, 20, layer_v_end, 'filled'); colormap(gca, 'spring'); colorbar;
    ylabel('v state', 'FontSize', 14); xlabel('Node ID', 'FontSize', 14); grid on;
    % title(sprintf('L%d v (Up:%d, Down:%d)', k, n_up_v, n_down_v), 'FontSize', 14);

    out_focus_snap = fullfile(plots_dir, [fname_base, sprintf('_L%d_focus_snap.png', k)]);
    exportgraphics(f5, out_focus_snap, 'Resolution', 300);

    fprintf('  [DONE] 层 %d: \t轨线(%s) \t快照(%s)\n', k, ...
        sprintf('L%d_focus_traj.png', k), sprintf('L%d_focus_snap.png', k));
end
