%% sigma_eta_connectivity_er.m: 噪声强度(eta)影响分析 (新版 find_thresholds 适配版)
clear; clc; close all;
addpath('simulations', 'networks');

RUN_SIMULATION = true;  % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)

%% 1. 实验参数配置
TOPO_TYPE = 'ER';
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta  = 0.1 * DYNA.alpha;
DYNA.noise = 0.01;      % 种子预热时的名义噪声背景
DYNA.T_END = 500;
DYNA.steps = 2;         % 同样的极速配置
DYNA.init_perturb = 0.1;

% 搜索范围设置
DYNA.sigma_min = 0;
DYNA.sigma_max = 30;

P_LIST   = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1];
ETA_LIST = linspace(0, 3, 301);

% 数据保存路径 (动态包含 beta 强度以防覆盖)
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end
data_path = fullfile(res_dir, sprintf('scan_eta_p_%s_N%d_K%d_a%.3f_b%.3f.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta));

%% 2. 核心并行计算区
num_p = length(P_LIST);
num_e = length(ETA_LIST);
if RUN_SIMULATION
    [E_idx, P_idx] = meshgrid(1:num_e, 1:num_p);
    E_idx = E_idx(:);
    P_idx = P_idx(:);
    total_tasks = length(E_idx);

    Sigma_F_Results = zeros(total_tasks, 1);
    Sigma_B_Results = zeros(total_tasks, 1);

    % --- A. 构造全局层间算子 ---
    adj_inter = (ones(DYNA.K) - eye(DYNA.K));
    DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

    % --- B. 【核心改进】加载唯一的全局强斑图种子 (由 pattern_evolution.m 生成) ---
    fprintf('[PRE] 正在加载全局强斑图种子 (sigma=100, p=0.03)...\n');
    seed_file = fullfile('results', 'evolution_er_N200_K5_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat');
    if ~exist(seed_file, 'file'), error('找不到全局种子文件，请先运行 pattern_evolution.m！'); end
    tmp_seed = load(seed_file);
    y_universal_seed = tmp_seed.Y(end, :)';

    fprintf('[PRE] 正在并行预加载 %d 组拓扑...\n', num_p);
    NET_LIBS = cell(num_p, 1);
    parfor p_i = 1:num_p
        p_val = P_LIST(p_i);
        net_path = fullfile('results', 'topology', 'ER', sprintf('N%d_p%.3f.mat', DYNA.N, p_val));
        tmp = load(net_path);
        NET_LIBS{p_i} = tmp.nets;
    end

    % --- C. 批量执行二分查找探测 ---
    fprintf('\n[SIM] 开始并行扫描 (eta x connectivity): %d 组任务... \n', total_tasks);
    simStart = tic;

    parfor t = 1:total_tasks
        loopStart = tic;
        dyna_local = DYNA;
        dyna_local.noise = ETA_LIST(E_idx(t));
        cur_p_idx = P_idx(t);
        p_val = P_LIST(cur_p_idx);

        % 手动注入网络和全局强种子，最大化压榨 find_thresholds 的性能
        dyna_local.L_intra = NET_LIBS{cur_p_idx};
        dyna_local.y_seed  = y_universal_seed; % 使用全局唯一的强种子

        res = find_thresholds(TOPO_TYPE, p_val, dyna_local);
        Sigma_F_Results(t) = res(1);
        Sigma_B_Results(t) = res(2);

        fprintf('  [PAR] 任务 %d/%d: p=%.3f, eta=%.3f 计算完成，耗时: %.2f 秒。\n', t, total_tasks, p_val, dyna_local.noise, toc(loopStart));
    end

    simTime = toc(simStart);
    fprintf('[DONE] 仿真任务全量完成！耗时: %.2f 秒。\n', simTime);

    % 持久化保存
    SF_Matrix = reshape(Sigma_F_Results, [num_p, num_e]);
    SB_Matrix = reshape(Sigma_B_Results, [num_p, num_e]);
    save(data_path, 'SF_Matrix', 'SB_Matrix', 'P_LIST', 'ETA_LIST');
end

%% 3. 数据加载与对比绘图
if exist(data_path, 'file'), load(data_path); end

h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.45]);
ax = axes('Position', [0.18, 0.18, 0.75, 0.75]); hold on;

colors = lines(num_p);
for i = 1:num_p
    % 正向扫描
    plot(ETA_LIST, SF_Matrix(i, :), '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('p=%.2f (\\sigma_f)', P_LIST(i)));
    % 反向扫描
    plot(ETA_LIST, SB_Matrix(i, :), '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
        'MarkerSize', 6, 'DisplayName', sprintf('p=%.2f (\\sigma_b)', P_LIST(i)));
end

% title('Impact of Noise Intensity \eta on Hysteresis Thresholds', 'FontSize', 12);
xlabel('Noise Intensity \eta', 'FontSize', 14);
ylabel('Critical Diffusion Ratio \sigma', 'FontSize', 14);
% xlim(ax, [0, 1])
legend('Location', 'bestoutside', 'FontSize', 14);
grid on; set(ax, 'Box', 'on', 'FontSize', 14);

plots_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
out_img = fullfile(plots_dir, sprintf('scan_eta_p_%s_N%d_K%d_a%.3f_b%.3f.png', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta));
exportgraphics(h, out_img, 'Resolution', 300);
fprintf('[DONE] 绘图已更新: %s\n', out_img);
