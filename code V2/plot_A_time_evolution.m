%% plot_A_time_evolution.m: 同一噪声下不同 sigma 的序参量时间演化
clear; clc; close all;

%% 1. 路径与参数配置
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'simulations'), fullfile(script_dir, 'networks'));

TOPO_TYPE = 'ER';
TOPO_PARAM = 0.01;  % ER：连接概率 p；WS：重连概率；BA：参数 m；SF：幂律指数 gamma

cfg.N = 200;
cfg.K = 5;
cfg.alpha = 0.05;
cfg.beta = 0.005;
cfg.noise = 1;
cfg.T_END = 500;
cfg.dt = 0.005;
cfg.steps = 2001;   % 保存状态数，包含 t = 0 和 T_END
cfg.init_perturb = 0.1;
cfg.detect_convergence = false;  % 本脚本需要完整时间轨迹，暂时关闭收敛提前停止

SIGMA_LIST = 10:0.1:11;  % 覆盖图 1b 的 sigma 扫描范围，并保留原设定值 16
INITIAL_SEED = 1;
cfg.noise_seed = 1;  % 使用同一噪声随机种子，保证各组噪声序列一致

%% 2. 加载网络并构造层间耦合
topology_dir = fullfile(script_dir, 'results', 'topology', upper(TOPO_TYPE));
switch upper(TOPO_TYPE)
    case 'ER'
        topology_file = fullfile(topology_dir, ...
            sprintf('N%d_p%.3f.mat', cfg.N, TOPO_PARAM));
    case 'WS'
        topology_file = fullfile(topology_dir, ...
            sprintf('N%d_pr%.2f.mat', cfg.N, TOPO_PARAM));
    case 'BA'
        topology_file = fullfile(topology_dir, ...
            sprintf('N%d_m%d.mat', cfg.N, TOPO_PARAM));
    case 'SF'
        topology_file = fullfile(topology_dir, ...
            sprintf('N%d_g%.1f.mat', cfg.N, TOPO_PARAM));
    otherwise
        error('不支持的拓扑类型：%s', TOPO_TYPE);
end

if ~exist(topology_file, 'file')
    error('找不到拓扑文件：%s', topology_file);
end

topology_data = load(topology_file, 'nets');
cfg.L_intra = topology_data.nets(1:cfg.K);
adj_inter = ones(cfg.K) - eye(cfg.K);
cfg.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

%% 3. 运行不同 sigma 下的仿真并计算 A(t)
% 每次仿真均重置初始条件和噪声种子，使用相同的随机扰动与噪声序列。
NK = cfg.N * cfg.K;
A_time = zeros(cfg.steps, length(SIGMA_LIST));
A_final = zeros(size(SIGMA_LIST));

for i = 1:length(SIGMA_LIST)
    cfg.sigma = SIGMA_LIST(i);
    rng(INITIAL_SEED, 'twister');
    [t, Y, cfg_out] = solve_multiplex(cfg);

    u = Y(:, 1:NK);
    v = Y(:, NK + 1:2 * NK);
    deviation_squared = (u - 5).^2 + (v - 10).^2;
    A_time(:, i) = sqrt(sum(deviation_squared, 2) / NK);
    A_final(i) = cfg_out.A_final;
end

%% 4. 绘制并保存不同 sigma 的 A(t)
h = figure('Visible', 'off', 'Color', 'w', 'Name', '不同 sigma 下的 A(t) 演化', ...
    'Units', 'normalized', 'Position', [0.2, 0.1, 0.5, 0.8]);
tiledlayout(length(SIGMA_LIST), 1, 'TileSpacing', 'compact');

for i = 1:length(SIGMA_LIST)
    ax = nexttile;
    plot(t, A_time(:, i), 'LineWidth', 1.6);
    xlabel('t', 'FontSize', 14);
    ylabel('$A(t)$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('序参量演化（\\sigma = %.4g）', SIGMA_LIST(i)), 'FontSize', 14);
    grid on; box on;
    set(ax, 'FontSize', 14);
end

fig_dir = fullfile(script_dir, 'fig');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
out_img = fullfile(fig_dir, sprintf( ...
    'A_time_evolution_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f.png', ...
    lower(TOPO_TYPE), cfg.N, cfg.K, TOPO_PARAM, cfg.alpha, cfg.beta, cfg.noise));
exportgraphics(h, out_img, 'Resolution', 300);

for i = 1:length(SIGMA_LIST)
    fprintf('sigma=%.4g 时的 A_final（末段时间窗均值，归一化）：%.8g\n', ...
        SIGMA_LIST(i), A_final(i));
end
fprintf('仿真参数：拓扑=%s，N=%d，K=%d，alpha=%.4g，beta=%.4g，噪声=%.4g\n', ...
    TOPO_TYPE, cfg.N, cfg.K, cfg.alpha, cfg.beta, cfg.noise);
fprintf('[完成] 图像已保存：%s\n', out_img);

% 本脚本关闭了自动收敛提前停止；NaN/Inf 发散检查仍由求解器执行。
