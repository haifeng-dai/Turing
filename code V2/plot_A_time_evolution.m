%% plot_A_time_evolution.m: 同一噪声下不同 sigma 的序参量时间演化
clear; clc; close all;

%% 1. 路径与参数配置
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(script_dir, 'simulations'), fullfile(script_dir, 'networks'));

TOPO_TYPE = 'ER';
TOPO_PARAM = 0.01;  % ER：连接概率 p；WS：重连概率；BA：参数 m；SF：幂律指数 gamma

cfg.N = 200;
cfg.K = 5;
cfg.alpha = 0.05;
cfg.beta = 0.005;
cfg.noise = 1;
cfg.T_END = 1000;
cfg.dt = 0.005;
cfg.steps = 2001;   % 保存状态数，包含 t = 0 和 T_END
cfg.init_perturb = 0.1;
cfg.detect_convergence = false;  % 本脚本需要完整时间轨迹，暂时关闭收敛提前停止

SIGMA_LIST = 20.237:0.0001:20.2375;
RUN_SIGMA_INDICES = 1:numel(SIGMA_LIST);  % 默认运行全部；可改为 [1 5 10] 等序号子集
FORCE_RERUN = false;  % true：重算所选 sigma；false：优先使用已有缓存
INITIAL_SEED = 1;
cfg.noise_seed = 1;  % 使用同一噪声随机种子，保证各组噪声序列一致

if isempty(RUN_SIGMA_INDICES) || any(RUN_SIGMA_INDICES < 1) || ...
        any(RUN_SIGMA_INDICES > numel(SIGMA_LIST)) || ...
        any(RUN_SIGMA_INDICES ~= fix(RUN_SIGMA_INDICES))
    error('RUN_SIGMA_INDICES 中的序号必须是 SIGMA_LIST 范围内的正整数。');
end

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

%% 3. 加载缓存并运行所选 sigma
% 先准备与当前网络和模型参数一致的无噪声 Turing 初态。
results_dir = fullfile(script_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
topology_info = dir(topology_file);

seed_file = fullfile(results_dir, sprintf( ...
    'A_time_turing_seed_%s_N%d_K%d_topo%.8g_netbytes%d_netmtime%.10f_a%.6g_b%.6g_T%.6g_dt%.6g_steps%d_init%.6g_iseed%d.mat', ...
    lower(TOPO_TYPE), cfg.N, cfg.K, TOPO_PARAM, topology_info.bytes, ...
    topology_info.datenum, cfg.alpha, cfg.beta, cfg.T_END, cfg.dt, ...
    cfg.steps, cfg.init_perturb, INITIAL_SEED));

if exist(seed_file, 'file')
    seed_data = load(seed_file, 'turing_seed');
    turing_seed = seed_data.turing_seed;
else
    fprintf('[种子] 未找到匹配的 Turing 初态，开始生成 sigma=100、无噪声种子。\n');
    seed_cfg = cfg;
    seed_cfg.sigma = 100;
    seed_cfg.noise = 0;
    seed_cfg.detect_convergence = true;
    seed_cfg.y0 = [];
    rng(INITIAL_SEED, 'twister');
    [~, Y_seed, ~] = solve_multiplex(seed_cfg);
    turing_seed = Y_seed(end, :)';

    seed_u = turing_seed(1:cfg.N * cfg.K);
    seed_v = turing_seed(cfg.N * cfg.K + 1:end);
    seed_A = sqrt(sum((seed_u - 5).^2 + (seed_v - 10).^2) / (cfg.N * cfg.K));
    if ~isfinite(seed_A) || seed_A <= 0.02
        error('sigma=100 生成的种子未形成明显 Turing 态（A=%.8g）。', seed_A);
    end

    temp_seed_file = strrep(seed_file, '.mat', '_tmp.mat');
    save(temp_seed_file, 'turing_seed', 'seed_A', '-v7');
    [saved, message] = movefile(temp_seed_file, seed_file, 'f');
    if ~saved
        error('Turing 初态种子保存失败：%s', message);
    end
    fprintf('[种子] 已保存 Turing 初态，A=%.8g：%s\n', seed_A, seed_file);
end

if numel(turing_seed) ~= 2 * cfg.N * cfg.K
    error('Turing 初态长度与当前 N、K 不匹配：%s', seed_file);
end
seed_info = dir(seed_file);

% 缓存按 sigma 和初态类型分别保存；每完成一条曲线就写盘，意外中断后可续跑。
cache_file = fullfile(results_dir, sprintf( ...
    'A_time_evolution_cache_%s_N%d_K%d_topo%.8g_netbytes%d_netmtime%.10f_seedbytes%d_seedmtime%.10f_a%.6g_b%.6g_n%.6g_T%.6g_dt%.6g_steps%d_init%.6g_iseed%d_nseed%d.mat', ...
    lower(TOPO_TYPE), cfg.N, cfg.K, TOPO_PARAM, ...
    topology_info.bytes, topology_info.datenum, seed_info.bytes, ...
    seed_info.datenum, cfg.alpha, cfg.beta, ...
    cfg.noise, cfg.T_END, ...
    cfg.dt, cfg.steps, cfg.init_perturb, ...
    INITIAL_SEED, cfg.noise_seed));

if exist(cache_file, 'file')
    cache_data = load(cache_file, 'result');
    result = cache_data.result;
else
    result.version = 2;
    result.sigma = [];
    result.A_time = cell(0, 2);  % 第 1 列：均匀态扰动；第 2 列：Turing 初态
    result.A_final = zeros(0, 2);
    result.t = [];
end

if ~isfield(result, 'version') || result.version ~= 2
    error('缓存格式不受支持，请移走或删除缓存文件后重试：%s', cache_file);
end

fprintf('本次选择运行 %d/%d 个 sigma，缓存文件：%s\n', ...
    numel(RUN_SIGMA_INDICES), numel(SIGMA_LIST), cache_file);
NK = cfg.N * cfg.K;

for run_idx = 1:numel(RUN_SIGMA_INDICES)
    i = RUN_SIGMA_INDICES(run_idx);
    sigma = SIGMA_LIST(i);
    cache_idx = find(abs(result.sigma - sigma) <= 1e-10 * max(1, abs(sigma)), 1);

    if isempty(cache_idx)
        cache_idx = numel(result.sigma) + 1;
        result.sigma(cache_idx) = sigma;
        result.A_time{cache_idx, 1} = [];
        result.A_time{cache_idx, 2} = [];
        result.A_final(cache_idx, 1:2) = NaN;
    end

    for initial_type = 1:2
        if ~FORCE_RERUN && ~isempty(result.A_time{cache_idx, initial_type})
            fprintf('[缓存] sigma=%.10g，初态=%s，跳过仿真。\n', ...
                sigma, initial_type_name(initial_type));
            continue;
        end

        fprintf('[仿真 %d/%d] sigma=%.10g，初态=%s\n', ...
            run_idx, numel(RUN_SIGMA_INDICES), sigma, initial_type_name(initial_type));
        cfg_run = cfg;
        cfg_run.sigma = sigma;
        if initial_type == 1
            cfg_run.y0 = [];
        else
            cfg_run.y0 = turing_seed;
        end
        % 每次仿真均重置随机流，确保各组使用相同初始扰动和噪声序列。
        rng(INITIAL_SEED, 'twister');
        [t, Y, cfg_out] = solve_multiplex(cfg_run);

        u = Y(:, 1:NK);
        v = Y(:, NK + 1:2 * NK);
        deviation_squared = (u - 5).^2 + (v - 10).^2;
        A_current = sqrt(sum(deviation_squared, 2) / NK);
        result.A_time{cache_idx, initial_type} = A_current;
        result.A_final(cache_idx, initial_type) = cfg_out.A_final;
        result.t = t;

        temp_cache_file = strrep(cache_file, '.mat', '_tmp.mat');
        save(temp_cache_file, 'result', '-v7');
        [saved, message] = movefile(temp_cache_file, cache_file, 'f');
        if ~saved
            error('缓存文件写入失败：%s', message);
        end
        fprintf('[已保存] sigma=%.10g，初态=%s，A_final=%.8g\n', ...
            sigma, initial_type_name(initial_type), cfg_out.A_final);
    end
end

%% 4. 绘制并保存两种初态下所选 sigma 的 A(t)
h = figure('Visible', 'off', 'Color', 'w', 'Name', '不同 sigma 和初态下的 A(t) 演化');
tiledlayout(numel(RUN_SIGMA_INDICES), 2, 'TileSpacing', 'compact');

for run_idx = 1:numel(RUN_SIGMA_INDICES)
    i = RUN_SIGMA_INDICES(run_idx);
    sigma = SIGMA_LIST(i);
    cache_idx = find(abs(result.sigma - sigma) <= 1e-10 * max(1, abs(sigma)), 1);
    if isempty(cache_idx)
        error('sigma=%.10g 没有仿真结果或缓存，请检查 RUN_SIGMA_INDICES。', sigma);
    end

    if isempty(result.A_time{cache_idx, 1}) || isempty(result.A_time{cache_idx, 2})
        error('sigma=%.10g 的两种初态结果未齐全，请重新运行所选序号。', sigma);
    end

    for initial_type = 1:2
        ax = nexttile;
        plot(result.t, result.A_time{cache_idx, initial_type}, 'LineWidth', 1.6);
        xlabel('t', 'FontSize', 14);
        ylabel('$A(t)$', 'FontSize', 14, 'Interpreter', 'latex');
        title(sprintf('\\sigma=%.10g，%s', sigma, initial_type_name(initial_type)), ...
            'FontSize', 14);
        grid on; box on;
        set(ax, 'FontSize', 14);
    end
end

fig_dir = fullfile(script_dir, 'fig');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
out_img = fullfile(fig_dir, sprintf( ...
    'A_time_evolution_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_T%.0f.png', ...
    lower(TOPO_TYPE), cfg.N, cfg.K, TOPO_PARAM, cfg.alpha, cfg.beta, ...
    cfg.noise, cfg.T_END));
exportgraphics(h, out_img, 'Resolution', 300);

for i = 1:length(SIGMA_LIST)
    cache_idx = find(abs(result.sigma - SIGMA_LIST(i)) <= ...
        1e-10 * max(1, abs(SIGMA_LIST(i))), 1);
    if ~isempty(cache_idx)
        fprintf('sigma=%.10g 时，均匀态初值 A_final=%.8g，Turing 初值 A_final=%.8g\n', ...
            SIGMA_LIST(i), result.A_final(cache_idx, 1), result.A_final(cache_idx, 2));
    end
end
fprintf('仿真参数：拓扑=%s，N=%d，K=%d，alpha=%.4g，beta=%.4g，噪声=%.4g\n', ...
    TOPO_TYPE, cfg.N, cfg.K, cfg.alpha, cfg.beta, cfg.noise);
fprintf('[完成] 图像已保存：%s\n', out_img);

% 本脚本关闭普通轨迹的自动收敛提前停止；生成 Turing 种子时允许收敛检测。

function name = initial_type_name(initial_type)
if initial_type == 1
    name = '均匀态扰动初值';
else
    name = 'Turing 初态';
end
end
