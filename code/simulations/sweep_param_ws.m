function sweep_param_ws(topo_type, k_val, p_rewire, cfg)
% SWEEP_PARAM_WS: WS 网络专用的参数扫描引擎（处理两个拓扑参数：平均度 K 和重连概率 p）

%% 1. 环境与路径初始化
rng(0, 'twister');                           % 锁定随机种子
topo_base = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'topology');
type_u = upper(topo_type);

%% 2. 为 WS 网络加载网络文件
% WS 网络文件命名规则：N{N}_K{k_val}_pr{p_rewire}.mat
switch type_u
    case 'WS'
        net_file = fullfile(topo_base, 'WS', sprintf('N%d_K%d_pr%.2f.mat', cfg.N, k_val, p_rewire));
    otherwise
        error('sweep_param_ws 专用于 WS 网络。若需其他拓扑，请使用 sweep_param。');
end

% 检查网络文件是否存在
if ~exist(net_file, 'file')
    error('[ERROR] 未找到 WS 网络文件: %s\n请先运行 generate_standard_networks.m 生成标准网络库！', net_file);
end

data = load(net_file, 'nets');
fprintf('[LOAD] 网络已从文件加载: %s\n', net_file);
cfg.L_intra = data.nets;

%% 3. 仿真实时参数配置
cfg.dt = 0.005;
cfg.steps = 2;

%% 4. 获取全局唯一的"斑图种子"
seed_name = sprintf('evolution_er_N%d_K%d_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', cfg.N, cfg.K);
seed_file = fullfile(fileparts(mfilename('fullpath')), '..', 'results', seed_name);

if ~exist(seed_file, 'file')
    error('【严重错误】未找到全局种子文件！请先运行 pattern_evolution.m 以生成标准参考种子。');
end

fprintf('[PRE] 正在加载全局强斑图种子 (sigma=100)...\n');
tmp_seed = load(seed_file);
y_seed = tmp_seed.Y(end, :)';

%% 5. 扫描序列生成与任务打包
sigma_range = linspace(cfg.sigma_min, cfg.sigma_max, cfg.sigma_npts);
npts = cfg.sigma_npts;
total_tasks = 2 * npts;

Sigma_Tasks = [sigma_range, fliplr(sigma_range)];
Y0_Tasks = cell(1, total_tasks);

for i = 1:npts
    Y0_Tasks{i} = [];
end

for i = npts+1:total_tasks
    Y0_Tasks{i} = y_seed;
end

A_results = zeros(1, total_tasks);

%% 6. 核心并行循环
simStart = tic;
parfor t = 1:total_tasks
    c = cfg;
    c.sigma = Sigma_Tasks(t);
    c.y0 = Y0_Tasks{t};

    [~, ~, cfg_out] = solve_multiplex(c);
    A_results(t) = cfg_out.A_final;
end

%% 7. 结果还原与数据拆解
A_fwd = A_results(1:npts);
A_bwd = fliplr(A_results(npts+1:end));

%% 8. 结果持久化与保存到磁盘
res_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end

results_name = sprintf('hysteresis_%s_N%d_K%d_pr%.2f_mK%d_a%.3f_b%.3f_n%.2f_results.mat', ...
    lower(type_u), cfg.N, cfg.K, p_rewire, k_val, cfg.alpha, cfg.beta, cfg.noise);
results_file = fullfile(res_dir, results_name);

save(results_file, 'sigma_range', 'A_fwd', 'A_bwd', 'npts', 'total_tasks');
fprintf('[DONE] 结果已保存: %s\n', results_file);

fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('[INFO] 任务总结\n');
fprintf('├─ 扫描点数: %d (双向 = %d)\n', npts, total_tasks);
fprintf('├─ 扫描范围: sigma ∈ [%.2f, %.2f]\n', cfg.sigma_min, cfg.sigma_max);
fprintf('├─ WS 拓扑: N=%d, K=%d, p_rewire=%.2f\n', cfg.N, k_val, p_rewire);
fprintf('├─ 动力学: alpha=%.3f, beta=%.3f, noise=%.2f\n', cfg.alpha, cfg.beta, cfg.noise);
fprintf('├─ 运行时间: %.2f 秒\n', toc(simStart));
fprintf('└─ 结果文件: %s\n', results_name);
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

end
