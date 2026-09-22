function sweep_param(topo_type, topo_val, cfg)
% SWEEP_PARAM: 能够识别滞后效应（Hysteresis）的高性能并行参数扫描引擎

%% 1. 环境与路径初始化
rng(0, 'twister');                           % 锁定随机种子，确保实验结果不同机器可重现
% 获取当前文件所在目录，定位到上级 results/topology 下的网络库
topo_base = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'topology');
type_u = upper(topo_type);                    % 将拓扑类型转为大写（如 'er' -> 'ER'）

%% 2. 拓扑库动态加载 (适配扁平化映射)
% 根据网络大小 N 和生成的概率 p（或 m 等）拼接对应的 .mat 文件名
switch type_u
    case 'ER', net_file = fullfile(topo_base, 'ER', sprintf('N%d_p%.3f.mat', cfg.N, topo_val));
    case 'WS', net_file = fullfile(topo_base, 'WS', sprintf('N%d_pr%.2f.mat', cfg.N, topo_val));
    case 'BA', net_file = fullfile(topo_base, 'BA', sprintf('N%d_m%d.mat', cfg.N, topo_val));
    case 'SF', net_file = fullfile(topo_base, 'SF', sprintf('N%d_g%.1f.mat', cfg.N, topo_val));
    otherwise, error('未知拓扑类型: %s', topo_type);
end

% 检查网络文件是否存在，不存在则终止报错
if ~exist(net_file, 'file'), error('找不到网络文件: %s', net_file); end
data = load(net_file, 'nets');                % 从磁盘加载保存好的网络结构
cfg.L_intra = data.nets;                      % 将拉普拉斯算子池存入配置结构体中

%% 3. 仿真实时参数配置
cfg.dt = 0.005;                               % 设置较大的步长，加快寻找稳态的速度
cfg.steps = 2;                                % 每个测试点仅取终态，节省计算内存

%% 4. 获取全局唯一的"斑图种子" (反向扫描的关键)
% 锁定策略：强制加载由 pattern_evolution.m 生成的全局种子，拒绝任何形式的实时训练以保证数据绝对对齐
seed_name = sprintf('evolution_er_N%d_K%d_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', cfg.N, cfg.K);
seed_file = fullfile(fileparts(mfilename('fullpath')), '..', 'results', seed_name);

if ~exist(seed_file, 'file')
    error('【严重错误】未找到全局种子文件！请先运行 pattern_evolution.m 以生成标准参考种子。');
end

fprintf('[PRE] 正在加载全局强斑图种子 (sigma=100)...\n');
tmp_seed = load(seed_file);
y_seed = tmp_seed.Y(end, :)';

%% 5. 扫描序列生成与任务打包
% 生成从小到大排列的扩散比序列序列
sigma_range = linspace(cfg.sigma_min, cfg.sigma_max, cfg.sigma_npts);
npts = cfg.sigma_npts;                         % 记录单向扫描点的总数
total_tasks = 2 * npts;                        % 合计任务点 = 正向扫描点 + 反向扫描点

% 构造包含正反双向的任务总表
Sigma_Tasks = [sigma_range, fliplr(sigma_range)];
% 为每个并行核心分配对应的初值 (Cell 数组存储)
Y0_Tasks = cell(1, total_tasks);

% 1..N: 正向点扫描任务：初值设为空（触发随机噪声，看什么时候能长出斑图）
for i = 1:npts
    Y0_Tasks{i} = [];
end

% N+1..2N: 反向点扫描任务：初值统一设定为斑图种子（看强斑图什么时候会崩塌）
for i = npts+1:total_tasks
    Y0_Tasks{i} = y_seed;
end
A_results = zeros(1, total_tasks);             % 预分配空间用于并行存储结果

%% 6. 核心并行循环 (采用高稳健直接打印模式)
simStart = tic;                                % 启动全局计时器
parfor t = 1:total_tasks
    c = cfg;                                   % 领取独立的配置副本
    c.sigma = Sigma_Tasks(t);                  % 配发该任务对应的扩散比参数
    c.y0 = Y0_Tasks{t};                        % 配发该任务对应的初值状态

    [~, ~, cfg_out] = solve_multiplex(c);               % 执行核心微分方程组求解

    % 取带时间窗口平均的最终全系统序参数 A
    A_results(t) = cfg_out.A_final;
end

%% 7. 结果还原与数据拆解
A_fwd = A_results(1:npts);                     % 对应 sigma_range 原序作为正向回线
A_bwd = fliplr(A_results(npts+1:end));         % 获取反向计算点，翻转回 sigma_range 的序

%% 8. 结果持久化与保存到磁盘
% 根据仿真背景生成唯一的输出文件名
% 如果配置中包含层间拓扑信息，则在文件名中包含该信息
if isfield(cfg, 'inter_topo') && ~isempty(cfg.inter_topo)
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_inter_%s_results.mat', ...
        lower(topo_type), cfg.N, cfg.K, topo_val, cfg.alpha, cfg.beta, cfg.noise, cfg.inter_topo);
else
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        lower(topo_type), cfg.N, cfg.K, topo_val, cfg.alpha, cfg.beta, cfg.noise);
end
% 指定到 results/ 目录下保存
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
out_path = fullfile(out_dir, mat_name);

% 将所有的分析结果存入 .mat 文件供后续绘图脚本使用
save(out_path, 'sigma_range', 'A_fwd', 'A_bwd');
fprintf('[DONE] 滞后环扫描已完成，总耗时: %.2f 秒。\n存入: %s\n', toc(simStart), out_path);

end
