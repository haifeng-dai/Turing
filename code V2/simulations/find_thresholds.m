function thresholds = find_thresholds(topo_type, topo_val, cfg)
% FIND_THRESHOLDS: 自动探测系统的正向（sigma_f）和反向（sigma_b）临界门槛点

%% 1. 环境与预热配置
rng(0, 'twister');        % 锁定随机种子，确保临界点探测结果的可重现性
topo_base = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'topology');
type_u = upper(topo_type); % 将拓扑名转为大写

%% 2. 加载底层网络拓扑库 (优化：若已经注入则跳过加载)
if ~isfield(cfg, 'L_intra') || isempty(cfg.L_intra)
    switch type_u
        case 'ER', net_file = fullfile(topo_base, 'ER', sprintf('N%d_p%.3f.mat', cfg.N, topo_val));
        case 'WS', net_file = fullfile(topo_base, 'WS', sprintf('N%d_pr%.2f.mat', cfg.N, topo_val));
        case 'BA', net_file = fullfile(topo_base, 'BA', sprintf('N%d_m%d.mat', cfg.N, topo_val));
        case 'SF', net_file = fullfile(topo_base, 'SF', sprintf('N%d_g%.1f.mat', cfg.N, topo_val));
        otherwise, error('未知拓扑类型: %s', topo_type);
    end

    if ~exist(net_file, 'file'), error('找不到网络文件: %s', net_file); end
    data = load(net_file, 'nets');
    cfg.L_intra = data.nets;
end

%% 3. 仿真实时控制参数 (针对阈值搜索的优化)
cfg.dt = 0.001;                 % 使用高效的大步长，加快逼近速度
cfg.steps = 2;                  % 二分搜索只关心终态稳定性，不记录过程数据

% 设置探测准则
thresh_A = 0.05;                % 序参数 A 高于此值视为“斑图存在”，低于此值视为“均一态”
epsilon  = 0.05;                % 二分法的搜索精度，当区间小于此值时停止探测

%% 4. 反向扫描用的“强斑图种子”获取 (优化：若已经注入则直接使用)
if isfield(cfg, 'y_seed') && ~isempty(cfg.y_seed)
    y_seed = cfg.y_seed;
else
    % --- 【修订】强制从磁盘读取全局唯一 100 种子以保证全项目一致性 (拒绝实时训练回退) ---
    seed_file = fullfile(fileparts(mfilename('fullpath')), '..', 'results', ...
        sprintf('evolution_%s_N200_K5_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', lower(topo_type)));

    if exist(seed_file, 'file')
        fprintf('[PRE] 正在从磁盘加载全局强斑图种子 (sigma=100)...\n');
        tmp = load(seed_file, 'Y');
        y_seed = tmp.Y(end, :)';
    else
        % --- 拒绝触发 solve_multiplex 生成随机种子，强制要求手动生成 ---
        error('【严重错误】未找到全局种子: %s\n请先以 USE_SEED=false, sigma=100 运行以生成标准斑图库。', seed_file);
    end
end

%% 5. 第一阶段：二分查找探测 sigma_f (正向扫描分叉点)
% 目标：从噪声启动，寻找能自发产生斑图的最小 sigma
low = cfg.sigma_min; high = cfg.sigma_max;
% fprintf('\n[STEP 1] 正在通过二分法探测正向临界点 (sigma_f)...\n');

while (high - low) > epsilon
    mid = (low + high) / 2;     % 取区间中点进行测试
    c = cfg; c.sigma = mid;
    c.y0 = [];                  % 核心：设置初值为空（触发随机噪声扰动）
    c.early_stop = 'forward';    % 启动早停优化 (一旦 A 达标则提前结束仿真)

    [~, ~, c_out] = solve_multiplex(c); % 执行反应扩散仿真
    A_final = c_out.A_final;

    if A_final > thresh_A       % 如果在中点处形成了斑图
        high = mid;             % 临界点在下方，收缩 high
        % fprintf('  > sigma=%.2f: %.2f能长出图 [收缩向左]\n', mid, A_final);
    else                        % 如果仍然是平庸均一态
        low = mid;              % 临界点在上方，收缩 low
        % fprintf('  > sigma=%.2f: %.2f还是无图 [收缩向右]\n', mid, A_final);
    end
end
sigma_f = (low + high) / 2;     % 获取最终探测到的正向临界值

%% 6. 第二阶段：二分查找探测 sigma_b (反向扫描坍塌点)
% 目标：从斑图种子启动，寻找能维持斑图活力的最小 sigma
low = cfg.sigma_min; high = cfg.sigma_max;
% fprintf('\n[STEP 2] 正在通过二分法探测反向临界点 (sigma_b)...\n');

while (high - low) > epsilon
    mid = (low + high) / 2;
    c = cfg; c.sigma = mid;
    c.y0 = y_seed;              % 核心：设置初值为预训练好的“种子”
    c.early_stop = 'backward';  % 启动早停优化 (一旦斑图崩溃提前退出)

    [~, ~, c_out] = solve_multiplex(c);
    A_final = c_out.A_final;

    if A_final > thresh_A       % 如果斑图在此处顽强存活
        high = mid;             % 继续访问更低的临界极限
        % fprintf('  > sigma=%.2f: %.2f 斑图存活 [向左探测]\n', mid, A_final);
    else                        % 如果斑图彻底崩坏消失了
        low = mid;              % 坍塌点在右侧，收缩 low
        % fprintf('  > sigma=%.2f: %.2f 斑图崩坏 [向右回退]\n', mid, A_final);
    end
end
sigma_b = (low + high) / 2;     % 获取最终探测到的反向临界值
thresholds = [sigma_f, sigma_b];
end
