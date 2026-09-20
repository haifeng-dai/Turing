%% sigma_beta_noise_er.m: 层间耦合(beta)与噪声(eta)对滞回阈值的影响分析
clear; clc; close all;
addpath('simulations', 'networks');

RUN_SIMULATION = false;  % 控制是否运行仿真 (true: 运行并保存, false: 直接读取并绘图)

%% 1. 实验参数配置
TOPO_TYPE = 'ER';
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
P_VAL_FIXED = 0.03;      % 固定连接概率，研究噪声影响
DYNA.T_END = 500;        % 临界探测不需要太长平衡时间
DYNA.steps = 2;          % 扫描任务仅取两张快照
DYNA.init_perturb = 0.1;

% 搜索范围设置
DYNA.sigma_min = 0;
DYNA.sigma_max = 80;

% 扫描变量定义
ETA_LIST     = 0:0.1:1;
BETA_FACTORS = 0:0.1:10;
BETA_LIST    = BETA_FACTORS * DYNA.alpha;

% 数据保存路径
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end
data_path = fullfile(res_dir, sprintf('scan_beta_noise_%s_N%d_K%d_a%.3f_p%.3f.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, P_VAL_FIXED));

%% 2. 核心并行计算区
num_e = length(ETA_LIST);
num_b = length(BETA_LIST);

if RUN_SIMULATION
    % === 并行池初始化：充分利用所有64个核 ===
    if isempty(gcp('nocreate'))  % 检查是否已有并行池
        parpool('Threads', 64);  % 创建64个线程的并行池
        fprintf('[POOL] 已初始化 64 线程并行池\n');
    else
        pool = gcp;
        fprintf('[POOL] 使用现有并行池 (%d 个worker)\n', pool.NumWorkers);
    end

    % *** 精细扫描已经在后面的 parfor 中充分并行化 ***
    [B_idx, E_idx] = meshgrid(1:num_b, 1:num_e);
    B_idx = B_idx(:);
    E_idx = E_idx(:);
    total_tasks = length(B_idx);  % = num_e * num_b = 11 * 101 = 1111 个任务

    Sigma_F_Results = zeros(total_tasks, 1);
    Sigma_B_Results = zeros(total_tasks, 1);

    q = parallel.pool.DataQueue;
    hWait = waitbar(0, '正在初始化...', 'Name', '并行任务进度');
    cleanupObj = onCleanup(@() delete(hWait)); % 确保代码结束或报错时关闭进度条
    tStart_sim = tic;
    afterEach(q, @(~) nUpdateProgress(total_tasks, hWait, tStart_sim));

    % --- A. 构造全局层间拉普拉斯算子 ---
    adj_inter = (ones(DYNA.K) - eye(DYNA.K));
    DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

    % --- B. 加载全局强斑图种子 (保持项目一致性) ---
    fprintf('[PRE] 正在加载全局强斑图种子 (sigma=100)...\n');
    seed_file = fullfile('results', 'evolution_er_N200_K5_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat');
    if ~exist(seed_file, 'file'), error('找不到全局种子文件，请先运行 pattern_evolution.m！'); end
    tmp_seed = load(seed_file);
    y_universal_seed = tmp_seed.Y(end, :)';

    % --- C. 加载固定拓扑结构 ---
    fprintf('[PRE] 正在加载固定拓扑 (p=%.3f)...\n', P_VAL_FIXED);
    net_path = fullfile('results', 'topology', 'ER', sprintf('N%d_p%.3f.mat', DYNA.N, P_VAL_FIXED));
    tmp_net = load(net_path);
    L_intra_lib = tmp_net.nets;

    % 1. 粗略扫描 (Rough Scan) - 快速参考版本（8个采样点）
    fprintf('\n[SIM] 执行粗略扫描以提取预测线 (Rough Scan)...[快速参考版本]\n');
    n_rough = 8; % 仅8个点快速获得趋势，粗扫描占总计算~1/12（88 vs 1111）
    idx_rough = round(linspace(1, length(BETA_LIST), n_rough));
    beta_rough = BETA_LIST(idx_rough);

    sf_rough = zeros(num_e, n_rough);
    sb_rough = zeros(num_e, n_rough);

    % 平坦化索引：(e_i, b_i) -> 线性索引，充分利用64个核
    [B_rough, E_rough] = meshgrid(1:n_rough, 1:num_e);  % 创建 (eta, beta) 网格
    B_rough = B_rough(:);   % 扁平化为列向量
    E_rough = E_rough(:);
    n_rough_tasks = length(E_rough);

    sf_rough_vec = zeros(n_rough_tasks, 1);
    sb_rough_vec = zeros(n_rough_tasks, 1);

    % 完全并行：num_e * n_rough 个任务在64核上并行执行
    parfor task_id = 1:n_rough_tasks
        c = DYNA;
        c.noise = ETA_LIST(E_rough(task_id));
        c.beta = beta_rough(B_rough(task_id));
        c.L_intra = L_intra_lib;
        c.y_seed = y_universal_seed;
        res = find_thresholds(TOPO_TYPE, P_VAL_FIXED, c);
        sf_rough_vec(task_id) = res(1);
        sb_rough_vec(task_id) = res(2);
    end

    % 重塑回矩阵形式
    sf_rough = reshape(sf_rough_vec, [num_e, n_rough]);
    sb_rough = reshape(sb_rough_vec, [num_e, n_rough]);

    % 2. 预测插值 (Prediction)
    sigma_pred_f = zeros(num_e, num_b);
    sigma_pred_b = zeros(num_e, num_b);
    for e_i = 1:num_e
        sigma_pred_f(e_i, :) = interp1(beta_rough, sf_rough(e_i, :), BETA_LIST, 'pchip', DYNA.sigma_max);
        sigma_pred_b(e_i, :) = interp1(beta_rough, sb_rough(e_i, :), BETA_LIST, 'pchip', DYNA.sigma_max);
    end

    % --- D. 批量执行二分查找探测 ---
    fprintf('\n[SIM] 开始局域精细并行扫描 (Guided Fine Scan): %d 组任务...[64核充分并行]\n', total_tasks);
    simStart = tic;

    parfor t = 1:total_tasks
        dyna_local = DYNA;
        dyna_local.beta  = BETA_LIST(B_idx(t));
        dyna_local.noise = ETA_LIST(E_idx(t));

        % 注入网络和种子
        dyna_local.L_intra = L_intra_lib;
        dyna_local.y_seed  = y_universal_seed;

        % 利用预测值缩小搜索区间 (+- 10 留出裕量)
        pred_center = max(sigma_pred_f(E_idx(t), B_idx(t)), sigma_pred_b(E_idx(t), B_idx(t)));
        if pred_center > DYNA.sigma_max - 5
            dyna_local.sigma_min = 0;
            dyna_local.sigma_max = DYNA.sigma_max;
        else
            dyna_local.sigma_min = max(0, pred_center - 10);
            dyna_local.sigma_max = min(DYNA.sigma_max, pred_center + 10);
        end

        res = find_thresholds(TOPO_TYPE, P_VAL_FIXED, dyna_local);
        Sigma_F_Results(t) = res(1);
        Sigma_B_Results(t) = res(2);

        % 发送信号以更新进度
        send(q, t);
    end

    simTime = toc(simStart);
    fprintf('[DONE] 仿真完成！总耗时: %.2f 秒。\n', simTime);

    % 整理并持久化保存
    SF_Matrix = reshape(Sigma_F_Results, [num_e, num_b]);
    SB_Matrix = reshape(Sigma_B_Results, [num_e, num_b]);
    save(data_path, 'SF_Matrix', 'SB_Matrix', 'ETA_LIST', 'BETA_LIST', 'P_VAL_FIXED');
end

%% 3. 数据可视化
fprintf('\n[VIS] 准备数据可视化...\n');
if exist(data_path, 'file')
    load(data_path);
    fprintf('[LOAD] 已加载结果文件: %s\n', data_path);
else
    fprintf('[WARN] 结果文件不存在！请先运行仿真 (RUN_SIMULATION=true)。\n');
    num_e = length(ETA_LIST);
    num_b = length(BETA_LIST);
    return;
end

% 重新定义num_e和num_b (从加载的列表长度获取)
num_e = length(ETA_LIST);
num_b = length(BETA_LIST);

h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.45, 0.5]);
ax = axes('Position', [0.15, 0.15, 0.75, 0.75]); hold on;

colors = jet(num_e); % 使用渐变色区分噪声强度
for i = 1:num_e
    % 正向扫描 (实线)
    plot(BETA_LIST, SF_Matrix(i, :), '-', 'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('\\eta=%.1f (\\sigma_f)', ETA_LIST(i)));
    % 反向扫描 (虚线)
    plot(BETA_LIST, SB_Matrix(i, :), '--', 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'HandleVisibility', 'off'); % 隐藏虚线的 legend 以免太乱
end

% 属性精修
xlabel('Inter-layer Coupling Strength \beta', 'FontSize', 14);
ylabel('Critical Diffusion Ratio \sigma', 'FontSize', 14);
title(sprintf('Impact of Noise (ER, p=%.3f)', P_VAL_FIXED), 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 12);
grid on; box on;
set(ax, 'FontSize', 14);

plots_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
out_img = fullfile(plots_dir, sprintf('scan_beta_noise_%s_p%.3f.png', lower(TOPO_TYPE), P_VAL_FIXED));
exportgraphics(h, out_img, 'Resolution', 300);
fprintf('[DONE] 绘图已生成: %s\n', out_img);

% ========================================================
% 进度监控子函数 (基于 DataQueue)
% ========================================================
function nUpdateProgress(total, hWait, tStart)
persistent count
if isempty(count), count = 0; end
count = count + 1;

% 更新 waitbar
progress = count / total;
waitbar(progress, hWait, sprintf('处理进度: %d/%d (%.0f%%) | 耗时: %.1f 秒', ...
    count, total, progress*100, toc(tStart)));
end
