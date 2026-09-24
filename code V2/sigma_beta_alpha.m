function sigma_beta_alpha(DYNA, TOPO_TYPE, P_VAL_FIXED, ALPHA_LIST, RATIO_LIST, RUN_SIMULATION)
%% sigma_beta_alpha.m: 层间耦合(beta)与层内扩散(alpha)对滞回阈值的影响分析
addpath('simulations', 'networks');

% 数据保存路径
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end
data_path = fullfile(res_dir, sprintf('scan_beta_alpha_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, P_VAL_FIXED));

%% 2. 核心并行计算区 [64核优化版]
num_a = length(ALPHA_LIST);
num_r = length(RATIO_LIST);

if RUN_SIMULATION
    % === 并行池初始化：充分利用所有64个核 ===
    if isempty(gcp('nocreate'))  % 检查是否已有并行池
        parpool('Threads', 64);  % 创建64个线程的并行池
        fprintf('[POOL] 已初始化 64 线程并行池\n');
    else
        pool = gcp;
        fprintf('[POOL] 使用现有并行池 (%d 个worker)\n', pool.NumWorkers);
    end

    % *** 精细扫描已经在后面的 parfor 中充分并行化，以下索引用于主循环 ***
    [R_idx, A_idx] = meshgrid(1:num_r, 1:num_a);
    R_idx = R_idx(:);
    A_idx = A_idx(:);
    total_tasks = length(R_idx);  % = num_a * num_r = 16 * 101 = 1616 个任务

    Sigma_F_Results = zeros(total_tasks, 1);
    Sigma_B_Results = zeros(total_tasks, 1);

    % --- A. 进度监控器初始化 (DataQueue 方案) ---
    q = parallel.pool.DataQueue;
    hWait = waitbar(0, '正在初始化...', 'Name', '并行任务进度');
    cleanupObj = onCleanup(@() delete(hWait)); % 确保代码结束或报错时关闭进度条
    tStart_sim = tic;
    afterEach(q, @(~) nUpdateProgress(total_tasks, hWait, tStart_sim));

    % --- B. 构造层间拉普拉斯算子 ---
    adj_inter = (ones(DYNA.K) - eye(DYNA.K));
    DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

    % --- C. 加载全局强制斑图种子 (保证实验一致性) ---
    fprintf('[PRE] 正在加载全局强斑图种子...\n');
    seed_file = fullfile('results', 'evolution_er_N200_K5_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat');
    if ~exist(seed_file, 'file'), error('找不到全局种子文件，请先运行 pattern_evolution.m！'); end
    tmp_seed = load(seed_file);
    y_universal_seed = tmp_seed.Y(end, :)';

    % --- D. 加载固定拓扑网络 ---
    net_path = fullfile('results', 'topology', 'ER', sprintf('N%d_p%.3f.mat', DYNA.N, P_VAL_FIXED));
    tmp_net = load(net_path);
    L_intra_lib = tmp_net.nets;

    % 1. 粗略扫描 (Rough Scan) - 快速参考版本（8个采样点）
    fprintf('\n[SIM] 执行粗略扫描以提取预测线 (Rough Scan)...[快速参考版本]\n');
    n_rough = 8; % 仅8个点快速获得趋势，粗扫描占总计算~1/12（128 vs 1616）
    idx_rough = round(linspace(1, length(RATIO_LIST), n_rough));
    ratio_rough = RATIO_LIST(idx_rough);

    % 平坦化索引：(a_i, r_i) -> 线性索引，充分利用64个核
    [R_rough, A_rough] = meshgrid(1:n_rough, 1:num_a);  % 创建 (alpha, ratio) 网格
    R_rough = R_rough(:);   % 扁平化为列向量
    A_rough = A_rough(:);
    n_rough_tasks = length(A_rough);

    sf_rough_vec = zeros(n_rough_tasks, 1);
    sb_rough_vec = zeros(n_rough_tasks, 1);

    % 完全并行：num_a * n_rough 个任务在64核上并行执行
    parfor task_id = 1:n_rough_tasks
        c = DYNA;
        c.alpha = ALPHA_LIST(A_rough(task_id));
        c.beta = c.alpha * ratio_rough(R_rough(task_id));
        c.L_intra = L_intra_lib;
        c.y_seed = y_universal_seed;
        res = find_thresholds(TOPO_TYPE, P_VAL_FIXED, c);
        sf_rough_vec(task_id) = res(1);
        sb_rough_vec(task_id) = res(2);
    end

    % 重塑回矩阵形式
    sf_rough = reshape(sf_rough_vec, [num_a, n_rough]);
    sb_rough = reshape(sb_rough_vec, [num_a, n_rough]);

    % 2. 预测插值 (Prediction)
    sigma_pred_f = zeros(num_a, num_r);
    sigma_pred_b = zeros(num_a, num_r);
    for a_i = 1:num_a
        sigma_pred_f(a_i, :) = interp1(ratio_rough, sf_rough(a_i, :), RATIO_LIST, 'pchip', DYNA.sigma_max);
        sigma_pred_b(a_i, :) = interp1(ratio_rough, sb_rough(a_i, :), RATIO_LIST, 'pchip', DYNA.sigma_max);
    end

    fprintf('\n[SIM] 开始局域精细并行扫描 (Guided Fine Scan): %d 组任务... \n', total_tasks);
    simStart = tic;

    parfor t = 1:total_tasks
        dyna_local = DYNA;
        cur_alpha = ALPHA_LIST(A_idx(t));
        dyna_local.alpha = cur_alpha;
        dyna_local.beta  = cur_alpha * RATIO_LIST(R_idx(t)); % 使用倍数关系

        % 注入网络和种子
        dyna_local.L_intra = L_intra_lib;
        dyna_local.y_seed  = y_universal_seed;

        % 利用预测值缩小搜索区间 (+- 10 留出裕量)
        pred_center = max(sigma_pred_f(A_idx(t), R_idx(t)), sigma_pred_b(A_idx(t), R_idx(t)));
        if pred_center > DYNA.sigma_max - 5
            dyna_local.sigma_min = 0;
            dyna_local.sigma_max = DYNA.sigma_max;
        else
            dyna_local.sigma_min = max(0, pred_center - 10);
            dyna_local.sigma_max = min(DYNA.sigma_max, pred_center + 10);
        end

        % 调用二分搜索探测阈值
        res = find_thresholds(TOPO_TYPE, P_VAL_FIXED, dyna_local);
        Sigma_F_Results(t) = res(1);
        Sigma_B_Results(t) = res(2);

        % 发送信号以更新进度
        send(q, t);
    end

    simTime = toc(simStart);
    fprintf('[DONE] 仿真完成！总耗时: %.2f 秒。\n', simTime);

    % 整理并持久化保存
    SF_Matrix = reshape(Sigma_F_Results, [num_a, num_r]);
    SB_Matrix = reshape(Sigma_B_Results, [num_a, num_r]);
    save(data_path, 'SF_Matrix', 'SB_Matrix', 'ALPHA_LIST', 'RATIO_LIST', 'P_VAL_FIXED');
    fprintf('[SAVE] 结果已保存到: %s\n', data_path);
end

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

end
