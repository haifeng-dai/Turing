function sigma_eta_connectivity_er(DYNA, TOPO_TYPE, P_LIST, ETA_LIST)
% SIGMA_ETA_CONNECTIVITY_ER: 计算并保存噪声与连接概率对应的临界扩散比。
function_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(function_dir, 'simulations'), fullfile(function_dir, 'networks'));

% 数据保存路径 (动态包含 beta 强度以防覆盖)
res_dir = fullfile(function_dir, 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end
data_name = sprintf('scan_eta_p_%s_N%d_K%d_a%.3f_b%.3f.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta);
if isfield(DYNA, 'noise_seed') && ~isempty(DYNA.noise_seed)
    [~, data_stem, data_ext] = fileparts(data_name);
    data_name = sprintf('%s_seed%.0f%s', data_stem, DYNA.noise_seed, data_ext);
end
data_path = fullfile(res_dir, data_name);

%% 2. 核心并行计算区
num_p = length(P_LIST);
num_e = length(ETA_LIST);
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
seed_file = fullfile(res_dir, sprintf( ...
    'evolution_er_N%d_K%d_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', DYNA.N, DYNA.K));
if ~exist(seed_file, 'file')
    seed_dyna = DYNA;
    seed_dyna.sigma = 100;
    seed_dyna.noise = 0;
    pattern_evolution(seed_dyna, TOPO_TYPE, 0.03, false);
end
if ~exist(seed_file, 'file'), error('标准种子生成后仍未找到文件: %s', seed_file); end
tmp_seed = load(seed_file);
y_universal_seed = tmp_seed.Y(end, :)';

fprintf('[PRE] 正在并行预加载 %d 组拓扑...\n', num_p);
NET_LIBS = cell(num_p, 1);
parfor p_i = 1:num_p
    p_val = P_LIST(p_i);
    net_path = fullfile(res_dir, 'topology', 'ER', sprintf('N%d_p%.3f.mat', DYNA.N, p_val));
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
