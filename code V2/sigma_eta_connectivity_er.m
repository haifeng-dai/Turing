function sigma_eta_connectivity_er(DYNA, TOPO_TYPE, P_VAL, ETA_LIST)
% SIGMA_ETA_CONNECTIVITY_ER: 计算并保存特定单个 p 下，噪声 eta 对应的临界扩散比。
function_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(function_dir, 'simulations'), fullfile(function_dir, 'networks'));

% 1. 数据保存路径：文件名中显式包含 p 值
res_dir = fullfile(function_dir, 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end

data_name = sprintf('scan_eta_p%.3f_%s_N%d_K%d_a%.3f_b%.3f.mat', ...
    P_VAL, lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta);

if isfield(DYNA, 'noise_seed') && ~isempty(DYNA.noise_seed)
    [~, data_stem, data_ext] = fileparts(data_name);
    data_name = sprintf('%s_seed%.0f%s', data_stem, DYNA.noise_seed, data_ext);
end
data_path = fullfile(res_dir, data_name);

%% 2. 环境与前置加载
num_e = length(ETA_LIST);
Sigma_F_Results = zeros(1, num_e);
Sigma_B_Results = zeros(1, num_e);

% A. 构造全局层间算子
adj_inter = (ones(DYNA.K) - eye(DYNA.K));
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

% B. 加载全局标准强斑图种子 (sigma=100)
seed_file = fullfile(res_dir, sprintf( ...
    'evolution_er_N%d_K%d_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', DYNA.N, DYNA.K));
if ~exist(seed_file, 'file')
    seed_dyna = DYNA;
    seed_dyna.sigma = 100;
    seed_dyna.noise = 0;
    pattern_evolution(seed_dyna, TOPO_TYPE, 0.03, false);
end
tmp_seed = load(seed_file);
y_universal_seed = tmp_seed.Y(end, :)';

% C. 预加载当前单个 p 对应的拓扑网络
net_path = fullfile(res_dir, 'topology', 'ER', sprintf('N%d_p%.3f.mat', DYNA.N, P_VAL));
if ~exist(net_path, 'file'), error('找不到网络拓扑文件: %s', net_path); end
tmp = load(net_path);
L_intra_current = tmp.nets;

%% 3. 并行扫描 (仅针对当前的 P_VAL 扫描所有 eta)
fprintf('\n[SIM] 开始并行扫描 (p=%.3f, 共 %d 个 eta 任务)... \n', P_VAL, num_e);
simStart = tic;

parfor e_i = 1:num_e
    loopStart = tic;
    dyna_local = DYNA;
    dyna_local.noise = ETA_LIST(e_i);

    % 注入网络与种子
    dyna_local.L_intra = L_intra_current;
    dyna_local.y_seed  = y_universal_seed;

    res = find_thresholds(TOPO_TYPE, P_VAL, dyna_local);
    Sigma_F_Results(e_i) = res(1);
    Sigma_B_Results(e_i) = res(2);

    fprintf('  [PAR] 任务 %d/%d: p=%.3f, eta=%.3f 完成，耗时: %.2f 秒。\n', ...
        e_i, num_e, P_VAL, dyna_local.noise, toc(loopStart));
end

simTime = toc(simStart);
fprintf('[DONE] p=%.3f 扫描完成！总耗时: %.2f 秒。\n', P_VAL, simTime);

% 持久化保存单个 p 的结果
sf_vec = Sigma_F_Results;
sb_vec = Sigma_B_Results;
p_val  = P_VAL;
save(data_path, 'p_val', 'sf_vec', 'sb_vec', 'ETA_LIST');
fprintf('[SAVE] 结果已成功保存到: %s\n', data_path);

end
