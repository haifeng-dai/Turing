function pattern_evolution(DYNA, TOPO_TYPE, TOPO_PARAM, USE_SEED)
%% pattern_evolution.m: 图样演化保存版 (中英混合对齐)
% clc; close all;
addpath('simulations', 'networks');

% USE_SEED 由调用方指定: true 从斑图种子开始反向模拟，false 从随机初值开始正向模拟。

%% 2. 文件路径定义 & 仿真执行区
% 自动映射分支扫描方向
if USE_SEED, scan_mode = 'bwd'; else, scan_mode = 'fwd'; end

mat_name = sprintf('evolution_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_s%.1f_n%.2f_%s_results.mat', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, DYNA.beta, DYNA.sigma, DYNA.noise, scan_mode);
out_path = fullfile(fileparts(mfilename('fullpath')), 'results', mat_name);

% --- A. 动态加载层内拓扑库 (适配扁平化映射) ---
topo_base = fullfile(fileparts(mfilename('fullpath')), 'results', 'topology');
type_u = upper(TOPO_TYPE);
switch type_u
    case 'ER', net_file = fullfile(topo_base, 'ER', sprintf('N%d_p%.3f.mat', DYNA.N, TOPO_PARAM));
    case 'WS', net_file = fullfile(topo_base, 'WS', sprintf('N%d_pr%.2f.mat', DYNA.N, TOPO_PARAM));
    case 'BA', net_file = fullfile(topo_base, 'BA', sprintf('N%d_m%d.mat', DYNA.N, TOPO_PARAM));
    case 'SF', net_file = fullfile(topo_base, 'SF', sprintf('N%d_g%.1f.mat', DYNA.N, TOPO_PARAM));
    otherwise, error('未知拓扑类型: %s', TOPO_TYPE);
end

if ~exist(net_file, 'file'), error('找不到网络文件: %s', net_file); end
data = load(net_file, 'nets');
DYNA.L_intra = data.nets; % 此处正式填入 L_intra 字段

% --- B. 构造层间拉普拉斯算子 L_inter ---
adj_inter = ones(DYNA.K) - eye(DYNA.K);
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

if USE_SEED
    fprintf('[PRE] 正在加载全局强斑图种子 (sigma=100) 作为演化初值 (反向模式)...\n');
    % 使用与 sweep_param.m 相同的标准种子路径逻辑
    seed_name = sprintf('evolution_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_s100.0_n0.00_fwd_results.mat', ...
        lower(TOPO_TYPE), DYNA.N, DYNA.K, 0.03, 0.05, 0.005);
    seed_file = fullfile(fileparts(mfilename('fullpath')), 'results', seed_name);

    if ~exist(seed_file, 'file')
        error('找不到标准种子文件: %s\n请先设置 USE_SEED=false, sigma=100.0 运行以生成种子。', seed_file);
    end
    tmp_seed = load(seed_file);
    DYNA.y0 = tmp_seed.Y(end, :)';
else
    fprintf('[PRE] 从随机噪声开始模拟 (正向模式)...\n');
    DYNA.y0 = []; % 触发 solve_multiplex 内部的随机初值生成
end

%% 4. 执行单点演化仿真 (主任务)
fprintf('[SIM] 正在执行单点演化仿真 (sigma=%.2f)...\n', DYNA.sigma);
[t, Y, cfg] = solve_multiplex(DYNA); % 核心解算

% 根据第一层节点的度进行排序 (用于后续轨迹分析)
degs = -diag(DYNA.L_intra{1});
[~, sort_idx] = sort(degs, 'descend');

% --- 结果持久化存档 ---
save(out_path, 't', 'Y', 'cfg', 'sort_idx');
fprintf('[DONE] 仿真分析已完成并存入: %s\n', mat_name);

end
