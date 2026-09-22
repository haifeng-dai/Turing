%% generate_standard_networks.m: 生成标准多层网络拓扑库 (ER/WS/BA/SF)
% 更新：已采用扁平化输出结构，每组参数仅对应一个包含了 20 层储备的 .mat 文件。
% 存储样例: results/topology/ER/N200_p0.030.mat

clear; clc;
addpath('networks');

% ==========================================
% 1. 全局配置 (单一实例模式)
% ==========================================
N_vals = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000];      % 节点规模
K = 10;                   % 每套实例内部包含 20 层储备 (方便按需 K 计算)

% --- 2. 拓扑参数扫频网格 ---
p_vals = 0.01:0.01:1.0;            % ER 连通概率
k_vals_ws = 4:2:20;                % WS 平均度 K (4,6,8,10,12,14,16,18,20)
p_rewire_vals = 0.01:0.01:1.0;     % WS 重连概率
% m_vals = 2:1:10;                   % BA 新增边数
m_vals = [20];                   % BA 新增边数
gamma_vals = 2.1:0.1:3.0;          % SF 幂律指数

fprintf('[INIT] 开始生成标准拓扑库 (扁平化存储, K=%d)...\n', K);
overallStart = tic;

script_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(script_dir, 'results', 'topology');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

% ==========================================
% 3. 嵌套生成逻辑 (并行化参数层级)
% ==========================================
for n_idx = 1:length(N_vals)
    N = N_vals(n_idx);
    fprintf('\n[SCALE] 正在处理规模 N = %d ...\n', N);

    % --- 3.1 ER 网络 (扁平存放在 ER 根目录) ---
    fprintf('  > 正在并行生成 ER 拓扑库...\n');
    er_root = fullfile(results_dir, 'ER');
    if ~exist(er_root, 'dir'), mkdir(er_root); end

    % parfor i = 1:length(p_vals)
    %     p_val = p_vals(i);
    %     save_name = fullfile(er_root, sprintf('N%d_p%.3f.mat', N, p_val));

    %     rng(N + i * 1000, 'twister');
    %     nets = cell(K, 1);
    %     for k = 1:K
    %         nets{k} = sparse(gen_er(N, p_val));
    %     end
    %     save_parfor(save_name, nets);
    % end

    % --- 3.2 WS 随机重连网络 (两层参数：K 和 p_rewire) ---
    fprintf('  > 正在并行生成 WS 拓扑库 (K=%d种 × p_rewire=%d种)...\n', length(k_vals_ws), length(p_rewire_vals));
    ws_root = fullfile(results_dir, 'WS');
    if ~exist(ws_root, 'dir'), mkdir(ws_root); end

    % 为避免 parfor 嵌套，先构造线性索引
    ws_tasks = [];
    for k_idx = 1:length(k_vals_ws)
        for pr_idx = 1:length(p_rewire_vals)
            ws_tasks = [ws_tasks; k_idx, pr_idx];
        end
    end
    n_ws_tasks = size(ws_tasks, 1);

    % parfor i = 1:n_ws_tasks
    %     k_idx = ws_tasks(i, 1);
    %     pr_idx = ws_tasks(i, 2);
    %     k_val = k_vals_ws(k_idx);
    %     pr = p_rewire_vals(pr_idx);
    %     save_name = fullfile(ws_root, sprintf('N%d_K%d_pr%.2f.mat', N, k_val, pr));

    %     rng(N + k_idx * 100 + pr_idx * 1000 + 10000, 'twister');
    %     nets = cell(K, 1);
    %     for k = 1:K
    %         nets{k} = sparse(gen_ws(N, k_val, pr));
    %     end
    %     save_parfor(save_name, nets);
    % end

    % --- 3.3 BA 无标度网络 ---
    fprintf('  > 正在并行生成 BA 拓扑库...\n');
    ba_root = fullfile(results_dir, 'BA');
    if ~exist(ba_root, 'dir'), mkdir(ba_root); end

    parfor i = 1:length(m_vals)
        m = m_vals(i);
        save_name = fullfile(ba_root, sprintf('N%d_m%d.mat', N, m));

        rng(N + i * 1000 + 20000, 'twister');
        nets = cell(K, 1);
        for k = 1:K
            nets{k} = sparse(gen_ba(N, m));
        end
        save_parfor(save_name, nets);
    end

    % --- 3.4 SF 幂律分布网络 ---
    fprintf('  > 正在并行生成 SF 拓扑库...\n');
    sf_root = fullfile(results_dir, 'SF');
    if ~exist(sf_root, 'dir'), mkdir(sf_root); end

    % parfor i = 1:length(gamma_vals)
    %     g = gamma_vals(i);
    %     save_name = fullfile(sf_root, sprintf('N%d_g%.1f.mat', N, g));

    %     rng(N + i * 1000 + 30000, 'twister');
    %     nets = cell(K, 1);
    %     for k = 1:K
    %         nets{k} = sparse(gen_sf_gamma(N, g));
    %     end
    %     save_parfor(save_name, nets);
    % end
end

fprintf('\n[DONE] 全套拓扑库已重构生成完毕！总耗时: %.1f 秒\n', toc(overallStart));

%% 辅助函数 (parfor 专用)
function save_parfor(path, nets)
save(path, 'nets');
end
