%% 固定五层拓扑的图 2c：谱分析 + 五次噪声重复的实际阈值
% 直接运行：读取匹配的缓存，仅补算缺失任务，然后绘制四张图。
% 只画图：设 RUN_SIMULATION=false；缓存不完整时明确报错，不启动扫描。
% 噪声种子不生成网络，也不改变初始扰动；误差带为均值 +/- 1 样本标准差。
% 沿用 find_thresholds 的二分法、A>0.05 判据和 0.05 搜索精度。
% 所有配置和 r 共享项目预生成 ER 斑图种子，不进行参考态训练。
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(script_dir, 'networks'), fullfile(script_dir, 'simulations'));

RUN_SIMULATION = true;
FORCE_RERUN = false;       % 主动重算当前实验；不会重新生成已保存的网络
REQUESTED_WORKERS = 64;
BATCH_SIZE = 64;           % 每批完成后保存；中断最多损失当前批
SPECTRUM_BATCH_SIZE = 64;

config.N = 200;
config.K = 5;
config.alpha = 0.10;
config.noise = 0.01;
config.T_END = 500;
config.dt = 0.001;
config.steps = 2;
config.init_perturb = 0.1;
config.detect_convergence = true;
config.sigma_min = 0;
config.sigma_max = 80;
config.ratios = 0:0.1:10;
config.noise_seeds = 1:5;
config.topology_seed = 20260930;
config.ws_degree = 6;
config.ws_rewire = 0.10;
config.ba_m = 3;            % N=200 时实际平均度 5.94
config.er_p = 0.030;        % N=200 时期望平均度 5.97；WS 为 6
config.schema_version = 2;

assert(config.K == 5, '本实验约定五层网络。');
assert(numel(config.ratios) >= 2 && all(diff(config.ratios) > 0) && ...
    all(config.ratios >= 0), 'ratios 必须为至少两个递增的非负数。');
assert(numel(config.noise_seeds) >= 2 && ...
    numel(unique(config.noise_seeds)) == numel(config.noise_seeds), '需至少两个不同噪声种子。');
assert(~FORCE_RERUN || RUN_SIMULATION, 'FORCE_RERUN=true 时必须允许运行仿真。');
[shared_seed, seed_file] = load_standard_pattern_seed(config, RUN_SIMULATION);

%% 1. 独立拓扑缓存：改变噪声种子或扫描范围不会改变网络
topology_fields = {'N','K','topology_seed','ws_degree','ws_rewire','ba_m','er_p'};
for i = 1:numel(topology_fields)
    topology_config.(topology_fields{i}) = config.(topology_fields{i});
end
topology_key = sprintf('N%d_K%d_ws%d_pr%.3f_ba%d_er%.4f_toposeed%d', ...
    config.N, config.K, config.ws_degree, config.ws_rewire, ...
    config.ba_m, config.er_p, config.topology_seed);
cache_dir = fullfile(script_dir, 'results', 'fig2c_topology_mean_std', topology_key);
topology_file = fullfile(cache_dir, 'topologies.mat');
if exist(topology_file, 'file')
    topology = load(topology_file);
    if ~isfield(topology, 'topology_config') || ...
            ~isequaln(topology.topology_config, topology_config)
        error('拓扑缓存配置不匹配：%s。请选用新的 topology_seed。', topology_file);
    end
    cases = topology.cases;
    L_inter = topology.L_inter;
    fprintf('[NETWORK] 读取固定网络：%s\n', topology_file);
elseif RUN_SIMULATION
    if ~exist(cache_dir, 'dir'), mkdir(cache_dir); end
    cases = fig2c_topology_networks(config);
    L_inter = sparse(config.K*eye(config.K)-ones(config.K));
    save(topology_file, 'cases', 'L_inter', 'topology_config');
    fprintf('[NETWORK] 已生成并保存固定网络：%s\n', topology_file);
else
    error('仅绘图模式缺少拓扑缓存：%s。请先设 RUN_SIMULATION=true。', topology_file);
end
validateNetworks(cases, L_inter, config);
for c = 1:3
    fprintf('  %s：%s；各层实际平均度 %s\n', cases(c).label, ...
        strjoin(cases(c).layer_types, '-'), mat2str(cases(c).mean_degrees, 4));
end

%% 2. 实验缓存：统一种子协议使用新文件，旧逐 r 训练结果不混用、不覆盖
source_signature.solver = fileread(which('solve_multiplex'));
source_signature.thresholds = fileread(which('find_thresholds'));
source_signature.spectrum = fileread(which('fig2c_topology_spectrum'));
source_signature.seed_loader = fileread(which('load_standard_pattern_seed'));
source_signature.shared_seed = shared_seed; % 文件内容改变时拒绝混用旧阈值
source_signature.protocol_version = 2;
experiment_tag = sprintf('a%.3f_eta%.3f', config.alpha, config.noise);
result_file = fullfile(cache_dir, ['scan_sharedseed_' experiment_tag '.mat']);
num_r = numel(config.ratios);
num_seeds = numel(config.noise_seeds);
if exist(result_file, 'file') && ~FORCE_RERUN
    state = load(result_file);
    if ~isfield(state, 'config') || ~isequaln(state.config, config) || ...
            ~isfield(state, 'source_signature') || ...
            ~isequaln(state.source_signature, source_signature)
        error(['实验缓存与当前参数/积分代码不一致，拒绝混用。\n%s\n' ...
            '需要重算时设 FORCE_RERUN=true；旧实验如需保留，请先另存该文件。'], result_file);
    end
    assert(isequal(size(state.sf_repeats), [3 num_seeds num_r]) && ...
        isequal(size(state.sb_repeats), [3 num_seeds num_r]), '缓存阈值维度错误。');
    fprintf('[CACHE] 已读取断点：%s\n', result_file);
else
    state.config = config;
    state.source_signature = source_signature;
    state.topology_config = topology_config;
    state.topology_file = topology_file;
    state.sf_repeats = nan(3, num_seeds, num_r);
    state.sb_repeats = nan(3, num_seeds, num_r);
    state.sigma_spec = nan(3, num_r);
    state.lambda_star = nan(3, num_r);
    state.shared_seed = shared_seed;
    state.seed_file = seed_file;
    state.spectrum_done = false(3, num_r);
    state.elapsed_seconds = 0;
    if RUN_SIMULATION, checkpoint(result_file, state); end
end
finished = isfinite(state.sf_repeats) & isfinite(state.sb_repeats);
fprintf('[PLAN] 阈值任务 = 3 拓扑 x %d 噪声种子 x %d 耦合比 = %d；已完成 %d。\n', ...
    num_seeds, num_r, numel(finished), nnz(finished));
fprintf('[PLAN] 每个任务分别二分查找正向、反向阈值，不是一次积分。\n');
if ~RUN_SIMULATION && (~all(finished(:)) || ~all(state.spectrum_done(:)))
    error('仅绘图模式：缓存尚不完整。设 RUN_SIMULATION=true 可继续补算，不会丢失已完成任务。');
end

%% 3. 仅计算确定性谱阈值；不积分、不生成新的斑图初值
run_clock = tic;
previous_elapsed = state.elapsed_seconds;
workers = 0;
if RUN_SIMULATION && (~all(state.spectrum_done(:)) || ~all(finished(:)))
    workers = setupPool(REQUESTED_WORKERS);
else
    fprintf('[CACHE] 无需扫描，不启动并行池。\n');
end
base_cfg = struct('N', config.N, 'K', config.K, 'alpha', config.alpha, ...
    'noise', config.noise, 'T_END', config.T_END, 'dt', config.dt, ...
    'steps', config.steps, 'init_perturb', config.init_perturb, ...
    'sigma_min', config.sigma_min, 'sigma_max', config.sigma_max, ...
    'detect_convergence', config.detect_convergence, 'L_inter', L_inter);
pending = find(~state.spectrum_done);
for first = 1:SPECTRUM_BATCH_SIZE:numel(pending)
    ids = pending(first:min(first+SPECTRUM_BATCH_SIZE-1, numel(pending)));
    fprintf('\n[SPECTRUM] 本批 %d 个谱分析任务（无积分），已缓存 %d/%d。\n', ...
        numel(ids), nnz(state.spectrum_done), numel(state.spectrum_done));
    batch = cell(1, numel(ids));
    parfor (j = 1:numel(ids), workers)
        [c, ri] = ind2sub([3 num_r], ids(j));
        batch{j} = computeSpectrum(cases(c), config.ratios(ri), base_cfg);
    end
    failures = {};
    for j = 1:numel(ids)
        [c, ri] = ind2sub([3 num_r], ids(j));
        item = batch{j};
        if isempty(item.error)
            state.sigma_spec(c,ri) = item.spectrum;
            state.lambda_star(c,ri) = item.lambda_star;
            state.spectrum_done(c,ri) = true;
        else
            failures{end+1} = sprintf('%s, r=%.3g：%s', cases(c).label, ...
                config.ratios(ri), item.error); %#ok<SAGROW>
        end
    end
    state.elapsed_seconds = previous_elapsed+toc(run_clock);
    checkpoint(result_file, state);
    fprintf('[SAVED] 谱分析已完成 %d/%d，累计 %.1f s。\n', ...
        nnz(state.spectrum_done), numel(state.spectrum_done), state.elapsed_seconds);
    stopOnFailures(failures, result_file);
end

%% 4. 拓扑 x 噪声种子 x 耦合比扁平并行；无需按 r 拆文件
finished = isfinite(state.sf_repeats) & isfinite(state.sb_repeats);
pending = find(~finished);
for first = 1:BATCH_SIZE:numel(pending)
    ids = pending(first:min(first+BATCH_SIZE-1, numel(pending)));
    fprintf('\n[SCAN] 本批 %d 个正/反阈值任务，已缓存 %d/%d。\n', ...
        numel(ids), nnz(finished), numel(finished));
    results = nan(numel(ids), 2);
    errors = cell(1, numel(ids));
    parfor (j = 1:numel(ids), workers)
        [c, si, ri] = ind2sub([3 num_seeds num_r], ids(j));
        cfg = base_cfg;
        cfg.beta = config.alpha*config.ratios(ri);
        cfg.L_intra = cases(c).L;
        cfg.y_seed = shared_seed;
        cfg.noise_seed = config.noise_seeds(si);
        % find_thresholds 会固定 rng(0) 初始化；noise_seed 仅控制局部噪声流。
        % 网络显式注入；所有配置和 r 复用同一个已读取的标准斑图初值。
        try
            thresholds = find_thresholds('ER', config.er_p, cfg);
            if any(~isfinite(thresholds)), error('阈值结果非有限值。'); end
            results(j,:) = thresholds;
            errors{j} = '';
            fprintf('[SCAN DONE] %s r=%.2f seed=%d -> sf=%.4f, sb=%.4f\n', ...
                cases(c).tag, config.ratios(ri), config.noise_seeds(si), ...
                thresholds(1), thresholds(2));
        catch err
            errors{j} = err.message;
        end
    end
    failures = {};
    for j = 1:numel(ids)
        if isempty(errors{j})
            state.sf_repeats(ids(j)) = results(j,1);
            state.sb_repeats(ids(j)) = results(j,2);
        else
            [c,si,ri] = ind2sub([3 num_seeds num_r], ids(j));
            failures{end+1} = sprintf('%s, r=%.3g, seed=%d：%s', ...
                cases(c).label, config.ratios(ri), config.noise_seeds(si), errors{j}); %#ok<SAGROW>
        end
    end
    finished = isfinite(state.sf_repeats) & isfinite(state.sb_repeats);
    state.elapsed_seconds = previous_elapsed+toc(run_clock);
    checkpoint(result_file, state);
    fprintf('[SAVED] 阈值已完成 %d/%d，累计 %.1f s。\n', ...
        nnz(finished), numel(finished), state.elapsed_seconds);
    stopOnFailures(failures, result_file);
end

%% 5. 仅沿噪声重复维计算统计；谱分析没有噪声重复误差
state.sf_mean = reshape(mean(state.sf_repeats,2), 3, num_r);
state.sb_mean = reshape(mean(state.sb_repeats,2), 3, num_r);
state.sf_var = reshape(var(state.sf_repeats,0,2), 3, num_r); % 样本方差，分母 n-1
state.sb_var = reshape(var(state.sb_repeats,0,2), 3, num_r);
state.sf_std = sqrt(state.sf_var);
state.sb_std = sqrt(state.sb_var);
near_boundary = state.sf_repeats <= config.sigma_min+0.05 | ...
    state.sf_repeats >= config.sigma_max-0.05 | ...
    state.sb_repeats <= config.sigma_min+0.05 | ...
    state.sb_repeats >= config.sigma_max-0.05;
if any(near_boundary(:))
    warning('%d 个任务阈值接近搜索边界，不能认定区间内存在分叉；请检查并扩大 sigma 区间。', nnz(near_boundary));
end
if any(state.sf_mean(:) < state.sb_mean(:))
    warning('部分正向均值低于反向均值；未交换实际数据，仅对相区填色边界排序。');
end
state.statistics_note = 'Independent noise realizations on fixed topology; shared precomputed ER seed for all cases and ratios; sample variance n-1; bands mean +/- 1 SD, not CI.';
if RUN_SIMULATION, checkpoint(result_file, state); end
output_tag = [topology_key '_sharedseed_' experiment_tag];
image_files = fig2c_topology_plot(state, cases, fullfile(script_dir,'fig'), output_tag);
fprintf('\n=== 图 2c 固定拓扑噪声重复分析完成 ===\n结果：%s\n', result_file);
for i = 1:numel(image_files), fprintf('图 %d：%s\n', i, image_files{i}); end

%% 局部工具：不修改项目现有求解器和阈值函数
function item = computeSpectrum(net, ratio, base_cfg)
item = struct('spectrum',NaN,'lambda_star',NaN,'error','');
try
    [item.spectrum,item.lambda_star] = fig2c_topology_spectrum(net.L, ...
        base_cfg.L_inter, base_cfg.alpha, ratio);
    fprintf('[SPECTRUM DONE] %s r=%.2f spectrum=%.4f\n',net.tag,ratio,item.spectrum);
catch err
    item.error = err.message;
end
end

function workers = setupPool(requested)
workers = 0;
if requested <= 0
    fprintf('[POOL] REQUESTED_WORKERS=0，串行补算。\n');
    return;
end
if ~license('test','Distrib_Computing_Toolbox')
    warning('无并行工具箱，将串行补算。');
    return;
end
try
    pool = gcp('nocreate');
    if isempty(pool)
        cluster = parcluster('local');
        pool = parpool(cluster, min(requested,cluster.NumWorkers));
    end
    workers = pool.NumWorkers;
    fprintf('[POOL] 实际并行 worker 数：%d。\n', workers);
catch err
    warning('并行池启动失败，将串行补算：%s', err.message);
end
end

function checkpoint(path, state)
temporary_file = [path '.tmp.mat'];
save(temporary_file, '-struct', 'state', '-v7');
[ok,msg] = movefile(temporary_file,path,'f');
if ~ok, error('断点文件替换失败：%s',msg); end
end

function stopOnFailures(failures, result_file)
if ~isempty(failures)
    error('本批有 %d 个失败任务，成功任务已经保存，可重新运行补算。\n%s\n断点：%s', ...
        numel(failures), strjoin(failures,newline), result_file);
end
end

function validateNetworks(cases, L_inter, config)
assert(numel(cases)==3 && isequal(size(L_inter),[config.K config.K]) && ...
    norm(L_inter-(config.K*eye(config.K)-ones(config.K)),'fro')<1e-10, '网络缓存层间结构错误。');
for c = 1:3
    assert(numel(cases(c).L)==config.K, '网络缓存层数错误。');
    for k = 1:config.K
        L = cases(c).L{k};
        assert(isequal(size(L),[config.N config.N]) && ...
            norm(L-L','fro')<1e-10 && norm(full(sum(L,2)))<1e-8, '网络缓存拉普拉斯矩阵错误。');
    end
end
assert(all(ismember({'WS','BA','ER'},cases(3).layer_types)), '混合配置必须包含 WS、BA、ER。');
end
