%% fig2c_noise_species_ablation.m：原图 2(c) 的 u/v 噪声消融
% 读取已存在的 ER 双组分噪声基准，只扫描 u-only 与 v-only 两种新条件。
% 固定单个噪声种子，不做随机种子重复、均值或方差估计。
% 原求解器、原图 2(c) 数据和标准斑图初态均保持不变。
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(script_dir,'networks'), fullfile(script_dir,'simulations'));

RUN_SIMULATION = true;
REQUESTED_WORKERS = 64;
BATCH_SIZE = 64;

config.N = 200;
config.K = 5;
config.topology = 'ER';
config.p = 0.030;
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
% 原 solve_multiplex 在未指定局部种子时以 rng(1024) 产生噪声。
% 用同一固定种子作单次消融对比；这不是多种子扫描。
config.noise_seed = 1024;
config.schema_version = 3;

conditions = {'u_only','v_only'};
noise_masks = [1 0; 0 1];
num_conditions = numel(conditions);
num_ratios = numel(config.ratios);

%% 1. 读取原图 2(c) 已有的双组分 ER 基准曲线
baseline_file = fullfile(script_dir,'results', ...
    'scan_beta_alpha_er_N200_K5_p0.030.mat');
if ~isfile(baseline_file)
    error('缺少原图 2(c) 的双组分基准文件：%s',baseline_file);
end
baseline = load(baseline_file,'SF_Matrix','SB_Matrix','ALPHA_LIST', ...
    'RATIO_LIST','P_VAL_FIXED');
required = {'SF_Matrix','SB_Matrix','ALPHA_LIST','RATIO_LIST','P_VAL_FIXED'};
for i = 1:numel(required)
    if ~isfield(baseline,required{i})
        error('原图 2(c) 基准文件缺少字段 %s：%s',required{i},baseline_file);
    end
end
if baseline.P_VAL_FIXED ~= config.p
    error('基准网络 p=%.6g 与当前配置 p=%.6g 不匹配。', ...
        baseline.P_VAL_FIXED,config.p);
end
alpha_index = find(abs(baseline.ALPHA_LIST-config.alpha) <= 1e-12);
if numel(alpha_index) ~= 1
    error('基准数据中未唯一找到 alpha=%.3f。',config.alpha);
end
if ~isequal(size(baseline.SF_Matrix),[numel(baseline.ALPHA_LIST),numel(baseline.RATIO_LIST)]) || ...
        ~isequal(size(baseline.SB_Matrix),size(baseline.SF_Matrix))
    error('原图 2(c) 基准阈值矩阵与 alpha/ratio 坐标尺寸不匹配。');
end
if numel(baseline.RATIO_LIST) ~= num_ratios || ...
        any(abs(baseline.RATIO_LIST(:)'-config.ratios) > 1e-12)
    error('原图 2(c) 基准的 beta/alpha 横轴与当前扫描范围不一致。');
end
baseline_sf = reshape(baseline.SF_Matrix(alpha_index,:),1,[]);
baseline_sb = reshape(baseline.SB_Matrix(alpha_index,:),1,[]);
if any(~isfinite(baseline_sf)) || any(~isfinite(baseline_sb))
    error('原图 2(c) 双组分基准曲线包含非有限阈值，拒绝绘图。');
end

%% 2. 读取与原图 2(c) 一致的固定 ER 网络和预生成斑图初态
network_file = fullfile(script_dir,'results','topology','ER', ...
    sprintf('N%d_p%.3f.mat',config.N,config.p));
if ~isfile(network_file)
    error('缺少原图 2(c) 的 ER 网络文件：%s',network_file);
end
network = load(network_file,'nets');
if ~isfield(network,'nets') || ~iscell(network.nets) || numel(network.nets) < config.K
    error('ER 网络文件没有至少 %d 层可用网络：%s',config.K,network_file);
end
for k = 1:config.K
    if ~isequal(size(network.nets{k}),[config.N config.N])
        error('ER 网络第 %d 层的矩阵尺寸不是 %d×%d。',k,config.N,config.N);
    end
end
[shared_seed,seed_file] = load_standard_pattern_seed(config,false);
L_inter = sparse(config.K*eye(config.K)-ones(config.K));

%% 3. 消融专用断点；与之前三拓扑/五种子缓存隔离
experiment_tag = sprintf('er_N%d_K%d_p%.3f_a%.3f_eta%.3f_seed%d', ...
    config.N,config.K,config.p,config.alpha,config.noise,config.noise_seed);
result_file = fullfile(script_dir,'results', ...
    ['scan_noise_species_ablation_' experiment_tag '.mat']);
image_file = fullfile(script_dir,'fig', ...
    ['fig2c_noise_species_ablation_' experiment_tag '.png']);

source_signature.thresholds = fileread(which('find_thresholds_species_noise'));
source_signature.solver = fileread(which('solve_multiplex_species_noise'));
source_signature.shared_seed = shared_seed;
source_signature.protocol_version = 3;

if isfile(result_file)
    state = load(result_file);
    if ~isfield(state,'config') || ~isequaln(state.config,config) || ...
            ~isfield(state,'conditions') || ~isequal(state.conditions,conditions) || ...
            ~isfield(state,'noise_masks') || ~isequal(state.noise_masks,noise_masks) || ...
            ~isfield(state,'source_signature') || ...
            ~isequaln(state.source_signature,source_signature)
        error('消融缓存与当前参数或求解代码不匹配，拒绝混用：%s',result_file);
    end
    assert(isequal(size(state.sf),[num_conditions num_ratios]) && ...
        isequal(size(state.sb),[num_conditions num_ratios]) && ...
        isequal(size(state.done),[num_conditions num_ratios]), ...
        '消融缓存数组尺寸错误。');
    fprintf('[缓存] 读取消融断点：%s\n',result_file);
else
    state.config = config;
    state.conditions = conditions;
    state.noise_masks = noise_masks;
    state.source_signature = source_signature;
    state.baseline_file = baseline_file;
    state.network_file = network_file;
    state.seed_file = seed_file;
    state.shared_seed = shared_seed;
    state.sf = nan(num_conditions,num_ratios);
    state.sb = nan(num_conditions,num_ratios);
    state.done = false(num_conditions,num_ratios);
    state.elapsed_seconds = 0;
end

pending = find(~state.done);
fprintf('[设置] 1 种 ER 拓扑；双组分基准直接读取，不重算。\n');
fprintf('[设置] 新扫描仅含 u-only、v-only；固定 noise_seed=%d；无种子重复。\n', ...
    config.noise_seed);
fprintf('[设置] 待算 %d/%d 个阈值任务（每个任务含正、反向阈值）。\n', ...
    numel(pending),numel(state.done));
if ~RUN_SIMULATION && ~isempty(pending)
    error('只读模式：还有 %d 个消融阈值任务未完成。',numel(pending));
end
if ~isfolder(fullfile(script_dir,'results')), mkdir(fullfile(script_dir,'results')); end
if ~isfile(result_file) && RUN_SIMULATION, checkpoint(result_file,state); end

%% 4. 并行补算两种单组分噪声条件
workers = 0;
if RUN_SIMULATION && ~isempty(pending)
    workers = setupPool(REQUESTED_WORKERS,min(BATCH_SIZE,numel(pending)));
else
    fprintf('[缓存] 无待算任务，不启动并行池。\n');
end
base_cfg = struct('N',config.N,'K',config.K,'alpha',config.alpha, ...
    'noise',config.noise,'T_END',config.T_END,'dt',config.dt, ...
    'steps',config.steps,'init_perturb',config.init_perturb, ...
    'detect_convergence',config.detect_convergence, ...
    'sigma_min',config.sigma_min,'sigma_max',config.sigma_max, ...
    'L_inter',L_inter,'L_intra',{network.nets(1:config.K)}, ...
    'y_seed',shared_seed,'noise_seed',config.noise_seed);
timer = tic;
previous_elapsed = state.elapsed_seconds;
for first = 1:BATCH_SIZE:numel(pending)
    ids = pending(first:min(first+BATCH_SIZE-1,numel(pending)));
    values = nan(numel(ids),2);
    messages = cell(numel(ids),1);
    fprintf('[并行] 本批 %d 个任务；已完成 %d/%d。\n', ...
        numel(ids),nnz(state.done),numel(state.done));
    parfor (j = 1:numel(ids),workers)
        [mi,ri] = ind2sub([num_conditions num_ratios],ids(j));
        cfg = base_cfg;
        cfg.beta = config.alpha*config.ratios(ri);
        cfg.noise_species_mask = noise_masks(mi,:);
        try
            result = find_thresholds_species_noise('ER',config.p,cfg);
            if any(~isfinite(result)), error('阈值结果非有限值。'); end
            values(j,:) = result;
            messages{j} = '';
        catch err
            messages{j} = err.message;
        end
    end
    failures = {};
    for j = 1:numel(ids)
        if isempty(messages{j})
            state.sf(ids(j)) = values(j,1);
            state.sb(ids(j)) = values(j,2);
            state.done(ids(j)) = true;
        else
            [mi,ri] = ind2sub([num_conditions num_ratios],ids(j));
            failures{end+1} = sprintf('%s，固定种子=%d，beta/alpha=%.3g：%s', ...
                conditions{mi},config.noise_seed,config.ratios(ri), ...
                messages{j}); %#ok<SAGROW>
        end
    end
    state.elapsed_seconds = previous_elapsed+toc(timer);
    checkpoint(result_file,state);
    fprintf('[保存] 已完成 %d/%d；累计 %.1f 秒。\n', ...
        nnz(state.done),numel(state.done),state.elapsed_seconds);
    if ~isempty(failures)
        error('fig2c_noise_species_ablation:TaskFailed', ...
            '已保存成功结果；重新运行可补算失败项。\n%s',strjoin(failures,newline));
    end
end

if any(~state.done(:))
    error('消融数据不完整，不生成对比图。');
end
state.baseline_sf = baseline_sf;
state.baseline_sb = baseline_sb;
state.baseline_file = baseline_file;
state.statistics_note = ['Single fixed noise seed; no seed averaging. ' ...
    'The both-species curve is read from the existing Fig. 2(c) ER scan.'];
if RUN_SIMULATION, checkpoint(result_file,state); end
if ~isfolder(fullfile(script_dir,'fig')), mkdir(fullfile(script_dir,'fig')); end
fig2c_noise_species_plot(state,image_file);
fprintf('\n=== 图 2(c) u/v 噪声消融完成 ===\n');
fprintf('双组分基准（复用）：%s\n',baseline_file);
fprintf('两种单组分新结果：%s\n',result_file);
fprintf('图：%s\n',image_file);

function workers = setupPool(requested,task_count)
workers = 0;
if requested <= 0, return; end
if ~license('test','Distrib_Computing_Toolbox')
    warning('未检测到并行工具箱，改为串行计算。');
    return;
end
try
    pool = gcp('nocreate');
    if isempty(pool)
        cluster = parcluster('local');
        pool = parpool(cluster,min([requested,cluster.NumWorkers,task_count]));
    end
    workers = min([requested,pool.NumWorkers,task_count]);
catch err
    warning('并行池启动失败，改为串行计算：%s',err.message);
end
end

function checkpoint(path,state)
temporary = [tempname(fileparts(path)) '.mat'];
cleanup = onCleanup(@() removeTemporary(temporary));
save(temporary,'-struct','state','-v7');
[ok,message] = movefile(temporary,path,'f');
assert(ok,'断点保存失败：%s',message);
end

function removeTemporary(path)
if isfile(path), delete(path); end
end
