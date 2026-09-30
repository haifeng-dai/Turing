%% scan_ba_lambda_star.m：BA 网络主导谱点参数扫描
% 层间为度 2 的最近邻环；alpha=0.05，beta/alpha=0.1。
% 扫描 600 组 (N,m,K)，仅读取项目已有 BA 网络，不生成或修改拓扑。
% lambda_star 是使谱阈值最小的正非均匀谱点，不是最大特征值。
% MAT 中各结果变量为等长列向量，每行对应一个 (N,m,K) 组合；支持断点续算。
% 热图保存到 fig。只读绘图时将 RUN_COMPUTATION=false；修改扫描范围请编辑下方配置。
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
RUN_COMPUTATION = true;
MAKE_FIGURE = true;
REQUESTED_WORKERS = 64;
BATCH_SIZE = 64;

config.N_values = 100:100:1000;
config.m_values = [2:10 20];
config.K_values = [3 4 5 6 8 10];
config.alpha = 0.05;
config.ratio = 0.1;
config.inter_topology = 'unweighted_degree2_ring';
config.schema_version = 1;
config.spectrum_options = struct('dense_limit', 600, ...
    'dense_verify_limit', 1200, 'eigs_counts', [16 32 64 128 256], ...
    'eigs_tol', 1e-10, 'eigs_maxit', 2000, ...
    'residual_tol', 1e-7, 'comparison_tol', 1e-8);

network_dir = fullfile(script_dir, 'results', 'topology', 'BA');
results_dir = fullfile(script_dir, 'results');
figure_dir = fullfile(script_dir, 'fig');
validateSettings(config, REQUESTED_WORKERS, BATCH_SIZE);
source_signature = fileread(fullfile(script_dir, 'ba_lambda_star_spectrum.m'));
file_tag = sprintf('N%s_m%s_K%s_a%.3f_r%.3f', ...
    vectorTag(config.N_values), vectorTag(config.m_values), ...
    vectorTag(config.K_values), config.alpha, config.ratio);
file_tag = strrep(file_tag, '.', 'p');
file_stem = ['scan_ba_lambda_star_' file_tag];
legacy_file_stem = ['scan_lambda_star_ba_ring_' file_tag];
result_file = fullfile(results_dir, [file_stem '.mat']);
legacy_result_file = fullfile(results_dir, [legacy_file_stem '.mat']);
figure_file = fullfile(figure_dir, ['ba_lambda_star_' file_tag '.png']);

%% 1. 严格读取：在启动并行池前检查全部网络和文件指纹
[pair_N, pair_m] = ndgrid(config.N_values, config.m_values);
manifest = struct('name', {}, 'sha256', {}, 'layers', {});
network_files = cell(numel(pair_N), 1);
for group_idx = 1:numel(pair_N)
    name = sprintf('N%d_m%d.mat', pair_N(group_idx), pair_m(group_idx));
    path = fullfile(network_dir, name);
    if ~isfile(path)
        error('scan_ba_lambda_star:MissingNetwork', ...
            '缺少已有 BA 网络：%s。不会自动生成，请先准备对应网络。', path);
    end
    variables = whos('-file', path);
    idx = find(strcmp({variables.name}, 'nets'), 1);
    assert(~isempty(idx) && strcmp(variables(idx).class, 'cell'), ...
        '%s 缺少 cell 类型 nets。', path);
    count = prod(variables(idx).size);
    assert(count >= max(config.K_values), '%s 只有 %d 层，需要至少 %d 层。', ...
        path, count, max(config.K_values));
    manifest(group_idx) = struct('name', name, 'sha256', sha256File(path), 'layers', count);
    network_files{group_idx} = path;
end
fprintf('[网络] 已检查 %d 个 BA 网络文件。\n', numel(manifest));

%% 2. 读取或初始化 MAT 断点缓存
if isfile(result_file)
    state = load(result_file);
    assert(isfield(state, 'config') && isequaln(state.config, config) && ...
        isfield(state, 'manifest') && isequaln(state.manifest, manifest) && ...
        isfield(state, 'source_signature') && strcmp(state.source_signature, source_signature), ...
        ['缓存与当前参数、算法或网络文件不同，拒绝混用。' ...
        '请调整文件名或保留对应实验配置：%s'], result_file);
    expected_rows = numel(config.N_values)*numel(config.m_values)*numel(config.K_values);
    assert(numel(state.N) == expected_rows && numel(state.m) == expected_rows && ...
        numel(state.K) == expected_rows && numel(state.lambda_star) == expected_rows && ...
        numel(state.done) == expected_rows, ...
        '缓存结果维度错误。');
    fprintf('[缓存] 已读取 %d/%d 个完成点。\n', nnz(state.done), numel(state.done));
    migrated_legacy = false;
elseif isfile(legacy_result_file)
    legacy_state = load(legacy_result_file);
    state = migrateLegacyState(legacy_state, config, manifest, source_signature);
    migrated_legacy = true;
    fprintf('[缓存] 已读取旧版 MAT 断点：%s\n', legacy_result_file);
    fprintf('[缓存] 后续断点将保存为项目命名的新 MAT 文件，不修改旧文件。\n');
else
    state = initialState(config, manifest, source_signature);
    migrated_legacy = false;
end
reported_result_file = result_file;
if migrated_legacy && ~RUN_COMPUTATION
    reported_result_file = legacy_result_file;
end
fprintf('[设置] 共 %d 组；alpha=%.3g，beta=%.3g，r=%.3g；层间为度 2 环。\n', ...
    numel(state.done), config.alpha, config.alpha*config.ratio, config.ratio);
if ~RUN_COMPUTATION && ~all(state.done)
    error('scan_ba_lambda_star:IncompleteCache', ...
        '只读模式：尚有 %d 个未完成点，不启动计算。', nnz(~state.done));
end
if migrated_legacy && RUN_COMPUTATION
    if ~isfolder(results_dir), mkdir(results_dir); end
    checkpoint(result_file, state);
end

%% 3. (N,m) 分组并行：每个任务只读取一次网络，计算该组尚缺的 K
jobs = struct('N', {}, 'm', {}, 'file', {}, 'rows', {}, 'K_values', {});
for group_idx = 1:numel(pair_N)
    rows = find(~state.done & state.N == pair_N(group_idx) & ...
        state.m == pair_m(group_idx));
    if ~isempty(rows)
        jobs(end+1) = struct('N', pair_N(group_idx), 'm', pair_m(group_idx), ...
            'file', network_files{group_idx}, 'rows', rows, ...
            'K_values', state.K(rows)); %#ok<SAGROW>
    end
end
run_clock = tic;
previous_elapsed = state.elapsed_seconds;
workers = 0;
if ~isempty(jobs)
    if ~isfolder(results_dir), mkdir(results_dir); end
    checkpoint(result_file, state);
    workers = setupPool(REQUESTED_WORKERS, min(numel(jobs), BATCH_SIZE));
    fprintf('[并行] 实际进程数 %d；待计算分组 %d；每批最多 %d 组。\n', ...
        workers, numel(jobs), BATCH_SIZE);
else
    fprintf('[缓存] 无需补算，不启动并行池。\n');
end
for first = 1:BATCH_SIZE:numel(jobs)
    batch_jobs = jobs(first:min(first+BATCH_SIZE-1, numel(jobs)));
    fprintf('[并行] 本批 %d 组；已完成 %d/%d 个谱点。\n', ...
        numel(batch_jobs), nnz(state.done), numel(state.done));
    batch_results = cell(numel(batch_jobs), 1);
    parfor (j = 1:numel(batch_jobs), workers)
        batch_results{j} = computeGroup(batch_jobs(j), config);
    end
    failures = {};
    for j = 1:numel(batch_jobs)
        group = batch_results{j};
        for ki = 1:numel(group)
            row = batch_jobs(j).rows(ki);
            if isempty(group(ki).error)
                item = group(ki).value;
                state.lambda_star(row) = item.lambda_star;
                state.sigma_c_spec(row) = item.sigma_c_spec;
                state.lambda_opt(row) = item.lambda_opt;
                state.lambda_offset(row) = item.lambda_offset;
                state.mean_degree(row) = item.mean_degree;
                state.relative_residual(row) = item.relative_residual;
                state.neighborhood_bound(row) = item.neighborhood_bound;
                state.computed_modes(row) = item.computed_modes;
                state.seconds(row) = item.seconds;
                state.dense_verified(row) = item.dense_verified;
                state.method{row} = item.method;
                state.status{row} = item.status;
                state.done(row) = true;
            else
                failures{end+1} = sprintf('N=%d m=%d K=%d: %s', ...
                    batch_jobs(j).N, batch_jobs(j).m, batch_jobs(j).K_values(ki), ...
                    group(ki).error); %#ok<SAGROW>
                state.status{row} = 'error';
            end
        end
    end
    state.elapsed_seconds = previous_elapsed+toc(run_clock);
    checkpoint(result_file, state);
    fprintf('[保存] 已完成 %d/%d 点；累计 %.1f 秒；%s\n', ...
        nnz(state.done), numel(state.done), state.elapsed_seconds, result_file);
    if ~isempty(failures)
        error('scan_ba_lambda_star:TaskFailed', ...
            '已保存成功结果；修正后重新运行会补算失败项。\n%s', strjoin(failures, '\n'));
    end
end

%% 4. 保存 MAT 结果并按需绘制统一色标热图
if MAKE_FIGURE && ~isfolder(figure_dir), mkdir(figure_dir); end
if MAKE_FIGURE, plotResults(state, config, figure_file); end
finite_points = isfinite(state.lambda_star);
valid_lambda_opt = state.lambda_opt(isfinite(state.lambda_opt));
fprintf('\n=== BA 主导谱点扫描完成 ===\n');
fprintf('完成 %d/%d 组；有效主导点 %d。\n', ...
    nnz(state.done), numel(state.done), nnz(finite_points));
fprintf('失败点 %d；无可接受非均匀模态点 %d。\n', ...
    nnz(strcmp(state.status, 'error')), nnz(strcmp(state.status, 'no_admissible_mode')));
if ~isempty(valid_lambda_opt)
    fprintf('连续阈值函数最优谱位置 lambda_opt=%.8f。\n', valid_lambda_opt(1));
end
if any(finite_points)
    fprintf('lambda_star 范围：[%.8f, %.8f]\n', ...
        min(state.lambda_star(finite_points)), max(state.lambda_star(finite_points)));
end
fprintf('MAT 结果文件：%s\n', reported_result_file);
if MAKE_FIGURE, fprintf('热图文件：%s\n', figure_file); end

function validateSettings(config, requested_workers, batch_size)
grids = {config.N_values, config.m_values, config.K_values};
for i = 1:numel(grids)
    v = grids{i};
    assert(isnumeric(v) && isvector(v) && ~isempty(v) && isreal(v) && ...
        all(isfinite(v)) && all(v == floor(v)) && all(diff(v) > 0), ...
        'N,m,K 网格须是严格递增的有限整数向量。');
end
assert(min(config.N_values) >= 2 && min(config.m_values) >= 1 && ...
    max(config.m_values) < min(config.N_values), 'BA 参数要求 1<=m<N。');
assert(min(config.K_values) >= 3, '简单无向度 2 环要求 K>=3，不允许 K=1,2。');
validateattributes(config.alpha, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(config.ratio, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(requested_workers, {'numeric'}, {'scalar','integer','nonnegative'});
validateattributes(batch_size, {'numeric'}, {'scalar','integer','positive'});
options = config.spectrum_options;
assert(options.dense_limit >= 0 && options.dense_verify_limit >= options.dense_limit, ...
    '完整谱验证上限必须不小于直接完整谱计算上限。');
assert(isvector(options.eigs_counts) && numel(options.eigs_counts) >= 2 && ...
    all(options.eigs_counts >= 2) && all(diff(options.eigs_counts) > 0) && ...
    all(options.eigs_counts == floor(options.eigs_counts)), ...
    'eigs_counts 至少含两个递增整数。');
validateattributes(options.eigs_tol, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(options.eigs_maxit, {'numeric'}, {'scalar','integer','positive'});
validateattributes(options.residual_tol, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(options.comparison_tol, {'numeric'}, {'scalar','real','finite','positive'});
end

function tag = vectorTag(values)
parts = arrayfun(@(value) sprintf('%g', value), values, 'UniformOutput', false);
tag = strjoin(parts, '-');
end

function state = initialState(config, manifest, source_signature)
[N, m, K] = ndgrid(config.N_values, config.m_values, config.K_values);
N = N(:); m = m(:); K = K(:); count = numel(N);
state.config = config;
state.manifest = manifest;
state.source_signature = source_signature;
state.N = N;
state.m = m;
state.K = K;
state.alpha = repmat(config.alpha, count, 1);
state.r = repmat(config.ratio, count, 1);
state.beta = repmat(config.alpha*config.ratio, count, 1);
state.inter_degree = repmat(2, count, 1);
state.lambda_star = nan(count, 1);
state.sigma_c_spec = nan(count, 1);
state.lambda_opt = nan(count, 1);
state.lambda_offset = nan(count, 1);
state.mean_degree = nan(count, 1);
state.relative_residual = nan(count, 1);
state.neighborhood_bound = nan(count, 1);
state.computed_modes = nan(count, 1);
state.seconds = nan(count, 1);
state.dense_verified = false(count, 1);
state.method = repmat({''}, count, 1);
state.status = repmat({'pending'}, count, 1);
state.network_file = arrayfun(@(n, a) sprintf('N%d_m%d.mat', n, a), ...
    N, m, 'UniformOutput', false);
state.done = false(count, 1);
state.elapsed_seconds = 0;
end

function state = migrateLegacyState(legacy_state, config, manifest, source_signature)
required_state = {'config','manifest','source_signature','results','done','elapsed_seconds'};
for i = 1:numel(required_state)
    assert(isfield(legacy_state, required_state{i}), ...
        '旧版 MAT 缓存缺少字段 %s。', required_state{i});
end
assert(isequaln(legacy_state.config, config) && ...
    isequaln(legacy_state.manifest, manifest) && ...
    strcmp(legacy_state.source_signature, source_signature), ...
    '旧版 MAT 缓存与当前参数、算法或 BA 网络文件不匹配，拒绝续算。');
assert(istable(legacy_state.results), '旧版缓存的 results 不是项目旧格式 table。');

state = initialState(config, manifest, source_signature);
expected_N = state.N;
expected_m = state.m;
expected_K = state.K;
expected_rows = numel(expected_N);
assert(height(legacy_state.results) == expected_rows && ...
    numel(legacy_state.done) == expected_rows, '旧版 MAT 缓存结果维度错误。');
fields = {'N','m','K','alpha','r','beta','inter_degree','lambda_star', ...
    'sigma_c_spec','lambda_opt','lambda_offset','mean_degree', ...
    'relative_residual','neighborhood_bound','computed_modes','seconds', ...
    'dense_verified','method','status','network_file'};
for i = 1:numel(fields)
    field = fields{i};
    assert(ismember(field, legacy_state.results.Properties.VariableNames), ...
        '旧版 MAT 缓存结果表缺少列 %s。', field);
    state.(field) = legacy_state.results.(field);
end
assert(isequaln(state.N, expected_N) && isequaln(state.m, expected_m) && ...
    isequaln(state.K, expected_K), ...
    '旧版 MAT 缓存中的 N、m、K 排列与当前扫描顺序不同。');
state.done = logical(legacy_state.done(:));
state.elapsed_seconds = legacy_state.elapsed_seconds;
end

function group = computeGroup(job, config)
group = repmat(struct('value', struct, 'error', ''), numel(job.K_values), 1);
try
    data = load(job.file, 'nets');
    layers = data.nets(1:max(job.K_values));
    degrees = zeros(numel(layers), 1);
    expected_degree_sum = 2*job.m*job.N-job.m*(job.m+1);
    for k = 1:numel(layers)
        L = sparse(layers{k});
        assert(isequal(size(L), [job.N job.N]) && isreal(L) && ...
            all(isfinite(nonzeros(L))), '第 %d 层矩阵尺寸或数值无效。', k);
        assert(norm(L-L', 1) < 1e-8 && max(abs(sum(L, 2))) < 1e-8, ...
            '第 %d 层不是对称拉普拉斯。', k);
        d = full(diag(L));
        offdiag = nonzeros(L-spdiags(d, 0, job.N, job.N));
        assert(all(d >= 0) && all(abs(offdiag+1) < 1e-8) && ...
            abs(sum(d)-expected_degree_sum) < 1e-8, ...
            '第 %d 层不是符合 N=%d,m=%d 的项目无权 BA 网络。', k, job.N, job.m);
        layers{k} = L;
        degrees(k) = mean(d);
    end
catch err
    for i = 1:numel(group), group(i).error = err.message; end
    return;
end
for i = 1:numel(job.K_values)
    clock = tic;
    try
        K = job.K_values(i);
        item = ba_lambda_star_spectrum(layers(1:K), config.alpha, ...
            config.ratio, config.spectrum_options);
        item.lambda_offset = item.lambda_star-item.lambda_opt;
        item.mean_degree = mean(degrees(1:K));
        item.seconds = toc(clock);
        group(i).value = item;
    catch err
        group(i).error = err.message;
    end
end
end

function count = setupPool(requested, task_count)
count = 0;
if requested == 0, return; end
assert(license('test', 'Distrib_Computing_Toolbox'), ...
    '并行计算需要 Parallel Computing Toolbox；串行可设 REQUESTED_WORKERS=0。');
pool = gcp('nocreate');
if isempty(pool)
    cluster = parcluster('local');
    target = min([requested, cluster.NumWorkers, task_count]);
    pool = parpool(cluster, target);
end
count = min([requested, pool.NumWorkers, task_count]);
end

function digest = sha256File(path)
fid = fopen(path, 'rb');
assert(fid >= 0, '无法读取网络文件：%s', path);
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
hasher = java.security.MessageDigest.getInstance('SHA-256');
hasher.update(bytes);
digest = lower(reshape(dec2hex(typecast(hasher.digest(), 'uint8'), 2)', 1, []));
end

function checkpoint(path, state)
parent = fileparts(path);
temporary = [tempname(parent) '.mat'];
cleanup = onCleanup(@() removeTemporary(temporary));
save(temporary, '-struct', 'state', '-v7');
[ok, message] = movefile(temporary, path, 'f');
assert(ok, '保存缓存失败：%s', message);
end

function removeTemporary(path)
if isfile(path), delete(path); end
end

function plotResults(state, config, path)
values = reshape(state.lambda_star, ...
    numel(config.N_values), numel(config.m_values), numel(config.K_values));
finite_values = values(isfinite(values));
assert(~isempty(finite_values), '没有有效主导点，不能绘制热图；请查看 MAT 文件中的 status 字段。');
limits = [min(finite_values), max(finite_values)];
if limits(1) == limits(2), limits = limits+[-1 1]*max(1, abs(limits(1)))*1e-6; end
panels = numel(config.K_values);
cols = min(3, panels);
h = figure('Color', 'w', 'Position', [100 100 1450 780]);
layout = tiledlayout(h, ceil(panels/cols), cols, 'TileSpacing', 'compact', 'Padding', 'compact');
for ki = 1:panels
    ax = nexttile(layout);
    matrix = values(:,:,ki)';
    pic = imagesc(ax, 1:numel(config.N_values), 1:numel(config.m_values), matrix);
    pic.AlphaData = isfinite(matrix);
    ax.Color = [0.85 0.85 0.85];
    set(ax, 'YDir', 'normal', 'FontName', 'Arial', ...
        'XTick', 1:numel(config.N_values), 'XTickLabel', config.N_values, ...
        'YTick', 1:numel(config.m_values), 'YTickLabel', config.m_values);
    xlabel(ax, 'N'); ylabel(ax, 'm'); title(ax, sprintf('K = %d', config.K_values(ki)));
    caxis(ax, limits);
end
colormap(h, parula(256));
bar = colorbar(ax);
bar.Layout.Tile = 'east';
bar.Label.String = '\lambda_*';
title(layout, sprintf('BA layers | degree-2 interlayer ring | alpha=%.3g, beta/alpha=%.3g', ...
    config.alpha, config.ratio), 'Interpreter', 'none');
exportgraphics(h, path, 'Resolution', 300);
close(h);
end
