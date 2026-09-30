%% scan_lambda_star_fig3_5.m：图 3–5 拓扑配置下的主导谱点变化
% 仅读取 results/topology 中已有的层内拉普拉斯；缺失文件或层数不足时报错。
% 图 3、4 使用完全图层间耦合；图 5 使用度为 2 的最近邻环。
% K=2 时以权重 2 的双边实现度为 2；不生成拓扑、不扫描随机种子、不保存 CSV。
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

RUN_COMPUTATION = true;
RUN_PARALLEL = true;
REQUESTED_WORKERS = 8;
BATCH_SIZE = 8;
MAKE_FIGURE = true;

config.alpha = 0.05;
config.ratio = 0.1;  % beta/alpha；对应 beta=0.005
config.spectrum_options = struct('dense_limit', 600, ...
    'dense_verify_limit', 1200, 'eigs_counts', [16 32 64 128 256], ...
    'eigs_tol', 1e-10, 'eigs_maxit', 2000, ...
    'residual_tol', 1e-7, 'comparison_tol', 1e-8);
config.interlayer_fig3_fig4 = 'complete_graph';
config.interlayer_fig5 = 'unweighted_degree2_ring';
config.ring_K2_policy = 'double_edge_weight_2';

results_dir = fullfile(script_dir, 'results');
topology_dir = fullfile(results_dir, 'topology');
figure_dir = fullfile(script_dir, 'fig');
result_file = fullfile(results_dir, ...
    'scan_lambda_star_fig3_5_a0p050_r0p100.mat');

definitions = buildCases(topology_dir);
manifest = inspectTopologyFiles(definitions);
config.spectrum_function_sha256 = sha256File( ...
    fullfile(script_dir, 'multiplex_lambda_star_spectrum.m'));
fprintf('[网络] 已验证 %d 个不同拓扑文件；只读取现有 nets。\n', numel(manifest));

if isfile(result_file)
    saved = load(result_file, 'state');
    assert(isfield(saved, 'state') && isfield(saved.state, 'schema_version') && ...
        saved.state.schema_version == 1 && ...
        isequaln(saved.state.config, config) && ...
        isequaln(saved.state.definitions, definitions) && ...
        isequaln(saved.state.manifest, manifest), ...
        ['已有 MAT 缓存与当前参数、谱算法或拓扑文件不同，拒绝混用：%s。' ...
        '保留旧文件并为新配置指定独立文件名。'], result_file);
    state = saved.state;
    assert(numel(state.done) == numel(definitions) && ...
        numel(state.values) == numel(definitions), '缓存结果数量不匹配。');
    fprintf('[缓存] 已完成 %d/%d 个谱点。\n', nnz(state.done), numel(state.done));
else
    state = initialState(config, definitions, manifest);
end

missing = find(~state.done);
if ~RUN_COMPUTATION && ~isempty(missing)
    error('scan_lambda_star_fig3_5:IncompleteCache', ...
        '只读计算模式下仍缺少 %d 个谱点；未启动计算。', numel(missing));
end

if RUN_COMPUTATION && ~isempty(missing)
    if ~isfolder(results_dir), mkdir(results_dir); end
    saveCheckpoint(result_file, state);
    workers = setupPool(RUN_PARALLEL, REQUESTED_WORKERS, numel(missing));
    fprintf('[并行] 实际 worker=%d，待计算谱点=%d，每批最多=%d。\n', ...
        workers, numel(missing), BATCH_SIZE);
    clock_all = tic;
    for first = 1:BATCH_SIZE:numel(missing)
        batch = missing(first:min(first+BATCH_SIZE-1, numel(missing)));
        output = cell(numel(batch), 1);
        if workers > 0
            parfor (j = 1:numel(batch), workers)
                output{j} = computeCase(definitions(batch(j)), config);
            end
        else
            for j = 1:numel(batch)
                output{j} = computeCase(definitions(batch(j)), config);
            end
        end

        failures = {};
        for j = 1:numel(batch)
            idx = batch(j);
            if isempty(output{j}.error)
                state.values(idx) = output{j}.value;
                state.done(idx) = true;
                state.errors{idx} = '';
            else
                state.errors{idx} = output{j}.error;
                failures{end+1} = sprintf('%s / %s / N=%d K=%d: %s', ...
                    definitions(idx).study, definitions(idx).topology, ...
                    definitions(idx).N, definitions(idx).K, output{j}.error); %#ok<SAGROW>
            end
        end
        state.elapsed_seconds = state.elapsed_seconds + toc(clock_all);
        clock_all = tic;
        saveCheckpoint(result_file, state);
        fprintf('[保存] 完成 %d/%d 点；累计 %.1f 秒：%s\n', ...
            nnz(state.done), numel(state.done), state.elapsed_seconds, result_file);
        if ~isempty(failures)
            error('scan_lambda_star_fig3_5:TaskFailed', ...
                '成功结果已保存；修复错误后重新运行会补算失败点。\n%s', ...
                strjoin(failures, '\n'));
        end
    end
elseif RUN_COMPUTATION
    fprintf('[缓存] 所有谱点均已完成，不启动并行池。\n');
end

assert(all(state.done), '仍有谱点未完成，不能生成最终图。');
summary_table = makeSummaryTable(definitions, state.values);
state.summary_table = summary_table;
if ~isfolder(results_dir), mkdir(results_dir); end
saveCheckpoint(result_file, state);

if MAKE_FIGURE
    if ~isfolder(figure_dir), mkdir(figure_dir); end
    plotStudy(definitions, state.values, 'fig3_connectivity', ...
        'Fig. 3: matched intralayer mean degree', ...
        'Mean intralayer degree', 'Complete interlayer graph', ...
        fullfile(figure_dir, 'lambda_star_fig3_connectivity_a0p050_r0p100.png'));
    plotStudy(definitions, state.values, 'fig4_network_size', ...
        'Fig. 4: network size', 'N', 'Complete interlayer graph', ...
        fullfile(figure_dir, 'lambda_star_fig4_network_size_a0p050_r0p100.png'));
    plotStudy(definitions, state.values, 'fig5_layer_count', ...
        'Fig. 5: layer count', 'K', ...
        'Degree-2 interlayer ring (K=2: double edge, weight 2)', ...
        fullfile(figure_dir, 'lambda_star_fig5_layer_count_a0p050_r0p100.png'));
end

fprintf('\n=== 图 3–5 主导谱点计算完成 ===\n');
fprintf('alpha=%.4g, beta/alpha=%.4g (beta=%.4g)；完成 %d/%d 点。\n', ...
    config.alpha, config.ratio, config.alpha*config.ratio, ...
    nnz(state.done), numel(state.done));
fprintf('结果 MAT：%s\n', result_file);
if MAKE_FIGURE
    fprintf('图片目录：%s\n', figure_dir);
end
disp(summary_table);

function definitions = buildCases(topology_dir)
definitions = struct('study', {}, 'topology', {}, 'N', {}, 'K', {}, ...
    'parameter_name', {}, 'parameter_value', {}, 'parameter_label', {}, ...
    'target_mean_degree', {}, 'interlayer_type', {}, 'topology_file', {});

% 图 3：平均度相近的 ER、WS、BA 拓扑配置；固定 N=200、K=5。
degree_targets = [6 8 20];
er_p = [0.03 0.04 0.10];
ws_k = [6 8 20];
ba_m = [3 4 10];
for i = 1:numel(degree_targets)
    definitions(end+1, 1) = makeDefinition('fig3_connectivity', 'ER', ...
        200, 5, 'p', er_p(i), sprintf('p=%.3f', er_p(i)), ...
        degree_targets(i), 'complete_graph', ...
        fullfile(topology_dir, 'ER', sprintf('N200_p%.3f.mat', er_p(i)))); %#ok<AGROW>
    definitions(end+1, 1) = makeDefinition('fig3_connectivity', 'WS', ...
        200, 5, 'k', ws_k(i), sprintf('k=%d, p_{rewire}=0.01', ws_k(i)), ...
        degree_targets(i), 'complete_graph', ...
        fullfile(topology_dir, 'WS', sprintf('N200_K%d_pr0.01.mat', ws_k(i)))); %#ok<AGROW>
    definitions(end+1, 1) = makeDefinition('fig3_connectivity', 'BA', ...
        200, 5, 'm', ba_m(i), sprintf('m=%d', ba_m(i)), ...
        degree_targets(i), 'complete_graph', ...
        fullfile(topology_dir, 'BA', sprintf('N200_m%d.mat', ba_m(i)))); %#ok<AGROW>
end

% 图 4：固定拓扑参数与 K=5，扫描网络规模 N。
n_values = [100 500 1000];
for N = n_values
    definitions(end+1, 1) = makeDefinition('fig4_network_size', 'ER', ...
        N, 5, 'p', 0.03, 'p=0.030', NaN, 'complete_graph', ...
        fullfile(topology_dir, 'ER', sprintf('N%d_p0.030.mat', N))); %#ok<AGROW>
    definitions(end+1, 1) = makeDefinition('fig4_network_size', 'WS', ...
        N, 5, 'k', 6, 'k=6, p_{rewire}=0.10', NaN, 'complete_graph', ...
        fullfile(topology_dir, 'WS', sprintf('N%d_K6_pr0.10.mat', N))); %#ok<AGROW>
    definitions(end+1, 1) = makeDefinition('fig4_network_size', 'BA', ...
        N, 5, 'm', 3, 'm=3', NaN, 'complete_graph', ...
        fullfile(topology_dir, 'BA', sprintf('N%d_m3.mat', N))); %#ok<AGROW>
end

% 图 5：固定 N=200 与各拓扑的基准参数，扫描层数 K；层间改为度 2 环。
k_values = [2 6 10];
for K = k_values
    definitions(end+1, 1) = makeDefinition('fig5_layer_count', 'ER', ...
        200, K, 'p', 0.03, 'p=0.030', NaN, 'degree2_ring', ...
        fullfile(topology_dir, 'ER', 'N200_p0.030.mat')); %#ok<AGROW>
    definitions(end+1, 1) = makeDefinition('fig5_layer_count', 'WS', ...
        200, K, 'k', 6, 'k=6, p_{rewire}=0.10', NaN, 'degree2_ring', ...
        fullfile(topology_dir, 'WS', 'N200_K6_pr0.10.mat')); %#ok<AGROW>
    definitions(end+1, 1) = makeDefinition('fig5_layer_count', 'BA', ...
        200, K, 'm', 3, 'm=3', NaN, 'degree2_ring', ...
        fullfile(topology_dir, 'BA', 'N200_m3.mat')); %#ok<AGROW>
end
end

function item = makeDefinition(study, topology, N, K, parameter_name, ...
    parameter_value, parameter_label, target_mean_degree, interlayer_type, topology_file)
item = struct('study', study, 'topology', topology, 'N', N, 'K', K, ...
    'parameter_name', parameter_name, 'parameter_value', parameter_value, ...
    'parameter_label', parameter_label, 'target_mean_degree', target_mean_degree, ...
    'interlayer_type', interlayer_type, 'topology_file', topology_file);
end

function manifest = inspectTopologyFiles(definitions)
files = unique({definitions.topology_file});
manifest = repmat(struct('file', '', 'sha256', '', ...
    'bytes', 0, 'available_layers', 0), numel(files), 1);
for i = 1:numel(files)
    path = files{i};
    if ~isfile(path)
        error('scan_lambda_star_fig3_5:MissingTopology', ...
            '缺少图 3–5 使用的已有拓扑，不会自动生成：%s', path);
    end
    variables = whos('-file', path);
    idx = find(strcmp({variables.name}, 'nets'), 1);
    assert(~isempty(idx) && strcmp(variables(idx).class, 'cell'), ...
        '拓扑文件缺少 cell 类型 nets：%s', path);
    available_layers = prod(variables(idx).size);
    required = max([definitions(strcmp({definitions.topology_file}, path)).K]);
    assert(available_layers >= required, ...
        '拓扑文件 %s 只有 %d 层，图 3–5 需要至少 %d 层。', ...
        path, available_layers, required);
    info = dir(path);
    manifest(i) = struct('file', path, 'sha256', sha256File(path), ...
        'bytes', info.bytes, 'available_layers', available_layers);
end
end

function state = initialState(config, definitions, manifest)
value_template = struct('x_value', NaN, 'mean_degree', NaN, ...
    'lambda_star', NaN, 'sigma_c_spec', NaN, 'lambda_opt', NaN, ...
    'relative_residual', NaN, 'neighborhood_bound', NaN, ...
    'computed_modes', 0, 'dense_verified', false, ...
    'method', 'pending', 'status', 'pending', 'seconds', NaN);
state = struct;
state.schema_version = 1;
state.config = config;
state.definitions = definitions;
state.manifest = manifest;
state.values = repmat(value_template, numel(definitions), 1);
state.done = false(numel(definitions), 1);
state.errors = repmat({''}, numel(definitions), 1);
state.elapsed_seconds = 0;
end

function count = setupPool(run_parallel, requested, task_count)
count = 0;
if ~run_parallel || requested == 0, return; end
assert(license('test', 'Distrib_Computing_Toolbox'), ...
    '并行计算需要 Parallel Computing Toolbox；可将 RUN_PARALLEL 设为 false。');
pool = gcp('nocreate');
if isempty(pool)
    cluster = parcluster('local');
    target = min([requested, cluster.NumWorkers, task_count]);
    assert(target >= 1, '本地并行集群没有可用 worker。');
    pool = parpool(cluster, target);
end
count = min([requested, pool.NumWorkers, task_count]);
end

function output = computeCase(definition, config)
value = struct('x_value', NaN, 'mean_degree', NaN, ...
    'lambda_star', NaN, 'sigma_c_spec', NaN, 'lambda_opt', NaN, ...
    'relative_residual', NaN, 'neighborhood_bound', NaN, ...
    'computed_modes', 0, 'dense_verified', false, ...
    'method', '', 'status', 'error', 'seconds', NaN);
output = struct('value', value, 'error', '');
clock = tic;
try
    data = load(definition.topology_file, 'nets');
    if ~isfield(data, 'nets') || ~iscell(data.nets) || ...
            numel(data.nets) < definition.K
        error('拓扑文件缺少所需层数：%s', definition.topology_file);
    end
    layers = data.nets(1:definition.K);
    degrees = zeros(definition.K, 1);
    for k = 1:definition.K
        layers{k} = sparse(layers{k});
        if ~isequal(size(layers{k}), [definition.N definition.N])
            error('第 %d 层尺寸与 N=%d 不符：%s', ...
                k, definition.N, definition.topology_file);
        end
        degrees(k) = mean(diag(layers{k}));
    end
    L_inter = makeInterlayerLaplacian(definition.K, definition.interlayer_type);
    spectrum = multiplex_lambda_star_spectrum(layers, L_inter, ...
        config.alpha, config.ratio, config.spectrum_options);
    value.lambda_star = spectrum.lambda_star;
    value.sigma_c_spec = spectrum.sigma_c_spec;
    value.lambda_opt = spectrum.lambda_opt;
    value.relative_residual = spectrum.relative_residual;
    value.neighborhood_bound = spectrum.neighborhood_bound;
    value.computed_modes = spectrum.computed_modes;
    value.dense_verified = spectrum.dense_verified;
    value.method = spectrum.method;
    value.status = spectrum.status;
    value.mean_degree = mean(degrees);
    switch definition.study
        case 'fig3_connectivity'
            value.x_value = value.mean_degree;
        case 'fig4_network_size'
            value.x_value = definition.N;
        case 'fig5_layer_count'
            value.x_value = definition.K;
        otherwise
            error('未知扫描类型：%s', definition.study);
    end
    value.seconds = toc(clock);
    output.value = value;
catch err
    output.error = err.message;
end
end

function L_inter = makeInterlayerLaplacian(K, topology)
switch topology
    case 'complete_graph'
        L_inter = K*speye(K)-sparse(ones(K));
    case 'degree2_ring'
        assert(K >= 2, '度 2 环至少需要两层。');
        adjacency = sparse(1:K, [2:K 1], 1, K, K);
        adjacency = adjacency+adjacency';
        L_inter = 2*speye(K)-adjacency;
    otherwise
        error('不支持的层间拓扑：%s', topology);
end
end

function saveCheckpoint(path, state)
parent = fileparts(path);
if ~isfolder(parent), mkdir(parent); end
temporary = [tempname(parent) '.mat'];
cleanup = onCleanup(@() removeTemporary(temporary));
save(temporary, 'state', '-v7');
[ok, message] = movefile(temporary, path, 'f');
assert(ok, '保存 MAT 检查点失败：%s', message);
end

function removeTemporary(path)
if isfile(path), delete(path); end
end

function digest = sha256File(path)
fid = fopen(path, 'rb');
assert(fid >= 0, '无法读取文件：%s', path);
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
hasher = java.security.MessageDigest.getInstance('SHA-256');
hasher.update(bytes);
digest = lower(reshape(dec2hex(typecast(hasher.digest(), 'uint8'), 2)', 1, []));
end

function summary = makeSummaryTable(definitions, values)
study = {definitions.study}';
topology = {definitions.topology}';
parameter_label = {definitions.parameter_label}';
interlayer_type = {definitions.interlayer_type}';
summary = table(study, topology, [definitions.N]', [definitions.K]', ...
    {definitions.parameter_name}', [definitions.parameter_value]', ...
    parameter_label, [values.x_value]', [values.mean_degree]', ...
    [values.lambda_star]', [values.sigma_c_spec]', ...
    interlayer_type, {values.method}', {values.status}', [values.seconds]', ...
    'VariableNames', {'study', 'topology', 'N', 'K', 'parameter_name', ...
    'parameter_value', 'parameter_label', 'x_value', 'mean_degree', ...
    'lambda_star', 'sigma_c_spec', 'interlayer_type', 'method', 'status', 'seconds'});
end

function plotStudy(definitions, values, study, title_text, x_label, ...
    subtitle_text, output_file)
topologies = {'ER', 'WS', 'BA'};
markers = {'o', 's', '^'};
colors = lines(numel(topologies));
figure_handle = figure('Color', 'w', 'Position', [100 100 850 620]);
axes_handle = axes(figure_handle);
hold(axes_handle, 'on');
for family_idx = 1:numel(topologies)
    selected = find(strcmp({definitions.study}, study) & ...
        strcmp({definitions.topology}, topologies{family_idx}));
    x = [values(selected).x_value]';
    y = [values(selected).lambda_star]';
    [x, order] = sort(x);
    y = y(order);
    plot(axes_handle, x, y, ['-' markers{family_idx}], ...
        'Color', colors(family_idx, :), 'LineWidth', 2, ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(family_idx, :), ...
        'DisplayName', topologies{family_idx});
end
xlabel(axes_handle, x_label, 'FontSize', 16);
ylabel(axes_handle, '\lambda_*', 'FontSize', 16, 'Interpreter', 'tex');
title(axes_handle, {title_text, subtitle_text}, ...
    'FontSize', 15, 'FontWeight', 'normal', 'Interpreter', 'tex');
legend(axes_handle, 'Location', 'best', 'FontSize', 13);
grid(axes_handle, 'on');
box(axes_handle, 'on');
set(axes_handle, 'FontSize', 13, 'LineWidth', 1.0);
exportgraphics(figure_handle, output_file, 'Resolution', 300);
close(figure_handle);
fprintf('[图] 已保存：%s\n', output_file);
end
