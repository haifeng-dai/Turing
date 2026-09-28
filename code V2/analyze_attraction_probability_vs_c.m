%% analyze_attraction_probability_vs_c.m：沿均匀态到斑图态方向估计经验吸引概率
% 复用 Figure 2 的网络、参数、随机动力学积分器和序参量定义。
% 不生成网络；不改写模型或噪声构造；不保存完整随机仿真时间序列。
clear; clc; close all;

%% 1. 定位 Figure 2 参数、阈值数据和已有网络
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(script_dir, 'simulations'), ...
    fullfile(script_dir, 'networks'));
fig2_script = fullfile(script_dir, 'fig2abc.m');
threshold_script = fullfile(script_dir, 'simulations', 'find_thresholds.m');
solver_script = fullfile(script_dir, 'simulations', 'solve_multiplex.m');
if ~exist(fig2_script, 'file') || ~exist(threshold_script, 'file') || ...
        ~exist(solver_script, 'file')
    error('找不到 Figure 2、阈值搜索或模型积分器脚本。');
end

fig2_source = fileread(fig2_script);
threshold_source = fileread(threshold_script);
topology_type = readTextAssignment(fig2_source, 'TOPO_TYPE');
N = readNumericAssignment(fig2_source, 'N');
K = readNumericAssignment(fig2_source, 'K');
topology_parameter = readNumericAssignment(fig2_source, 'P_VAL_FIXED');
noise_intensity = readNumericAssignment(fig2_source, 'DYNA.noise');
T_END = readNumericAssignment(fig2_source, 'DYNA.T_END');
figure2_steps = readNumericAssignment(fig2_source, 'DYNA.steps');
init_perturb = readNumericAssignment(fig2_source, 'DYNA.init_perturb');
dt = readNumericAssignment(threshold_source, 'cfg.dt');
classification_threshold = readNumericAssignment(threshold_source, 'thresh_A');

results_dir = fullfile(script_dir, 'results');
fig_dir = fullfile(script_dir, 'fig');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
threshold_file = fullfile(results_dir, sprintf( ...
    'scan_beta_alpha_%s_N%d_K%d_p%.3f.mat', lower(topology_type), ...
    N, K, topology_parameter));
if ~exist(threshold_file, 'file')
    error('找不到 Figure 2 阈值数据：%s', threshold_file);
end

threshold_data = load(threshold_file, 'SF_Matrix', 'SB_Matrix', ...
    'ALPHA_LIST', 'RATIO_LIST', 'BETA_LIST', 'P_VAL_FIXED');
if ~isfield(threshold_data, 'SF_Matrix') || ...
        ~isfield(threshold_data, 'SB_Matrix') || ...
        ~isfield(threshold_data, 'ALPHA_LIST')
    error('阈值 MAT 文件缺少 SF_Matrix、SB_Matrix 或 ALPHA_LIST。');
end
if isfield(threshold_data, 'RATIO_LIST')
    ratio_axis = threshold_data.RATIO_LIST(:)';
elseif isfield(threshold_data, 'BETA_LIST')
    % Figure 2 当前文件应保存 beta/alpha；此分支兼容旧版文件。
    ratio_axis = threshold_data.BETA_LIST(:)' ./ 0.1;
else
    error('阈值 MAT 文件中找不到 RATIO_LIST 或 BETA_LIST。');
end

alpha_requested = 0.10;
[alpha_error, alpha_idx] = min(abs(threshold_data.ALPHA_LIST(:)' - alpha_requested));
if alpha_error > 1e-8
    error('Figure 2 阈值文件中找不到 alpha=%.3f。', alpha_requested);
end
alpha = threshold_data.ALPHA_LIST(alpha_idx);
if size(threshold_data.SF_Matrix, 1) < alpha_idx || ...
        size(threshold_data.SB_Matrix, 1) < alpha_idx || ...
        size(threshold_data.SF_Matrix, 2) ~= numel(ratio_axis) || ...
        size(threshold_data.SB_Matrix, 2) ~= numel(ratio_axis)
    error('SF_Matrix/SB_Matrix 尺寸与 ALPHA_LIST、ratio 轴不一致。');
end

requested_R_LIST = [0.1 4.5 8.0];
R_LIST = zeros(size(requested_R_LIST));
sigma_f = zeros(size(requested_R_LIST));
sigma_b = zeros(size(requested_R_LIST));
for r_idx = 1:numel(requested_R_LIST)
    [~, ratio_idx] = min(abs(ratio_axis - requested_R_LIST(r_idx)));
    R_LIST(r_idx) = ratio_axis(ratio_idx);
    if abs(R_LIST(r_idx)-requested_R_LIST(r_idx)) > 1e-10
        fprintf('请求 r=%.6g，阈值数据中采用最近点 r=%.6g。\n', ...
            requested_R_LIST(r_idx), R_LIST(r_idx));
    end
    sigma_f(r_idx) = threshold_data.SF_Matrix(alpha_idx, ratio_idx);
    sigma_b(r_idx) = threshold_data.SB_Matrix(alpha_idx, ratio_idx);
end
if any(~isfinite(sigma_f) | ~isfinite(sigma_b))
    error('所选 beta/alpha 点的 forward/backward threshold 含非有限值。');
end
if any(sigma_b >= sigma_f)
    error(['所选参数中至少一个不满足 sigma_b < sigma_f，无法将其作为滞回区间分析：' ...
        '\nsigma_b=%s\nsigma_f=%s'], mat2str(sigma_b, 6), mat2str(sigma_f, 6));
end

sigma_common_low = max(sigma_b);
sigma_common_high = min(sigma_f);
if sigma_common_low < sigma_common_high
    sigma_selection_mode = 'common hysteresis intersection';
    sigma_test = repmat((sigma_common_low + sigma_common_high) / 2, ...
        size(R_LIST));
    matched_midpoint_note = '';
else
    sigma_selection_mode = 'individual hysteresis midpoints';
    sigma_test = (sigma_b + sigma_f) / 2;
    matched_midpoint_note = ...
        'Each coupling ratio is evaluated at the midpoint of its own hysteresis interval.';
end

topology_file = fullfile(results_dir, 'topology', upper(topology_type), ...
    sprintf('N%d_p%.3f.mat', N, topology_parameter));
if ~exist(topology_file, 'file')
    error('找不到 Figure 2 使用的已有网络文件：%s', topology_file);
end
network_data = load(topology_file, 'nets');
if ~isfield(network_data, 'nets') || numel(network_data.nets) < K
    error('网络文件的层数不足：%s', topology_file);
end
if ~iscell(network_data.nets)
    error('网络文件中的 nets 不是 cell 数组，和现有积分器输入格式不符。');
end
L_intra = network_data.nets(1:K);
for layer_idx = 1:K
    if ~isequal(size(L_intra{layer_idx}), [N N])
        error('网络第 %d 层尺寸不是 N×N。', layer_idx);
    end
    L_intra{layer_idx} = sparse(L_intra{layer_idx});
end
adj_inter = ones(K) - eye(K);
L_inter = diag(sum(adj_inter, 2)) - adj_inter;

%% 2. 仿真设置与初值方向
% 以下参数与 Figure 2(c) 及其阈值搜索入口一致。
R_LIST = R_LIST(:)';
C_LIST = 0:0.05:1;
N_REALIZATIONS = 100;
SEED_BASE = 20260927;
N_WORKERS_REQUESTED = 64;
INITIAL_PERTURBATION_RELATIVE_RMS = 1e-3;
CONTINUATION_MAX_STEP = 5.0;
CONTINUATION_MIN_STEP = 0.05;
Z_WILSON = 1.96;

% Figure 2 的积分器采用 Euler-Maruyama；steps=2 仅保存起点和终点，
% A_final 仍由积分器对末尾 20% 时间步求平均。
base_cfg = struct();
base_cfg.N = N;
base_cfg.K = K;
base_cfg.alpha = alpha;
base_cfg.noise = noise_intensity;
base_cfg.dt = dt;
base_cfg.T_END = T_END;
base_cfg.steps = 2;
base_cfg.init_perturb = init_perturb;
base_cfg.L_intra = L_intra;
base_cfg.L_inter = L_inter;
base_cfg.detect_convergence = false;

NK = N * K;
x_H = [5 * ones(NK, 1); 10 * ones(NK, 1)];
direction_stream = RandStream('mt19937ar', 'Seed', SEED_BASE);
zeta_u = randn(direction_stream, NK, 1);
zeta_v = randn(direction_stream, NK, 1);
zeta_u = zeta_u - mean(zeta_u);
zeta_v = zeta_v - mean(zeta_v);
zeta = [zeta_u; zeta_v];
zeta = zeta / norm(zeta);
epsilon_initial = INITIAL_PERTURBATION_RELATIVE_RMS * norm(x_H);

% 优先选同网络、同 alpha/beta/noise 且 sigma 不低于目标值的已存 backward 状态；
% 若不存在，则使用同网络的最高 sigma 已存斑图作为延拓起点。
pattern_seed_state = cell(size(R_LIST));
pattern_seed_sigma = NaN(size(R_LIST));
pattern_seed_source = cell(size(R_LIST));
for r_idx = 1:numel(R_LIST)
    beta = alpha * R_LIST(r_idx);
    [pattern_seed_state{r_idx}, pattern_seed_sigma(r_idx), ...
        pattern_seed_source{r_idx}] = findPatternSeed(results_dir, ...
        topology_type, N, K, topology_parameter, alpha, beta, ...
        noise_intensity, sigma_test(r_idx), L_intra, classification_threshold);
end

%% 3. 生成并确认各个 coupling ratio 的斑图参考状态
timer_total = tic;
pool = startLocalPool(N_WORKERS_REQUESTED);
actual_workers = pool.NumWorkers;
fprintf('Requested workers: %d\nActual workers: %d\n', ...
    N_WORKERS_REQUESTED, actual_workers);

n_ratio = numel(R_LIST);
pattern_reference_states = cell(1, n_ratio);
reference_A_final = NaN(1, n_ratio);
reference_continuation_steps = zeros(1, n_ratio);
reference_source = pattern_seed_source;
parfor r_idx = 1:n_ratio
    cfg_r = base_cfg;
    cfg_r.beta = alpha * R_LIST(r_idx);
    [x_pattern_local, A_reference_local, continuation_count_local] = buildPatternReference( ...
        cfg_r, pattern_seed_state{r_idx}, pattern_seed_sigma(r_idx), ...
        sigma_test(r_idx), classification_threshold, SEED_BASE, r_idx, ...
        CONTINUATION_MAX_STEP, CONTINUATION_MIN_STEP);
    pattern_reference_states{r_idx} = x_pattern_local;
    reference_A_final(r_idx) = A_reference_local;
    reference_continuation_steps(r_idx) = continuation_count_local;
end
if any(reference_A_final <= classification_threshold)
    error('至少一个斑图参考状态未通过现有 A_cut=%.6g 分类阈值。', ...
        classification_threshold);
end
for r_idx = 1:n_ratio
    if reference_A_final(r_idx) <= 5 * classification_threshold
        warning('r=%.3f 的参考斑图振幅仅为 A_cut 的 %.2f 倍，请检查状态分离。', ...
            R_LIST(r_idx), reference_A_final(r_idx) / classification_threshold);
    end
end

%% 4. 扁平任务表与独立随机 realization
n_c = numel(C_LIST);
n_jobs = n_ratio * n_c * N_REALIZATIONS;
[job_ratio_grid, job_c_grid, job_realization_grid] = ndgrid( ...
    1:n_ratio, 1:n_c, 1:N_REALIZATIONS);
job_ratio = job_ratio_grid(:);
job_c = job_c_grid(:);
job_realization = job_realization_grid(:);
job_seed = SEED_BASE + job_realization;
A_final_flat = NaN(n_jobs, 1);

parfor job_id = 1:n_jobs
    r_idx = job_ratio(job_id);
    c_idx = job_c(job_id);
    c_value = C_LIST(c_idx);
    x_P = pattern_reference_states{r_idx};
    y0 = x_H + c_value * (x_P - x_H) + epsilon_initial * zeta;

    cfg = base_cfg;
    cfg.beta = alpha * R_LIST(r_idx);
    cfg.sigma = sigma_test(r_idx);
    cfg.y0 = y0;
    cfg.noise_seed = job_seed(job_id);
    cfg.detect_convergence = false;
    [~, ~, cfg_out] = solve_multiplex(cfg);
    A_final_flat(job_id) = cfg_out.A_final;
end
A_final_all = reshape(A_final_flat, [n_ratio, n_c, N_REALIZATIONS]);
total_wall_time = toc(timer_total);
n_continuation_integrations = sum(reference_continuation_steps);
n_total_simulations = n_jobs + n_continuation_integrations;

%% 5. 分类概率、Wilson 区间与 c50
n_pattern = sum(A_final_all > classification_threshold, 3);
P_pattern = n_pattern / N_REALIZATIONS;
CI_low = zeros(size(P_pattern));
CI_high = zeros(size(P_pattern));
for r_idx = 1:n_ratio
    for c_idx = 1:n_c
        [CI_low(r_idx, c_idx), CI_high(r_idx, c_idx)] = wilsonInterval( ...
            n_pattern(r_idx, c_idx), N_REALIZATIONS, Z_WILSON);
    end
end

c50 = NaN(1, n_ratio);
for r_idx = 1:n_ratio
    c50(r_idx) = interpolateHalfProbability(C_LIST, P_pattern(r_idx, :));
    if isnan(c50(r_idx))
        fprintf('r=%.3f 的概率曲线未穿过 0.5，c50 未定义。\n', R_LIST(r_idx));
    end
end

%% 6. 绘制经验吸引概率曲线
curve_colors = lines(n_ratio);
main_figure = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [100 100 900 620]);
main_axes = axes(main_figure);
hold(main_axes, 'on');
for r_idx = 1:n_ratio
    lower_error = P_pattern(r_idx, :) - CI_low(r_idx, :);
    upper_error = CI_high(r_idx, :) - P_pattern(r_idx, :);
    errorbar(main_axes, C_LIST, P_pattern(r_idx, :), lower_error, upper_error, ...
        '-o', 'Color', curve_colors(r_idx, :), ...
        'MarkerFaceColor', curve_colors(r_idx, :), 'MarkerSize', 5, ...
        'LineWidth', 1.5, 'CapSize', 5, ...
        'DisplayName', sprintf('r = %.3g', R_LIST(r_idx)));
end
yline(main_axes, 0.5, '--', 'Color', [0.55 0.55 0.55], ...
    'HandleVisibility', 'off');
xlim(main_axes, [0 1]);
ylim(main_axes, [0 1]);
xlabel(main_axes, 'Initial pattern fraction c');
ylabel(main_axes, 'Probability of patterned final state');
title(main_axes, 'Empirical attraction probability along a representative direction');
if ~isempty(matched_midpoint_note)
    subtitle(main_axes, matched_midpoint_note, 'Interpreter', 'none');
end
legend(main_axes, 'Location', 'best');
grid(main_axes, 'on');
box(main_axes, 'on');
main_image = fullfile(fig_dir, 'attraction_probability_vs_c.png');
exportgraphics(main_figure, main_image, 'Resolution', 300);
close(main_figure);

c50_image = '';
if all(isfinite(c50))
    c50_figure = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [120 120 720 520]);
    plot(R_LIST, c50, '-o', 'LineWidth', 1.7, 'MarkerSize', 7, ...
        'MarkerFaceColor', [0.20 0.48 0.78]);
    xlabel('Inter-layer coupling ratio r = beta / alpha');
    ylabel('c_{50}');
    title('Empirical half-probability point');
    xlim([min(R_LIST) max(R_LIST)]);
    ylim([0 1]);
    grid on; box on;
    c50_image = fullfile(fig_dir, 'c50_vs_coupling.png');
    exportgraphics(c50_figure, c50_image, 'Resolution', 300);
    close(c50_figure);
end

%% 7. 保存 MAT、CSV 和运行摘要
row_count = n_ratio * n_c;
[csv_ratio_grid, csv_c_grid] = ndgrid(1:n_ratio, 1:n_c);
csv_ratio_idx = csv_ratio_grid(:);
csv_c_idx = csv_c_grid(:);
csv_r = reshape(R_LIST(csv_ratio_idx), [], 1);
csv_c = reshape(C_LIST(csv_c_idx), [], 1);
csv_sigma = reshape(sigma_test(csv_ratio_idx), [], 1);
csv_table = table(csv_r, csv_c, csv_sigma, n_pattern(:), ...
    repmat(N_REALIZATIONS, row_count, 1), P_pattern(:), CI_low(:), CI_high(:), ...
    'VariableNames', {'r', 'c', 'sigma', 'n_pattern', 'n_total', ...
    'P_pattern', 'CI_low', 'CI_high'});
csv_file = fullfile(results_dir, 'attraction_probability_vs_c.csv');
writetable(csv_table, csv_file);

simulation_parameters = struct( ...
    'topology_type', topology_type, 'topology_parameter', topology_parameter, ...
    'network_file', topology_file, 'alpha', alpha, 'noise_intensity', noise_intensity, ...
    'N', N, 'K', K, 'dt', dt, 'T_END', T_END, 'steps_saved', base_cfg.steps, ...
    'figure2_steps_setting', figure2_steps, 'init_perturb', init_perturb, ...
    'order_parameter_definition', ...
    'sqrt(sum((u-5).^2+(v-10).^2)/(N*K)), final 20 percent time-step mean', ...
    'noise_definition', ...
    'solve_multiplex.m edge-wise independent increments; same edge increment shared by u/v');
mat_file = fullfile(results_dir, 'attraction_probability_vs_c.mat');
save(mat_file, 'requested_R_LIST', 'R_LIST', 'C_LIST', 'N_REALIZATIONS', ...
    'sigma_b', 'sigma_f', 'sigma_test', 'sigma_selection_mode', ...
    'sigma_common_low', 'sigma_common_high', 'P_pattern', 'CI_low', 'CI_high', ...
    'n_pattern', 'c50', 'A_final_all', 'classification_threshold', ...
    'epsilon_initial', 'zeta', 'x_H', 'pattern_reference_states', ...
    'reference_A_final', 'reference_continuation_steps', 'reference_source', ...
    'pattern_seed_sigma', 'pattern_seed_source', 'SEED_BASE', 'job_seed', ...
    'simulation_parameters', 'threshold_file', 'topology_file', ...
    'N_WORKERS_REQUESTED', 'actual_workers', 'n_jobs', ...
    'n_continuation_integrations', 'n_total_simulations', 'total_wall_time', ...
    'matched_midpoint_note', 'main_image', 'csv_file', 'c50_image', '-v7');

fprintf('\n=== Attraction probability analysis complete ===\n');
fprintf('Figure 2 parameters: alpha=%.6g, eta=%.6g, N=%d, K=%d\n', ...
    alpha, noise_intensity, N, K);
fprintf('Network file: %s\n', topology_file);
fprintf('Threshold data: %s\n', threshold_file);
fprintf('Actual r values: %s\n', mat2str(R_LIST, 6));
for r_idx = 1:n_ratio
    fprintf('r=%.6g: sigma_b=%.8g, sigma_f=%.8g, sigma_test=%.8g, c50=', ...
        R_LIST(r_idx), sigma_b(r_idx), sigma_f(r_idx), sigma_test(r_idx));
    if isfinite(c50(r_idx))
        fprintf('%.6g\n', c50(r_idx));
    else
        fprintf('未定义（曲线未穿过 0.5）\n');
    end
    fprintf('  斑图参考态：A_final=%.6g，延拓积分次数=%d\n  种子文件：%s\n', ...
        reference_A_final(r_idx), reference_continuation_steps(r_idx), ...
        reference_source{r_idx});
end
if strcmp(sigma_selection_mode, 'common hysteresis intersection')
    fprintf('Sigma selection: common hysteresis value %.8g\n', sigma_test(1));
else
    fprintf(['Each coupling ratio is evaluated at the midpoint of its own ' ...
        'hysteresis interval.\n']);
end
fprintf('Requested workers: %d; actual workers: %d\n', ...
    N_WORKERS_REQUESTED, actual_workers);
fprintf(['Total simulations: %d (Monte Carlo: %d; continuation integrations: %d); ' ...
    'elapsed: %.2f s\n'], n_total_simulations, n_jobs, ...
    n_continuation_integrations, total_wall_time);
fprintf('MAT: %s\nCSV: %s\nFigure: %s\n', mat_file, csv_file, main_image);
if ~isempty(c50_image)
    fprintf('c50 figure: %s\n', c50_image);
end

%% 本地辅助函数
function value = readNumericAssignment(source_text, variable_name)
pattern = ['(?m)^\s*' regexptranslate('escape', variable_name) ...
    '\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*;'];
token = regexp(source_text, pattern, 'tokens', 'once');
if isempty(token)
    error('无法从当前 MATLAB 源文件中读取参数赋值：%s', variable_name);
end
value = str2double(token{1});
end

function value = readTextAssignment(source_text, variable_name)
pattern = ['(?m)^\s*' regexptranslate('escape', variable_name) ...
    '\s*=\s*''([^'']+)''\s*;'];
token = regexp(source_text, pattern, 'tokens', 'once');
if isempty(token)
    error('无法从当前 MATLAB 源文件中读取文本参数：%s', variable_name);
end
value = token{1};
end

function pool = startLocalPool(requested_workers)
pool = gcp('nocreate');
if ~isempty(pool)
    fprintf('检测到已有并行池，将继续使用，不中断现有 worker。\n');
    return;
end
cluster = parcluster('local');
profile_limit = cluster.NumWorkers;
first_attempt = min(requested_workers, profile_limit);
last_error = '';
for worker_count = first_attempt:-1:1
    try
        pool = parpool(cluster, worker_count);
        if worker_count < requested_workers
            fprintf('64 个 worker 未能启动，已降级到当前可启动的最大数量。\n');
        end
        return;
    catch pool_error
        last_error = pool_error.message;
        existing_pool = gcp('nocreate');
        if ~isempty(existing_pool)
            delete(existing_pool);
        end
    end
end
error('无法启动本地并行池。最后一次错误：%s', last_error);
end

function [seed_state, seed_sigma, seed_path] = findPatternSeed(results_dir, ...
        topology_type, N, K, topology_parameter, alpha, beta, noise_intensity, ...
        sigma_target, L_intra, A_cut)
pattern = ['^evolution_' lower(topology_type) '_N(\d+)_K(\d+)_p([0-9.]+)' ...
    '_a([0-9.]+)_b([0-9.]+)_s([0-9.]+)_n([0-9.]+)_(fwd|bwd)_results\.mat$'];
files = dir(fullfile(results_dir, sprintf('evolution_%s_N%d_K%d_p*_results.mat', ...
    lower(topology_type), N, K)));
candidates = struct('path', {}, 'sigma', {}, 'exact', {}, 'state', {});
for file_idx = 1:numel(files)
    token = regexp(files(file_idx).name, pattern, 'tokens', 'once');
    if isempty(token)
        continue;
    end
    file_N = str2double(token{1});
    file_K = str2double(token{2});
    file_param = str2double(token{3});
    file_alpha = str2double(token{4});
    file_beta = str2double(token{5});
    file_sigma = str2double(token{6});
    file_noise = str2double(token{7});
    if file_N ~= N || file_K ~= K || ...
            abs(file_param-topology_parameter) > 5.1e-4
        continue;
    end
    candidate_path = fullfile(files(file_idx).folder, files(file_idx).name);
    saved = load(candidate_path, 'Y', 'cfg');
    if ~isfield(saved, 'Y') || isempty(saved.Y) || ...
            numel(saved.Y(end, :)) ~= 2*N*K || ~isfield(saved, 'cfg') || ...
            ~isfield(saved.cfg, 'L_intra') || ...
            ~sameNetwork(saved.cfg.L_intra, L_intra, K)
        continue;
    end
    state = saved.Y(end, :)';
    if instantaneousPatternAmplitude(state, N*K) <= A_cut
        continue;
    end
    exact_parameters = abs(file_alpha-alpha) <= 5.1e-4 && ...
        abs(file_beta-beta) <= 5.1e-4 && ...
        abs(file_noise-noise_intensity) <= 5.1e-3;
    candidate_idx = numel(candidates) + 1;
    candidates(candidate_idx).path = candidate_path;
    candidates(candidate_idx).sigma = file_sigma;
    candidates(candidate_idx).exact = exact_parameters && strcmp(token{8}, 'bwd');
    candidates(candidate_idx).state = state;
end
if isempty(candidates)
    error(['找不到与当前网络 realization 对应的已有斑图状态。' ...
        '请检查 results/evolution_* 中是否有同网络的高 sigma 斑图文件。']);
end

exact_idx = find([candidates.exact] & ([candidates.sigma] >= sigma_target));
if ~isempty(exact_idx)
    [~, order_idx] = min([candidates(exact_idx).sigma]);
    chosen = exact_idx(order_idx);
else
    [~, chosen] = max([candidates.sigma]);
end
seed_state = candidates(chosen).state;
seed_sigma = candidates(chosen).sigma;
seed_path = candidates(chosen).path;
if seed_sigma < sigma_target
    error('所选斑图种子的 sigma=%.6g 低于目标 sigma=%.6g。', ...
        seed_sigma, sigma_target);
end
end

function is_same = sameNetwork(saved_layers, current_layers, K)
is_same = false;
if ~iscell(saved_layers) || numel(saved_layers) < K
    return;
end
for layer_idx = 1:K
    if ~isequal(sparse(saved_layers{layer_idx}), sparse(current_layers{layer_idx}))
        return;
    end
end
is_same = true;
end

function [x_pattern, A_reference, n_integrations] = buildPatternReference( ...
        base_cfg, initial_seed, seed_sigma, sigma_target, A_cut, seed_base, ...
        ratio_index, max_step, min_step)
n_integrations = 0;
sigma_current = seed_sigma;
y_current = initial_seed(:);
step_size = min(max_step, max(min_step, (sigma_current-sigma_target)/20));
attempt_index = 0;

% 先在目标模型参数下重新成熟高 sigma 参考态。
attempt_index = attempt_index + 1;
[y_current, A_current] = integrateState(base_cfg, sigma_current, y_current, ...
    continuationSeed(seed_base, ratio_index, attempt_index), true);
n_integrations = n_integrations + 1;
if A_current <= A_cut
    error('r 索引 %d 的高 sigma 参考态在目标模型下未形成斑图。', ratio_index);
end

while sigma_current > sigma_target + 1e-10
    attempt_index = attempt_index + 1;
    sigma_next = max(sigma_target, sigma_current-step_size);
    [y_trial, A_trial] = integrateState(base_cfg, sigma_next, y_current, ...
        continuationSeed(seed_base, ratio_index, attempt_index), true);
    n_integrations = n_integrations + 1;
    if A_trial > A_cut
        sigma_current = sigma_next;
        y_current = y_trial;
        if A_trial > 5*A_cut
            step_size = min(max_step, 1.5*step_size);
        else
            step_size = max(min_step, 0.75*step_size);
        end
    else
        step_size = step_size / 2;
        if step_size < min_step
            error(['r 索引 %d 的 backward continuation 在 sigma=%.8g 附近' ...
                '跌破 A_cut，无法确认目标斑图分支。'], ratio_index, sigma_next);
        end
    end
    if attempt_index > 2000
        error('r 索引 %d 的 backward continuation 超过最大步数。', ratio_index);
    end
end

% 目标 sigma 下关闭提前收敛，完整运行 Figure 2 的积分时长。
% solve_multiplex 内部 A_final 即末尾 20% 时间步的序参量均值。
attempt_index = attempt_index + 1;
cfg = base_cfg;
cfg.sigma = sigma_target;
cfg.y0 = y_current;
cfg.noise_seed = continuationSeed(seed_base, ratio_index, attempt_index);
cfg.detect_convergence = false;
[~, Y, cfg_out] = solve_multiplex(cfg);
n_integrations = n_integrations + 1;
x_pattern = Y(end, :)';
A_reference = cfg_out.A_final;
if ~isfinite(A_reference) || A_reference <= A_cut
    error('r 索引 %d 的目标 sigma 斑图参考态未通过 A_cut。', ratio_index);
end
end

function [y_end, A_final] = integrateState(base_cfg, sigma, y0, noise_seed, detect)
cfg = base_cfg;
cfg.sigma = sigma;
cfg.y0 = y0;
cfg.noise_seed = noise_seed;
cfg.detect_convergence = detect;
[~, Y, cfg_out] = solve_multiplex(cfg); %#ok<ASGLU>
y_end = Y(end, :)';
A_final = instantaneousPatternAmplitude(y_end, cfg.N*cfg.K);
end

function seed = continuationSeed(seed_base, ratio_index, integration_index)
seed = mod(double(seed_base) + 1000000*ratio_index + integration_index, 2^32-1);
seed = floor(seed);
end

function amplitude = instantaneousPatternAmplitude(state, NK)
u = state(1:NK);
v = state(NK+1:2*NK);
amplitude = sqrt(sum((u-5).^2 + (v-10).^2) / NK);
end

function [lower, upper] = wilsonInterval(success_count, sample_count, z)
p_hat = success_count / sample_count;
denominator = 1 + z^2/sample_count;
center = (p_hat + z^2/(2*sample_count)) / denominator;
half_width = z/denominator * sqrt(p_hat*(1-p_hat)/sample_count + ...
    z^2/(4*sample_count^2));
lower = max(0, center-half_width);
upper = min(1, center+half_width);
end

function value = interpolateHalfProbability(c_values, probabilities)
value = NaN;
exact_idx = find(abs(probabilities-0.5) < 1e-12, 1, 'first');
if ~isempty(exact_idx)
    value = c_values(exact_idx);
    return;
end
crossing_idx = find((probabilities(1:end-1) < 0.5 & ...
    probabilities(2:end) > 0.5) | (probabilities(1:end-1) > 0.5 & ...
    probabilities(2:end) < 0.5), 1, 'first');
if isempty(crossing_idx)
    return;
end
p1 = probabilities(crossing_idx);
p2 = probabilities(crossing_idx+1);
value = c_values(crossing_idx) + (0.5-p1) * ...
    (c_values(crossing_idx+1)-c_values(crossing_idx)) / (p2-p1);
end
