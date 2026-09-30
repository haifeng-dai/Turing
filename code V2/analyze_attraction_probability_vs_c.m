%% analyze_attraction_probability_vs_c.m：沿均匀态到斑图态方向估计经验吸引概率
% 复用 Figure 2 的网络、参数、随机动力学积分器和序参量定义。
% 不生成网络或斑图参考态；严格读取目标参数已有斑图，缺失报错。
% 不改写模型或噪声构造；不保存完整随机仿真时间序列。
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
    'ALPHA_LIST', 'RATIO_LIST', 'P_VAL_FIXED');
if ~isfield(threshold_data, 'SF_Matrix') || ...
        ~isfield(threshold_data, 'SB_Matrix') || ...
        ~isfield(threshold_data, 'ALPHA_LIST')
    error('阈值 MAT 文件缺少 SF_Matrix、SB_Matrix 或 ALPHA_LIST。');
end
if ~isfield(threshold_data, 'RATIO_LIST')
    error('阈值 MAT 文件中缺少 RATIO_LIST。');
end
ratio_axis = threshold_data.RATIO_LIST(:)';

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

requested_R_LIST = [0.1, 4.5, 8.0]; % [0.1 4.5 8.0]
% 每次运行将这些 r 视为同一组实验，在共同滞回区间的中点取同一个 sigma。
% 结果按整组 r 保存；更改此列表会生成另一份整组结果。
if isempty(requested_R_LIST)
    error('requested_R_LIST 至少需要选择一个 r。');
end
FORCE_RERUN = false; % true 时重新计算当前整组实验
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
if numel(unique(R_LIST)) ~= numel(R_LIST)
    error('所选 r 映射到了重复的阈值网格点，请调整 requested_R_LIST。');
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
if sigma_common_low >= sigma_common_high
    error(['所选 r 的滞回区间没有共同区域，无法在同一个 sigma 下比较：' ...
        '\n共同下界=%.8g，共同上界=%.8g'], ...
        sigma_common_low, sigma_common_high);
end
sigma_test = (sigma_common_low + sigma_common_high) / 2;

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
SEED_BASE = 1024;
N_WORKERS_REQUESTED = 64;
INITIAL_PERTURBATION_RELATIVE_RMS = 1e-3;
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

%% 3. 整组 r 共用一个缓存文件和同一个 sigma
% 旧版各 r 独立文件采用了不同的 sigma，不属于当前共同 sigma 实验。
n_ratio = numel(R_LIST);
n_c = numel(C_LIST);
ratio_tag = strjoin(arrayfun(@ratioFileTag, R_LIST, ...
    'UniformOutput', false), '_');
result_file = fullfile(results_dir, sprintf( ...
    'attraction_probability_vs_c_common_r%s.mat', ratio_tag));
JOB_BATCH_SIZE = 256; % 每批完成后写入同一个文件，支持中断后续算

cache_config = struct( ...
    'version', 1, 'R_LIST', R_LIST, 'C_LIST', C_LIST, ...
    'N_REALIZATIONS', N_REALIZATIONS, 'sigma_test', sigma_test, ...
    'sigma_f', sigma_f, 'sigma_b', sigma_b, ...
    'topology_type', topology_type, 'topology_parameter', topology_parameter, ...
    'N', N, 'K', K, 'alpha', alpha, 'noise_intensity', noise_intensity, ...
    'T_END', T_END, 'dt', dt, 'init_perturb', init_perturb, ...
    'classification_threshold', classification_threshold, ...
    'SEED_BASE', SEED_BASE, ...
    'INITIAL_PERTURBATION_RELATIVE_RMS', INITIAL_PERTURBATION_RELATIVE_RMS);
A_final_all = NaN(n_ratio, n_c, N_REALIZATIONS);
pattern_reference_state = NaN(2*NK, n_ratio);
reference_A_final = NaN(1, n_ratio);
if exist(result_file, 'file') && ~FORCE_RERUN
    cached = load(result_file);
    if ~isfield(cached, 'cache_config') || ...
            ~isequaln(cached.cache_config, cache_config) || ...
            ~isfield(cached, 'A_final_all') || ...
            ~isequal(size(cached.A_final_all), [n_ratio, n_c, N_REALIZATIONS]) || ...
            ~isfield(cached, 'pattern_reference_state') || ...
            ~isequal(size(cached.pattern_reference_state), [2*NK, n_ratio]) || ...
            ~isfield(cached, 'reference_A_final') || ...
            ~isequal(size(cached.reference_A_final), [1, n_ratio])
        error(['整组缓存与当前实验配置不一致：%s\n' ...
            '如需重新计算，请设置 FORCE_RERUN=true；旧文件不会被自动覆盖。'], ...
            result_file);
    end
    A_final_all = cached.A_final_all;
    pattern_reference_state = cached.pattern_reference_state;
    reference_A_final = cached.reference_A_final;
    fprintf('[CACHE] 已读取整组结果：%s，完成 %d/%d。\n', ...
        result_file, sum(isfinite(A_final_all(:))), numel(A_final_all));
end

timer_total = tic;
pool = [];
actual_workers = 0;
n_jobs_run = 0;
n_continuation_integrations = 0;

%% 4. 为需要补算的 r 读取目标斑图参考态，绝不自动训练或延拓
for r_idx = 1:n_ratio
    A_r = reshape(A_final_all(r_idx, :, :), [n_c, N_REALIZATIONS]);
    if all(isfinite(A_r(:)))
        continue;
    end
    if all(isfinite(pattern_reference_state(:, r_idx))) && ...
            isfinite(reference_A_final(r_idx)) && ...
            reference_A_final(r_idx) > classification_threshold
        continue;
    end
    fprintf('[REFERENCE] 正在读取 r=%.6g、sigma=%.8g 的已有目标斑图...\n', ...
        R_LIST(r_idx), sigma_test);
    reference_timer = tic;
    cfg_r = base_cfg;
    cfg_r.beta = alpha * R_LIST(r_idx);
    cfg_r.sigma = sigma_test;
    [reference_state_r, reference_A_r] = load_existing_pattern_reference( ...
        results_dir, topology_type, topology_parameter, cfg_r, classification_threshold);
    pattern_reference_state(:, r_idx) = reference_state_r;
    reference_A_final(r_idx) = reference_A_r;
    if reference_A_r <= 5 * classification_threshold
        warning('r=%.3f 的参考斑图振幅仅为 A_cut 的 %.2f 倍，请检查状态分离。', ...
            R_LIST(r_idx), reference_A_r / classification_threshold);
    end
    save(result_file, 'cache_config', 'A_final_all', ...
        'pattern_reference_state', 'reference_A_final', '-v7');
    fprintf('[REFERENCE] r=%.6g 读取完成，新增积分 0 次，耗时 %.1f s。\n', ...
        R_LIST(r_idx), toc(reference_timer));
end

%% 5. 所有 (r, c, realization) 一起并行计算，按批保存整组缓存
pending_linear = find(~isfinite(A_final_all(:)));
if ~isempty(pending_linear)
    fprintf('[POOL] 正在启动并行池，准备 %d 个待计算任务。\n', ...
        numel(pending_linear));
    pool = startLocalPool(N_WORKERS_REQUESTED);
    actual_workers = pool.NumWorkers;
    fprintf('Requested workers: %d\nActual workers: %d\n', ...
        N_WORKERS_REQUESTED, actual_workers);
    fprintf('[SCAN] 共同 sigma=%.8g，待计算 %d/%d 个随机演化。\n', ...
        sigma_test, numel(pending_linear), numel(A_final_all));
    for batch_start = 1:JOB_BATCH_SIZE:numel(pending_linear)
        batch_end = min(batch_start+JOB_BATCH_SIZE-1, numel(pending_linear));
        batch_linear = pending_linear(batch_start:batch_end);
        [batch_r, batch_c, batch_realization] = ind2sub( ...
            [n_ratio, n_c, N_REALIZATIONS], batch_linear);
        batch_A = NaN(size(batch_linear));
        parfor job_idx = 1:numel(batch_linear)
            r_idx = batch_r(job_idx);
            c_idx = batch_c(job_idx);
            realization_idx = batch_realization(job_idx);
            x_P = pattern_reference_state(:, r_idx);
            y0 = x_H + C_LIST(c_idx) * (x_P - x_H) + ...
                epsilon_initial * zeta;
            cfg = base_cfg;
            cfg.beta = alpha * R_LIST(r_idx);
            cfg.sigma = sigma_test;
            cfg.y0 = y0;
            cfg.noise_seed = SEED_BASE + realization_idx;
            cfg.detect_convergence = false;
            [~, ~, cfg_out] = solve_multiplex(cfg);
            batch_A(job_idx) = cfg_out.A_final;
        end
        if any(~isfinite(batch_A))
            error('当前批次有非有限的 A_final，已完成的前序批次仍保存在 %s。', ...
                result_file);
        end
        A_final_all(batch_linear) = batch_A;
        n_jobs_run = n_jobs_run + numel(batch_linear);
        save(result_file, 'cache_config', 'A_final_all', ...
            'pattern_reference_state', 'reference_A_final', '-v7');
        fprintf('[SCAN] 已完成 %d/%d，整组结果已保存。\n', ...
            sum(isfinite(A_final_all(:))), numel(A_final_all));
    end
else
    fprintf('整组缓存已完成，无需补算，也未启动并行池。\n');
end

if any(~isfinite(A_final_all(:)))
    error('整组结果尚未完成，无法绘制完整比较图：%s', result_file);
end

P_pattern = NaN(n_ratio, n_c);
CI_low = NaN(n_ratio, n_c);
CI_high = NaN(n_ratio, n_c);
c50 = NaN(1, n_ratio);
for r_idx = 1:n_ratio
    A_r = reshape(A_final_all(r_idx, :, :), [n_c, N_REALIZATIONS]);
    n_pattern_r = sum(A_r > classification_threshold, 2)';
    P_pattern(r_idx, :) = n_pattern_r / N_REALIZATIONS;
    for c_idx = 1:n_c
        [CI_low(r_idx, c_idx), CI_high(r_idx, c_idx)] = wilsonInterval( ...
            n_pattern_r(c_idx), N_REALIZATIONS, Z_WILSON);
    end
    c50(r_idx) = interpolateHalfProbability(C_LIST, P_pattern(r_idx, :));
end
save(result_file, 'cache_config', 'A_final_all', ...
    'pattern_reference_state', 'reference_A_final', ...
    'P_pattern', 'CI_low', 'CI_high', 'c50', '-v7');
total_wall_time = toc(timer_total);
n_total_simulations = n_jobs_run + n_continuation_integrations;

%% 6. 绘制经验吸引概率曲线
curve_colors = lines(n_ratio);
main_figure = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [100 100 900 620]);
main_axes = axes(main_figure);
hold(main_axes, 'on');
has_curve = false;
for r_idx = 1:n_ratio
    if any(~isfinite(P_pattern(r_idx, :)))
        continue;
    end
    lower_error = P_pattern(r_idx, :) - CI_low(r_idx, :);
    upper_error = CI_high(r_idx, :) - P_pattern(r_idx, :);
    errorbar(main_axes, C_LIST, P_pattern(r_idx, :), lower_error, upper_error, ...
        '-o', 'Color', curve_colors(r_idx, :), ...
        'MarkerFaceColor', curve_colors(r_idx, :), 'MarkerSize', 5, ...
        'LineWidth', 1.5, 'CapSize', 5, ...
        'DisplayName', sprintf('r = %.3g', R_LIST(r_idx)));
    has_curve = true;
end
yline(main_axes, 0.5, '--', 'Color', [0.55 0.55 0.55], ...
    'HandleVisibility', 'off');
xlim(main_axes, [0 1]);
ylim(main_axes, [0 1]);
xlabel(main_axes, 'Initial-condition interpolation parameter c');
ylabel(main_axes, 'Probability of patterned final state');
title(main_axes, 'Empirical attraction probability along a representative direction');
subtitle(main_axes, sprintf('Common $\\sigma = %.6g$ for all selected $r$', ...
    sigma_test), 'Interpreter', 'latex');
if has_curve
    legend(main_axes, 'Location', 'best');
end
grid(main_axes, 'on');
box(main_axes, 'on');
main_image = fullfile(fig_dir, 'attraction_probability_vs_c.png');
exportgraphics(main_figure, main_image, 'Resolution', 300);
close(main_figure);

%% 7. 输出运行摘要
fprintf('\n=== Attraction probability analysis complete ===\n');
fprintf('Figure 2 parameters: alpha=%.6g, eta=%.6g, N=%d, K=%d\n', ...
    alpha, noise_intensity, N, K);
fprintf('Network file: %s\n', topology_file);
fprintf('Threshold data: %s\n', threshold_file);
fprintf('Actual r values: %s\n', mat2str(R_LIST, 6));
for r_idx = 1:n_ratio
    fprintf('r=%.6g: sigma_b=%.8g, sigma_f=%.8g, sigma_test=%.8g, c50=', ...
        R_LIST(r_idx), sigma_b(r_idx), sigma_f(r_idx), sigma_test);
    if isfinite(c50(r_idx))
        fprintf('%.6g\n', c50(r_idx));
    else
        fprintf('未定义（曲线未穿过 0.5）\n');
    end
end
fprintf('Sigma selection: midpoint of common hysteresis interval %.8g\n', ...
    sigma_test);
fprintf('Combined result file: %s\n', result_file);
fprintf('Requested workers: %d; actual workers: %d\n', ...
    N_WORKERS_REQUESTED, actual_workers);
fprintf(['Total simulations: %d (Monte Carlo: %d; continuation integrations: %d); ' ...
    'elapsed: %.2f s\n'], n_total_simulations, n_jobs_run, ...
    n_continuation_integrations, total_wall_time);
fprintf('Figure: %s\n', main_image);

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

function tag = ratioFileTag(ratio)
% 将 beta/alpha 转为稳定的文件名片段，例如 4.5 -> 4p5。
tag = sprintf('%.10g', ratio);
tag = strrep(tag, '-', 'm');
tag = strrep(tag, '.', 'p');
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
