%% compare_forward_threshold_spectral_reference.m
% 比较 Fig. 2 随机仿真 forward 阈值与确定性超拉普拉斯谱参考阈值。
%
% 这里的谱分析只使用均匀稳态附近的确定性线性化，不包含乘性噪声 eta。
% 因此 sigma_c_spec 只能称为 deterministic linear spectral reference，
% 不能称为随机系统的精确理论阈值，也不用于预测 backward threshold。
clear; clc; close all;

%% 1. Fig. 2 参数与输出配置
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(script_dir, 'simulations'), fullfile(script_dir, 'networks'));

TOPO_TYPE = 'ER';
N = 200;
K = 5;
P_VAL_FIXED = 0.030;
ALPHA_FIXED = 0.1;
ETA_REFERENCE = 0.01;  % Fig. 2 仿真配置中的噪声参数；谱计算不使用它

ALPHA_TOL = 1e-10;
NONZERO_TOL = 1e-10;
ROOT_TOL = 1e-6;
COARSE_STEP = 0.1;
SIGMA_INITIAL_MAX = 100;
SIGMA_MAX_LIMIT = 2000;

res_dir = fullfile(script_dir, 'results');
fig_dir = fullfile(script_dir, 'fig');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Fig. 2 的 beta/alpha 扫描结果由 sigma_beta_alpha.m 写入。
threshold_file = fullfile(res_dir, sprintf( ...
    'scan_beta_alpha_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), N, K, P_VAL_FIXED));
topology_file = fullfile(res_dir, 'topology', upper(TOPO_TYPE), ...
    sprintf('N%d_p%.3f.mat', N, P_VAL_FIXED));

if ~exist(threshold_file, 'file')
    error('找不到 Fig. 2 forward 阈值文件：%s', threshold_file);
end
if ~exist(topology_file, 'file')
    error('找不到 Fig. 2 网络文件：%s', topology_file);
end

%% 2. 读取 Fig. 2 forward 阈值和 beta/alpha 横坐标
threshold_data = load(threshold_file);
required_threshold_fields = {'SF_Matrix', 'ALPHA_LIST'};
for field_idx = 1:numel(required_threshold_fields)
    field_name = required_threshold_fields{field_idx};
    if ~isfield(threshold_data, field_name)
        error('Fig. 2 阈值文件缺少字段：%s', field_name);
    end
end

if isfield(threshold_data, 'RATIO_LIST')
    RATIO_LIST = threshold_data.RATIO_LIST(:)';
elseif isfield(threshold_data, 'BETA_LIST')
    if ~isfield(threshold_data, 'ALPHA_LIST')
        error('阈值文件只有 BETA_LIST，但缺少 ALPHA_LIST，无法恢复 beta/alpha。');
    end
    RATIO_LIST = threshold_data.BETA_LIST(:)' / ALPHA_FIXED;
else
    error('Fig. 2 阈值文件缺少 RATIO_LIST 或 BETA_LIST。');
end

alpha_list = threshold_data.ALPHA_LIST(:)';
[alpha_error, alpha_idx] = min(abs(alpha_list - ALPHA_FIXED));
if alpha_error > ALPHA_TOL
    error('Fig. 2 结果中没有 alpha=%.8g，最近值为 %.8g。', ...
        ALPHA_FIXED, alpha_list(alpha_idx));
end

sigma_c_sim = threshold_data.SF_Matrix(alpha_idx, :);
sigma_c_sim = sigma_c_sim(:)';
if numel(sigma_c_sim) ~= numel(RATIO_LIST)
    error('RATIO_LIST 与 SF_Matrix 的列数不一致。');
end

if isfield(threshold_data, 'P_VAL_FIXED')
    if abs(threshold_data.P_VAL_FIXED - P_VAL_FIXED) > 1e-12
        error('阈值文件的 p=%.8g，与当前配置 p=%.8g 不一致。', ...
            threshold_data.P_VAL_FIXED, P_VAL_FIXED);
    end
end

fprintf('Fig. 2 forward 阈值文件：%s\n', threshold_file);
fprintf('Fig. 2 网络 realization：%s\n', topology_file);
fprintf('参数：N=%d，K=%d，alpha=%.8g，eta=%.8g，p=%.8g\n', ...
    N, K, ALPHA_FIXED, ETA_REFERENCE, P_VAL_FIXED);
fprintf('读取到 %d 个 beta/alpha 点，alpha 行索引为 %d。\n', ...
    numel(RATIO_LIST), alpha_idx);

%% 3. 读取与 Fig. 2 完全相同的 ER 网络 realization
topology_data = load(topology_file, 'nets');
if ~isfield(topology_data, 'nets') || numel(topology_data.nets) < K
    error('网络文件中的层数不足：需要 %d 层。', K);
end

L_intra = cell(K, 1);
for layer_idx = 1:K
    Lk = topology_data.nets{layer_idx};
    % 当前项目的 gen_er/gen_ws/gen_ba/gen_sf_gamma 保存的就是 D-A。
    % 若误传邻接矩阵，行和通常不为零，这里直接报错，避免静默得到错误谱。
    if max(abs(full(sum(Lk, 2)))) > 1e-8 || any(diag(Lk) < -1e-12)
        error('nets{%d} 不是拉普拉斯矩阵。若源文件保存邻接矩阵，应先构造 D-A。', layer_idx);
    end
    if norm(full(Lk - Lk'), 'fro') > 1e-8
        error('nets{%d} 不是对称拉普拉斯矩阵，无法使用当前实谱分析。', layer_idx);
    end
    L_intra{layer_idx} = sparse(Lk);
end

adj_inter = ones(K) - eye(K);
L_inter = diag(sum(adj_inter, 2)) - adj_inter;
L_intra_all = blkdiag(L_intra{:});
L_inter_big = kron(L_inter, speye(N));

%% 4. 计算每个 beta/alpha 的谱阈值
% Mimura-Murray 稳态 (u*,v*)=(5,10) 处的 Jacobian。
J = [10 / 3, -5; 10, -4];
lambda_all = cell(size(RATIO_LIST));
lambda_nonzero_all = cell(size(RATIO_LIST));
sigma_c_spec = NaN(size(RATIO_LIST));
lambda_star = NaN(size(RATIO_LIST));
q_star = NaN(size(RATIO_LIST));
mu_star = NaN(size(RATIO_LIST));
mu_at_sim_threshold = NaN(size(RATIO_LIST));
phi_star = cell(size(RATIO_LIST));

for ratio_idx = 1:numel(RATIO_LIST)
    ratio = RATIO_LIST(ratio_idx);

    % 论文记号：L = L^L + (beta/alpha) L^I。
    L_super = L_intra_all + ratio * L_inter_big;
    [eigenvectors, eigenvalues_matrix] = eig(full(L_super));
    [lambda_sorted, order] = sort(real(diag(eigenvalues_matrix)), 'ascend');
    eigenvectors = eigenvectors(:, order);

    if lambda_sorted(1) < -1e-8 || any(lambda_sorted < -1e-8)
        error('ratio=%.8g 的 supra-Laplacian 出现明显负特征值。', ratio);
    end
    if abs(lambda_sorted(1)) > 1e-8
        warning('ratio=%.8g 的最小特征值不是数值零：lambda_1=%.8g。', ...
            ratio, lambda_sorted(1));
    end

    nonzero_mask = lambda_sorted > NONZERO_TOL;
    lambda_nonzero = lambda_sorted(nonzero_mask);
    if isempty(lambda_nonzero)
        warning('ratio=%.8g 没有非零空间模态，谱阈值设为 NaN。', ratio);
        lambda_all{ratio_idx} = lambda_sorted;
        lambda_nonzero_all{ratio_idx} = lambda_nonzero;
        continue;
    end

    lambda_all{ratio_idx} = lambda_sorted;
    lambda_nonzero_all{ratio_idx} = lambda_nonzero;

    finite_sim_thresholds = sigma_c_sim(isfinite(sigma_c_sim));
    if isempty(finite_sim_thresholds)
        sigma_upper = SIGMA_INITIAL_MAX;
    else
        sigma_upper = max(SIGMA_INITIAL_MAX, max(finite_sim_thresholds) * 1.25);
    end

    [sigma_root, mu_root, root_found] = find_spectral_threshold( ...
        lambda_nonzero, ALPHA_FIXED, J, 0, sigma_upper, ...
        SIGMA_MAX_LIMIT, COARSE_STEP, ROOT_TOL);

    if ~root_found
        warning(['ratio=%.8g 在 sigma <= %.8g 范围内没有从负到正的线性谱失稳，' ...
            'sigma_c_spec 设为 NaN。'], ratio, SIGMA_MAX_LIMIT);
        continue;
    end

    sigma_c_spec(ratio_idx) = sigma_root;
    mu_star(ratio_idx) = mu_root;

    [~, q_local, mu_modes] = spectral_growth( ...
        sigma_root, lambda_nonzero, ALPHA_FIXED, J);
    lambda_star(ratio_idx) = lambda_nonzero(q_local);
    nonzero_indices = find(nonzero_mask);
    q_star(ratio_idx) = nonzero_indices(q_local);
    phi_star{ratio_idx} = eigenvectors(:, q_star(ratio_idx));

    if isfinite(sigma_c_sim(ratio_idx))
        [mu_at_sim_threshold(ratio_idx), ~, ~] = spectral_growth( ...
            sigma_c_sim(ratio_idx), lambda_nonzero, ALPHA_FIXED, J);
    end

    % 记录根附近的符号检查，确认是从稳定侧穿过 mu=0。
    delta_sigma = max(1e-4, 10 * ROOT_TOL);
    [mu_below, ~, ~] = spectral_growth( ...
        max(0, sigma_root - delta_sigma), lambda_nonzero, ALPHA_FIXED, J);
    [mu_above, ~, ~] = spectral_growth( ...
        sigma_root + delta_sigma, lambda_nonzero, ALPHA_FIXED, J);
    if ~(mu_below < 0 && mu_above > 0)
        warning(['ratio=%.8g 的根附近符号检查未满足 mu_below<0<mu_above：' ...
            'below=%.8g，root=%.8g，above=%.8g。'], ...
            ratio, mu_below, mu_root, mu_above);
    end
end

%% 5. 汇总结果并保存 MAT、CSV
difference = sigma_c_sim - sigma_c_spec;
relative_difference = difference ./ sigma_c_sim;

threshold_table = table(RATIO_LIST(:), sigma_c_sim(:), sigma_c_spec(:), ...
    difference(:), relative_difference(:), lambda_star(:), q_star(:), ...
    mu_star(:), mu_at_sim_threshold(:), ...
    'VariableNames', {'ratio', 'sigma_c_sim', 'sigma_c_spec', 'difference', ...
    'relative_difference', 'lambda_star', 'q_star', 'mu_star', ...
    'mu_at_sim_threshold'});

base_name = sprintf('compare_forward_threshold_spectral_%s_N%d_K%d_p%.3f_a%.3f', ...
    lower(TOPO_TYPE), N, K, P_VAL_FIXED, ALPHA_FIXED);
out_mat = fullfile(res_dir, [base_name, '.mat']);
out_csv = fullfile(res_dir, [base_name, '.csv']);
save(out_mat, 'threshold_table', 'RATIO_LIST', 'sigma_c_sim', 'sigma_c_spec', ...
    'difference', 'relative_difference', 'lambda_all', 'lambda_nonzero_all', ...
    'lambda_star', 'q_star', 'mu_star', 'mu_at_sim_threshold', 'phi_star', ...
    'J', 'ALPHA_FIXED', 'ETA_REFERENCE', 'P_VAL_FIXED', ...
    'NONZERO_TOL', 'ROOT_TOL', 'threshold_file', 'topology_file', '-v7');
writetable(threshold_table, out_csv);

%% 6. 主图：simulation forward threshold 与 spectral reference
main_fig = figure('Visible', 'off', 'Color', 'w', ...
    'Name', 'Forward threshold comparison', ...
    'Units', 'normalized', 'Position', [0.18, 0.18, 0.58, 0.62]);
main_ax = axes(main_fig);
hold(main_ax, 'on');
plot(main_ax, RATIO_LIST, sigma_c_sim, 'o-', ...
    'Color', [0.00, 0.447, 0.741], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.00, 0.447, 0.741], 'MarkerSize', 5, ...
    'DisplayName', 'Simulation forward threshold');
plot(main_ax, RATIO_LIST, sigma_c_spec, 's--', ...
    'Color', [0.850, 0.325, 0.098], 'LineWidth', 1.8, ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
    'DisplayName', 'Deterministic spectral reference');
xlabel(main_ax, '\beta/\alpha', 'FontSize', 15, 'Interpreter', 'tex');
ylabel(main_ax, '\sigma_c', 'FontSize', 15, 'Interpreter', 'tex');
title(main_ax, sprintf('Forward threshold comparison (ER, p=%.3f, N=%d, K=%d, \\alpha=%.3g)', ...
    P_VAL_FIXED, N, K, ALPHA_FIXED), 'FontSize', 14);
legend(main_ax, 'Location', 'best', 'FontSize', 10);
grid(main_ax, 'on');
box(main_ax, 'on');
set(main_ax, 'FontSize', 12);
main_img = fullfile(fig_dir, [base_name, '.png']);
exportgraphics(main_fig, main_img, 'Resolution', 300);

%% 7. 诊断图：lambda_star 与仿真阈值处的 mu_max
diagnostic_fig = figure('Visible', 'off', 'Color', 'w', ...
    'Name', 'Spectral threshold diagnostics', ...
    'Units', 'normalized', 'Position', [0.18, 0.14, 0.58, 0.72]);
diagnostic_layout = tiledlayout(diagnostic_fig, 2, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

ax_lambda = nexttile(diagnostic_layout);
plot(ax_lambda, RATIO_LIST, lambda_star, 'o-', ...
    'Color', [0.494, 0.184, 0.556], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.494, 0.184, 0.556], 'MarkerSize', 5);
xlabel(ax_lambda, '\beta/\alpha', 'FontSize', 14, 'Interpreter', 'tex');
ylabel(ax_lambda, '\lambda_*', 'FontSize', 14, 'Interpreter', 'tex');
title(ax_lambda, '最危险超拉普拉斯模态特征值', 'FontSize', 13);
grid(ax_lambda, 'on');
box(ax_lambda, 'on');
set(ax_lambda, 'FontSize', 12);

ax_mu = nexttile(diagnostic_layout);
plot(ax_mu, RATIO_LIST, mu_at_sim_threshold, 'o-', ...
    'Color', [0.466, 0.674, 0.188], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.466, 0.674, 0.188], 'MarkerSize', 5);
yline(ax_mu, 0, 'k--', 'LineWidth', 0.9);
xlabel(ax_mu, '\beta/\alpha', 'FontSize', 14, 'Interpreter', 'tex');
ylabel(ax_mu, '\mu_{max}(\\sigma_c^{sim})', 'FontSize', 14, 'Interpreter', 'tex');
title(ax_mu, '仿真 forward 阈值处的确定性谱增长率', 'FontSize', 13);
grid(ax_mu, 'on');
box(ax_mu, 'on');
set(ax_mu, 'FontSize', 12);

diagnostic_img = fullfile(fig_dir, [base_name, '_diagnostics.png']);
exportgraphics(diagnostic_fig, diagnostic_img, 'Resolution', 300);

%% 8. 一致性检查与结果摘要
lambda_one = zeros(size(RATIO_LIST));
for ratio_idx = 1:numel(RATIO_LIST)
    lambda_one(ratio_idx) = lambda_all{ratio_idx}(1);
end
fprintf('\n[检查] max(abs(lambda_1)) = %.8g\n', max(abs(lambda_one)));
fprintf('[输出] 主图：%s\n', main_img);
fprintf('[输出] 诊断图：%s\n', diagnostic_img);
fprintf('[输出] MAT：%s\n', out_mat);
fprintf('[输出] CSV：%s\n', out_csv);
fprintf('\n ratio       sigma_sim       sigma_spec      difference      lambda_star       mu_at_sim\n');
for ratio_idx = 1:numel(RATIO_LIST)
    fprintf('%8.4g  %13.8g  %13.8g  %13.8g  %13.8g  %13.8g\n', ...
        threshold_table.ratio(ratio_idx), threshold_table.sigma_c_sim(ratio_idx), ...
        threshold_table.sigma_c_spec(ratio_idx), threshold_table.difference(ratio_idx), ...
        threshold_table.lambda_star(ratio_idx), threshold_table.mu_at_sim_threshold(ratio_idx));
end

finite_diff_idx = find(isfinite(difference));
if ~isempty(finite_diff_idx)
    [max_abs_diff, local_diff_idx] = max(abs(difference(finite_diff_idx)));
    max_diff_idx = finite_diff_idx(local_diff_idx);
    fprintf('\n最大绝对差值：%.8g，ratio=%.8g。\n', ...
        max_abs_diff, RATIO_LIST(max_diff_idx));
    fprintf('谱参考与仿真 forward 阈值的整体趋势需结合上述曲线和差值判断。\n');
end

%% 局部函数
function [mu_max, q_star, mu_modes] = spectral_growth(sigma, lambda_nonzero, alpha, J)
% 计算所有非零空间模态的增长率，并返回最危险模态。
a = J(1, 1) - alpha * lambda_nonzero;
b = J(1, 2);
c = J(2, 1);
d = J(2, 2) - alpha * sigma * lambda_nonzero;
trace_modes = a + d;
discriminant = (a - d).^2 + 4 * b * c;
root_discriminant = sqrt(complex(discriminant, 0));
eigenvalue_plus = (trace_modes + root_discriminant) / 2;
eigenvalue_minus = (trace_modes - root_discriminant) / 2;
mu_modes = max(real(eigenvalue_plus), real(eigenvalue_minus));
[mu_max, q_star] = max(mu_modes);
end

function [sigma_root, mu_root, root_found] = find_spectral_threshold( ...
    lambda_nonzero, alpha, J, sigma_min, sigma_initial_max, ...
    sigma_max_limit, coarse_step, root_tol)
% 粗扫描寻找最早的负到正区间，再用二分法精确求根。
sigma_root = NaN;
mu_root = NaN;
root_found = false;
sigma_upper = max(sigma_initial_max, sigma_min + coarse_step);

while true
    sigma_grid = sigma_min:coarse_step:sigma_upper;
    mu_grid = zeros(size(sigma_grid));
    for sigma_idx = 1:numel(sigma_grid)
        [mu_grid(sigma_idx), ~, ~] = spectral_growth( ...
            sigma_grid(sigma_idx), lambda_nonzero, alpha, J);
    end

    crossing_idx = find(mu_grid(1:end-1) <= 0 & mu_grid(2:end) >= 0, 1, 'first');
    if ~isempty(crossing_idx)
        lo = sigma_grid(crossing_idx);
        hi = sigma_grid(crossing_idx + 1);
        mu_lo = mu_grid(crossing_idx);

        for iteration = 1:100
            mid = 0.5 * (lo + hi);
            [mu_mid, ~, ~] = spectral_growth(mid, lambda_nonzero, alpha, J);
            if mu_mid >= 0
                hi = mid;
            else
                lo = mid;
                mu_lo = mu_mid;
            end
            if hi - lo <= root_tol
                break;
            end
        end

        sigma_root = 0.5 * (lo + hi);
        [mu_root, ~, ~] = spectral_growth(sigma_root, lambda_nonzero, alpha, J);
        root_found = isfinite(mu_root) && mu_lo <= 0;
        return;
    end

    if sigma_upper >= sigma_max_limit
        return;
    end
    sigma_upper = min(2 * sigma_upper, sigma_max_limit);
end
end
