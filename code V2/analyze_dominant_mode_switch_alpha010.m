%% 诊断 beta/alpha=4.3--4.6 附近的主导超拉普拉斯模态
% 本脚本只读取已有谱分析结果和 ER 网络实例；不生成网络，
% 也不运行随机仿真或动力学仿真。简并检查会对已有网络做确定性特征分解。
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
res_dir = fullfile(script_dir, 'results');
fig_dir = fullfile(script_dir, 'fig');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

input_file = fullfile(res_dir, ...
    'compare_forward_threshold_spectral_er_N200_K5_p0.030_a0.100.mat');
if ~exist(input_file, 'file')
    error('找不到已有的谱分析结果文件：%s', input_file);
end
S = load(input_file);

RATIO_LIST = getField(S, {'RATIO_LIST','ratio_list'});
lambda_star = getField(S, {'lambda_star','LAMBDA_STAR'});
q_star = getField(S, {'q_star','Q_STAR'});
phi_star = getField(S, {'phi_star','PHI_STAR'});
lambda_all = getField(S, {'lambda_all','LAMBDA_ALL'});
sigma_c_spec = getField(S, {'sigma_c_spec','SIGMA_C_SPEC'});
sigma_c_sim = getField(S, {'sigma_c_sim','SIGMA_C_SIM'});

RATIO_LIST = RATIO_LIST(:)';
lambda_star = lambda_star(:)';
q_star = q_star(:)';
sigma_c_spec = sigma_c_spec(:)';
sigma_c_sim = sigma_c_sim(:)';
if isfield(S, 'N'), N = S.N; else, N = 200; end
if isfield(S, 'K'), K = S.K; else, K = 5; end
if isfield(S, 'ALPHA_FIXED'), alpha = S.ALPHA_FIXED; else, alpha = 0.1; end
if isfield(S, 'P_VAL_FIXED'), p_val = S.P_VAL_FIXED; else, p_val = 0.030; end
if N ~= 200 || K ~= 5 || abs(alpha - 0.1) > 1e-12
    error('本诊断脚本要求 N=200、K=5、alpha=0.1，但 MAT 文件元数据不匹配。');
end
if any([numel(lambda_star), numel(q_star), numel(sigma_c_spec), ...
        numel(sigma_c_sim), numel(lambda_all), numel(phi_star)] ~= numel(RATIO_LIST))
    error('谱结果字段长度与 RATIO_LIST 不一致。');
end
if any(diff(RATIO_LIST) < 0)
    error('RATIO_LIST 必须按升序排列，才能判断相邻参数点之间的跳变。');
end

% 在连续色散关系的有效开区间内数值寻找最小 sigma。
lambda_upper = (10/3) / alpha;
lambda_lo = max(realmin, lambda_upper * 1e-9);
lambda_hi = lambda_upper * (1 - 1e-9);
sigma_of_lambda = @(lambda) (110/3 + 4*alpha.*lambda) ./ ...
    ((alpha.*lambda) .* (10/3 - alpha.*lambda));
lambda_grid = linspace(lambda_lo, lambda_hi, 100001);
sigma_grid = sigma_of_lambda(lambda_grid);
[~, grid_min_idx] = min(sigma_grid);
bracket_lo = lambda_grid(max(1, grid_min_idx - 1));
bracket_hi = lambda_grid(min(numel(lambda_grid), grid_min_idx + 1));
lambda_opt = fminbnd(sigma_of_lambda, bracket_lo, bracket_hi);
sigma_min = sigma_of_lambda(lambda_opt);
fprintf('连续图灵色散关系最优点：lambda_opt=%.12g，sigma_min=%.12g（alpha=%.6g）\n', ...
    lambda_opt, sigma_min, alpha);

selected_ratios = [4.3 4.4 4.5 4.6];
n_selected = numel(selected_ratios);
selected_indices = zeros(1, n_selected);
selected_ratio_actual = zeros(1, n_selected);
selected_lambda = zeros(1, n_selected);
selected_q = zeros(1, n_selected);
selected_sigma_spec = zeros(1, n_selected);
selected_sigma_sim = zeros(1, n_selected);
gap_left = NaN(1, n_selected);
gap_right = NaN(1, n_selected);
degeneracy_tol = NaN(1, n_selected);
degeneracy_flags = false(1, n_selected);
degeneracy_multiplicity = ones(1, n_selected);
selected_phi = cell(1, n_selected);
selected_mode_components = cell(1, n_selected);
selected_eigenspace_energy = cell(1, n_selected);
selected_heatmap = cell(1, n_selected);
layer_energy = zeros(K, n_selected);
layer_average_profile = cell(1, n_selected);
layer_norm_profile = cell(1, n_selected);
synchronous_fraction = NaN(1, n_selected);
transverse_fraction = NaN(1, n_selected);
effective_nodes = NaN(1, n_selected);
node_participation_fraction = NaN(1, n_selected);

% 使用此前生成的同一份网络实例，不重新生成网络。
topology_file = fullfile(res_dir, 'topology', 'ER', ...
    sprintf('N%d_p%.3f.mat', N, p_val));
if isfield(S, 'topology_file')
    saved_topology_file = char(S.topology_file);
    if exist(saved_topology_file, 'file')
        topology_file = saved_topology_file;
    elseif exist(fullfile(script_dir, saved_topology_file), 'file')
        topology_file = fullfile(script_dir, saved_topology_file);
    end
end
if ~exist(topology_file, 'file')
    error('找不到检查简并子空间所需的已有网络实例：%s', topology_file);
end
T = load(topology_file, 'nets');
if ~isfield(T, 'nets') || numel(T.nets) < K
    error('网络文件中的层数不足，需要 %d 层：%s', K, topology_file);
end
L_intra = cell(K, 1);
for k = 1:K
    L_intra{k} = sparse(T.nets{k});
    if ~isequal(size(L_intra{k}), [N N])
        error('第 %d 层网络拉普拉斯矩阵的尺寸不是 N×N。', k);
    end
end
L_intra_all = blkdiag(L_intra{:});
L_inter = diag(sum(ones(K) - eye(K), 2)) - (ones(K) - eye(K));
L_inter_big = kron(L_inter, speye(N));

for j = 1:n_selected
    [ratio_error, idx] = min(abs(RATIO_LIST - selected_ratios(j)));
    if ratio_error > 1e-7
        error('已有 RATIO_LIST 中没有 %.1f 这个点（最近点相差 %.4g）。', ...
            selected_ratios(j), ratio_error);
    end
    selected_indices(j) = idx;
    selected_ratio_actual(j) = RATIO_LIST(idx);
    selected_lambda(j) = lambda_star(idx);
    selected_q(j) = q_star(idx);
    selected_sigma_spec(j) = sigma_c_spec(idx);
    selected_sigma_sim(j) = sigma_c_sim(idx);
    if ~isfinite(lambda_star(idx)) || ~isfinite(q_star(idx))
    error('beta/alpha=%.3f 处的 lambda_star 或 q_star 不是有限值。', RATIO_LIST(idx));
    end

    spectrum = lambda_all{idx};
    spectrum = real(spectrum(:));
    if any(diff(spectrum) < -1e-10)
        error('beta/alpha=%.3f 处 lambda_all 不是升序排列。', RATIO_LIST(idx));
    end
    if abs(q_star(idx) - round(q_star(idx))) > 1e-8
        error('beta/alpha=%.3f 处 q_star 不是整数索引。', RATIO_LIST(idx));
    end
    q = round(q_star(idx));
    if q < 1 || q > numel(spectrum)
        error('beta/alpha=%.3f 处 q_star=%.8g 超出 lambda_all 的索引范围。', ...
            RATIO_LIST(idx), q_star(idx));
    end
    % 上游脚本把 q_star 定义为完整升序超拉普拉斯谱中的索引。
    if abs(spectrum(q) - lambda_star(idx)) > 1e-6 * max(1, abs(lambda_star(idx)))
        error('beta/alpha=%.3f 处 lambda_all(q_star) 与 lambda_star 不一致。', RATIO_LIST(idx));
    end
    if q > 1, gap_left(j) = spectrum(q) - spectrum(q-1); end
    if q < numel(spectrum), gap_right(j) = spectrum(q+1) - spectrum(q); end
    degeneracy_tol(j) = max(1e-8, 1e-6 * max(1, abs(lambda_star(idx))));
    near_indices = find(abs(spectrum - lambda_star(idx)) < degeneracy_tol(j));
    local_gaps = [gap_left(j), gap_right(j)];
    degeneracy_flags(j) = numel(near_indices) > 1 || ...
        any(local_gaps(isfinite(local_gaps)) < 1e-6);
    degeneracy_multiplicity(j) = numel(near_indices);

    if degeneracy_flags(j)
        warning(['dominant eigenvalue is degenerate or nearly degenerate; ' ...
            'an individual eigenvector is not uniquely defined.']);
        % 源 MAT 只保存 phi_star，没有保存完整特征向量；因此仅在近简并点，
        % 使用已有网络实例重新做确定性特征分解，以取出整个子空间。
        L_super = L_intra_all + RATIO_LIST(idx) * L_inter_big;
        [Vfull, Dfull] = eig(full(L_super));
        [evals_full, order] = sort(real(diag(Dfull)), 'ascend');
        Vfull = Vfull(:, order);
        eig_match_error = min(abs(evals_full - lambda_star(idx)));
        if eig_match_error > 1e-7 * max(1, abs(lambda_star(idx)))
            error(['beta/alpha=%.3f 处从当前网络实例重算的谱与已保存 lambda_star 不匹配；' ...
                '可能网络 realization 已改变。'], RATIO_LIST(idx));
        end
        cluster = abs(evals_full - lambda_star(idx)) < degeneracy_tol(j);
        if ~any(cluster)
            [~, nearest_eig] = min(abs(evals_full - lambda_star(idx)));
            cluster(nearest_eig) = true;
        end
        basis = Vfull(:, cluster);
        basis = orth(basis);
        degeneracy_multiplicity(j) = size(basis, 2);
        H = zeros(N*K, 1);
        layer_raw = zeros(K, 1);
        m_basis = zeros(K, size(basis, 2));
        s_basis = zeros(K, size(basis, 2));
        sync_sum = 0;
        node_raw = zeros(N, 1);
        for b = 1:size(basis, 2)
            M = reshape(basis(:, b), N, K)'; % 行对应层，列对应节点
            H = H + abs(basis(:, b)).^2;
            layer_raw = layer_raw + sum(abs(M).^2, 2);
            m_basis(:, b) = mean(M, 2);
            s_basis(:, b) = sqrt(sum(abs(M).^2, 2));
            node_raw = node_raw + sum(abs(M).^2, 1)';
            M_sync = repmat(mean(M, 1), K, 1);
            sync_sum = sync_sum + norm(M_sync, 'fro')^2;
        end
        selected_phi{j} = basis;
        selected_mode_components{j} = [];
        selected_eigenspace_energy{j} = reshape(H, N, K)'; % 基底无关的简并子空间能量 P_j
        layer_energy(:, j) = layer_raw / sum(layer_raw);
        layer_average_profile{j} = m_basis; % 简并时各列随 eig 返回的基底而变，不单独作物理解读
        layer_norm_profile{j} = s_basis;
        synchronous_fraction(j) = sync_sum / size(basis, 2);
        effective_nodes(j) = sum(node_raw)^2 / sum(node_raw.^2);
    else
        basis = phi_star{idx}(:);
        if numel(basis) ~= N*K
            error('beta/alpha=%.3f 处 phi_star 长度为 %d，应为 N*K=%d。', ...
                RATIO_LIST(idx), numel(basis), N*K);
        end
        basis = basis / norm(basis);
        [~, imax] = max(abs(basis));
        if real(basis(imax)) < 0, basis = -basis; end
        % 状态顺序是每层连续 N 个节点；reshape 成 N×K 后转置为 K×N。
        M = reshape(basis, N, K)';
        selected_phi{j} = basis;
        selected_mode_components{j} = M;
        selected_eigenspace_energy{j} = abs(M).^2;
        layer_raw = sum(abs(M).^2, 2);
        layer_energy(:, j) = layer_raw / sum(layer_raw);
        layer_average_profile{j} = mean(M, 2);
        layer_norm_profile{j} = sqrt(layer_raw);
        M_sync = repmat(mean(M, 1), K, 1);
        synchronous_fraction(j) = norm(M_sync, 'fro')^2 / norm(M, 'fro')^2;
        node_raw = sum(abs(M).^2, 1)';
        effective_nodes(j) = sum(node_raw)^2 / sum(node_raw.^2);
    end
    transverse_fraction(j) = max(0, 1 - synchronous_fraction(j));
    node_participation_fraction(j) = effective_nodes(j) / N;

    fprintf(['beta/alpha=%.1f，主导索引 q_star=%d，主导特征值 lambda_star=%.12g，' ...
        '谱参考阈值 sigma_c_spec=%.12g，左间隙=%.6g，右间隙=%.6g，' ...
        '是否简并=%d（重数=%d）\n'], ...
        selected_ratio_actual(j), selected_q(j), selected_lambda(j), ...
        selected_sigma_spec(j), gap_left(j), gap_right(j), ...
        degeneracy_flags(j), degeneracy_multiplicity(j));
end

%% 图 1：主导特征值与升序谱索引
switch_img = fullfile(fig_dir, 'dominant_mode_switch_alpha010.png');
f1 = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 1100 780]);
tl = tiledlayout(f1, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(tl);
plot(ax1, RATIO_LIST, lambda_star, '-o', 'LineWidth', 1.5, ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.0 0.447 0.741]);
yline(ax1, lambda_opt, '--', sprintf('\\lambda_{opt}=%.4g', lambda_opt), ...
    'LineWidth', 1.4, 'LabelHorizontalAlignment', 'left');
hold(ax1, 'on');
for j = 1:n_selected
    plot(ax1, selected_ratio_actual(j), selected_lambda(j), 'kp', ...
        'MarkerFaceColor', [1.0 0.75 0.0], 'MarkerSize', 10, 'HandleVisibility', 'off');
end
xlim(ax1, [min(RATIO_LIST) max(RATIO_LIST)]);
ylabel(ax1, '\lambda_*', 'Interpreter', 'tex');
title(ax1, '最危险超拉普拉斯模态的特征值');
grid(ax1, 'on'); box(ax1, 'on');
ax2 = nexttile(tl);
plot(ax2, RATIO_LIST, q_star, '-o', 'LineWidth', 1.3, ...
    'MarkerSize', 3.5, 'Color', [0.850 0.325 0.098]);
hold(ax2, 'on');
for j = 1:n_selected
    xline(ax2, selected_ratio_actual(j), ':', 'Color', [0.55 0.55 0.55], ...
        'HandleVisibility', 'off');
    plot(ax2, selected_ratio_actual(j), selected_q(j), 'kp', ...
        'MarkerFaceColor', [1.0 0.75 0.0], 'MarkerSize', 10, 'HandleVisibility', 'off');
    text(ax2, selected_ratio_actual(j), selected_q(j), sprintf('  %.1f', selected_ratio_actual(j)), ...
        'FontSize', 8, 'VerticalAlignment', 'bottom');
end
xlim(ax2, [min(RATIO_LIST) max(RATIO_LIST)]);
xlabel(ax2, '\beta/\alpha', 'Interpreter', 'tex');
ylabel(ax2, 'q_*（完整升序谱中的索引）');
title(ax2, '最危险离散模态的索引');
grid(ax2, 'on'); box(ax2, 'on');
exportgraphics(f1, switch_img, 'Resolution', 300);

%% 图 2：主导特征向量热图；若有简并，则统一显示子空间能量
heatmap_img = fullfile(fig_dir, 'dominant_mode_heatmaps_alpha010.png');
use_energy_heatmap = any(degeneracy_flags);
if use_energy_heatmap
    selected_heatmap = selected_eigenspace_energy;
else
    selected_heatmap = selected_mode_components;
end
heat_vals = cellfun(@(x) x(:), selected_heatmap, 'UniformOutput', false);
if use_energy_heatmap
    heat_cmax = max(cellfun(@max, heat_vals));
    if heat_cmax <= 0, heat_cmax = 1; end
else
    heat_cmax = max(cellfun(@(x) max(abs(x)), selected_heatmap));
    if heat_cmax <= 0, heat_cmax = 1; end
end
f2 = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [80 100 1500 500]);
tl2 = tiledlayout(f2, 1, n_selected, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:n_selected
    ax = nexttile(tl2);
    imagesc(ax, 1:N, 1:K, selected_heatmap{j});
    set(ax, 'YDir', 'normal', 'YTick', 1:K, 'FontSize', 9);
    xlabel(ax, '节点索引 i');
    if j == 1, ylabel(ax, '层 k'); end
    if use_energy_heatmap
        title(ax, sprintf('r=%.1f, \\lambda_*=%.4g, q_*=%d, \\sigma_c=%.4g', ...
            selected_ratio_actual(j), selected_lambda(j), selected_q(j), selected_sigma_spec(j)), ...
            'Interpreter', 'tex', 'FontSize', 9);
        clim(ax, [0 heat_cmax]);
    else
        title(ax, sprintf('r=%.1f, \\lambda_*=%.4g, q_*=%d, \\sigma_c=%.4g', ...
            selected_ratio_actual(j), selected_lambda(j), selected_q(j), selected_sigma_spec(j)), ...
            'Interpreter', 'tex', 'FontSize', 9);
        clim(ax, [-heat_cmax heat_cmax]);
    end
    box(ax, 'on');
end
if use_energy_heatmap
    colormap(f2, parula(256));
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = '主导特征子空间能量';
else
    colormap(f2, blueWhiteRed(256));
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = '主导特征向量分量';
end
exportgraphics(f2, heatmap_img, 'Resolution', 300);

%% 图 3：疑似切换区间内各层模态能量
energy_img = fullfile(fig_dir, 'dominant_mode_layer_energy_alpha010.png');
f3 = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 850 500]);
bar(1:K, layer_energy, 'grouped');
xlabel('层 k'); ylabel('模态能量占比');
title('beta/alpha=4.3--4.6 区间内主导模态的层间能量分布');
legend(arrayfun(@(r) sprintf('\\beta/\\alpha=%.1f', r), ...
    selected_ratio_actual, 'UniformOutput', false), 'Location', 'best');
grid on; box on;
exportgraphics(f3, energy_img, 'Resolution', 300);

%% 比较相邻参数点的跳变，重点检查 4.0--4.8 区间
local_edge_mask = RATIO_LIST(1:end-1) >= 4.0 & RATIO_LIST(2:end) <= 4.8 & ...
    isfinite(lambda_star(1:end-1)) & isfinite(lambda_star(2:end)) & ...
    isfinite(q_star(1:end-1)) & isfinite(q_star(2:end));
lambda_diffs = abs(diff(lambda_star));
q_diffs = abs(diff(q_star));
if any(local_edge_mask)
    local_edges = find(local_edge_mask);
    [~, jmax] = max(lambda_diffs(local_edge_mask));
    lambda_jump_edge = local_edges(jmax);
    [~, jmax] = max(q_diffs(local_edge_mask));
    q_jump_edge = local_edges(jmax);
else
    lambda_jump_edge = NaN;
    q_jump_edge = NaN;
end
if isfinite(lambda_jump_edge)
    lambda_jump_pair = RATIO_LIST([lambda_jump_edge, lambda_jump_edge + 1]);
    q_jump_pair = RATIO_LIST([q_jump_edge, q_jump_edge + 1]);
    fprintf('局部最大 |delta lambda_star|：%.6g -> %.6g（变化 %.6g）。\n', ...
        lambda_jump_pair(1), lambda_jump_pair(2), lambda_diffs(lambda_jump_edge));
    fprintf('局部最大 |delta q_star|：%.6g -> %.6g（变化 %.6g）。\n', ...
        q_jump_pair(1), q_jump_pair(2), q_diffs(q_jump_edge));
else
    lambda_jump_edge = NaN;
    q_jump_edge = NaN;
    lambda_jump_pair = [NaN NaN];
    q_jump_pair = [NaN NaN];
    fprintf('4.0--4.8 区间内没有可比较的有限相邻数据点。\n');
end

sync_labels = cell(1, n_selected);
localized_few_nodes = node_participation_fraction <= 0.1;
for j = 1:n_selected
    if synchronous_fraction(j) >= 0.9
        sync_labels{j} = '层间以同步成分为主';
    elseif synchronous_fraction(j) <= 0.1
        sync_labels{j} = '层间以横向成分为主';
    else
        sync_labels{j} = '层间结构混合';
    end
end
energy_tv_43_46 = 0.5 * sum(abs(layer_energy(:, 1) - layer_energy(:, end)));
fprintf('\n模态结构摘要（同步/横向标签使用 0.9/0.1 的判别界限）：\n');
fprintf('“少数节点局域”按有效参与节点数不超过 N 的 10%% 标记，仅作为便于阅读的判据。\n');
for j = 1:n_selected
    fprintf(['r=%.1f：%s；同步成分=%.4f，横向成分=%.4f，' ...
        '有效参与节点数=%.2f/%d（%.1f%%），少数节点局域=%d，简并=%d。\n'], ...
        selected_ratio_actual(j), sync_labels{j}, synchronous_fraction(j), ...
        transverse_fraction(j), effective_nodes(j), N, ...
        100*node_participation_fraction(j), localized_few_nodes(j), degeneracy_flags(j));
    fprintf('  各层能量占比 [L1..LK]：%s\n', mat2str(layer_energy(:, j)', 5));
    fprintf('  各层平均分量 m_k：%s\n', mat2str(layer_average_profile{j}, 5));
    fprintf('  各层模态范数 s_k：%s\n', mat2str(layer_norm_profile{j}', 5));
end
fprintf('r=4.3 与 r=4.6 的层能量总变差距离：%.6g。\n', ...
    energy_tv_43_46);
overlap_singular_values = svd(selected_phi{1}' * selected_phi{end});
mode_subspace_overlap = mean(overlap_singular_values.^2);
fprintf('r=4.3 与 r=4.6 的主导模态子空间重叠度（平方余弦均值）：%.6g。\n', ...
    mode_subspace_overlap);
fprintf('简并时 m_k 各列依赖 eig 返回的基底；此时以子空间能量热图和层能量占比为准。\n');

eigenvalue_gaps = table(selected_ratio_actual(:), selected_q(:), selected_lambda(:), ...
    selected_sigma_spec(:), selected_sigma_sim(:), gap_left(:), gap_right(:), ...
    degeneracy_tol(:), degeneracy_flags(:), degeneracy_multiplicity(:), ...
    'VariableNames', {'ratio','q_star','lambda_star','sigma_c_spec','sigma_c_sim', ...
    'gap_left','gap_right','degeneracy_tolerance','degeneracy_flag','multiplicity'});

base_out = fullfile(res_dir, 'dominant_mode_switch_alpha010.mat');
save(base_out, 'RATIO_LIST', 'lambda_star', 'q_star', 'sigma_c_spec', 'sigma_c_sim', ...
    'lambda_opt', 'sigma_min', 'selected_ratios', 'selected_ratio_actual', ...
    'selected_indices', 'selected_phi', 'selected_heatmap', 'layer_energy', ...
    'layer_average_profile', 'layer_norm_profile', 'eigenvalue_gaps', ...
    'gap_left', 'gap_right', 'degeneracy_tol', 'degeneracy_flags', ...
    'degeneracy_multiplicity', 'synchronous_fraction', 'transverse_fraction', ...
    'effective_nodes', 'node_participation_fraction', 'energy_tv_43_46', ...
    'localized_few_nodes', 'mode_subspace_overlap', 'lambda_jump_edge', ...
    'q_jump_edge', 'lambda_jump_pair', 'q_jump_pair', 'topology_file', 'input_file', ...
    'use_energy_heatmap', '-v7');

fprintf('\n输出文件：\n%s\n%s\n%s\n%s\n', ...
    switch_img, heatmap_img, energy_img, base_out);

function value = getField(S, candidates)
for i = 1:numel(candidates)
    if isfield(S, candidates{i})
        value = S.(candidates{i});
        fprintf('读取 MAT 字段：%s\n', candidates{i});
        return;
    end
end
error('MAT 文件缺少必需字段，可接受的字段名为：%s', strjoin(candidates, ', '));
end

function cmap = blueWhiteRed(n)
% 构造以白色为中点的蓝-白-红发散色图。
if nargin < 1, n = 256; end
anchors = [0.10 0.25 0.75; 1 1 1; 0.75 0.10 0.10];
x = linspace(0, 1, n);
cmap = interp1([0 0.5 1], anchors, x, 'linear');
end
