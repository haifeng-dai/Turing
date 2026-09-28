%% diagnose_ba_cross_layer_alignment.m: 已有 BA 多层网络的跨层枢纽节点对齐诊断
% 只读取已有网络文件，不重新生成或修改网络，也不运行动力学仿真。
% 诊断结果仅保存为图像，不输出终端统计信息或汇总数据文件。
clear; close all;

%% 1. 参数配置与 BA 文件扫描
script_dir = fileparts(mfilename('fullpath'));
ba_dir = fullfile(script_dir, 'results', 'topology', 'BA');
fig_dir = fullfile(script_dir, 'fig');
if ~exist(ba_dir, 'dir')
    error('找不到 BA 网络目录：%s', ba_dir);
end
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

N_PERMUTATIONS = 100;
N_HUB_BASELINE = 10000;
HUB_FRACTIONS = [0.05 0.10 0.20];
RNG_SEED = 20260927; % 仅控制诊断中的随机置换，不影响网络文件
rng(RNG_SEED, 'twister');

all_mat_files = dir(fullfile(ba_dir, '*.mat'));
file_records = struct('path', {}, 'name', {}, 'N_name', {}, 'm_name', {});
for file_idx = 1:numel(all_mat_files)
    token = regexp(all_mat_files(file_idx).name, '^N(\d+)_m(\d+)\.mat$', 'tokens', 'once');
    if isempty(token)
        continue;
    end
    file_records(end+1).path = fullfile(ba_dir, all_mat_files(file_idx).name); %#ok<SAGROW>
    file_records(end).name = all_mat_files(file_idx).name;
    file_records(end).N_name = str2double(token{1});
    file_records(end).m_name = str2double(token{2});
end
if isempty(file_records)
    error('目录中没有符合 N*_m*.mat 命名形式的 BA 网络文件：%s', ba_dir);
end
[~, file_order] = sortrows([[file_records.N_name]' [file_records.m_name]']);
file_records = file_records(file_order);

%% 2. 逐文件统计节点度相关、枢纽节点重叠与置换基线
n_files = numel(file_records);
file_diagnostics = cell(n_files, 1);
summary_filename = strings(n_files, 1);
summary_matrix_type = strings(n_files, 1);
summary_N = zeros(n_files, 1);
summary_m = zeros(n_files, 1);
summary_K = zeros(n_files, 1);
summary_mean_pearson = NaN(n_files, 1);
summary_mean_spearman = NaN(n_files, 1);
summary_median_spearman = NaN(n_files, 1);
summary_min_spearman = NaN(n_files, 1);
summary_max_spearman = NaN(n_files, 1);
summary_jaccard5 = NaN(n_files, 1);
summary_jaccard10 = NaN(n_files, 1);
summary_jaccard20 = NaN(n_files, 1);
summary_mean_rho_index_degree = NaN(n_files, 1);
summary_median_rho_index_degree = NaN(n_files, 1);
summary_min_rho_index_degree = NaN(n_files, 1);
summary_max_rho_index_degree = NaN(n_files, 1);
summary_early_hub5 = NaN(n_files, 1);
summary_early_hub10 = NaN(n_files, 1);
summary_early_hub20 = NaN(n_files, 1);
summary_perm_mean_spearman = NaN(n_files, 1);
summary_perm_spearman_lo = NaN(n_files, 1);
summary_perm_spearman_hi = NaN(n_files, 1);
summary_perm_spearman_percentile = NaN(n_files, 1);
summary_perm_p_spearman = NaN(n_files, 1);
summary_perm_mean_jaccard10 = NaN(n_files, 1);
summary_perm_jaccard10_lo = NaN(n_files, 1);
summary_perm_jaccard10_hi = NaN(n_files, 1);
summary_perm_jaccard10_percentile = NaN(n_files, 1);
summary_perm_p_jaccard10 = NaN(n_files, 1);
summary_mc_jaccard10_mean = NaN(n_files, 1);
summary_mc_jaccard10_lo = NaN(n_files, 1);
summary_mc_jaccard10_hi = NaN(n_files, 1);
summary_mc_jaccard10_95 = NaN(n_files, 1);
summary_mc_jaccard5_mean = NaN(n_files, 1);
summary_mc_jaccard5_lo = NaN(n_files, 1);
summary_mc_jaccard5_hi = NaN(n_files, 1);
summary_mc_jaccard5_95 = NaN(n_files, 1);
summary_mc_jaccard20_mean = NaN(n_files, 1);
summary_mc_jaccard20_lo = NaN(n_files, 1);
summary_mc_jaccard20_hi = NaN(n_files, 1);
summary_mc_jaccard20_95 = NaN(n_files, 1);
summary_degree_alignment = false(n_files, 1);
summary_hub_alignment = false(n_files, 1);
summary_age_order_effect = false(n_files, 1);
summary_conclusion = strings(n_files, 1);

for file_idx = 1:n_files
    file_path = file_records(file_idx).path;
    saved_data = load(file_path, 'nets');
    if ~isfield(saved_data, 'nets') || ~iscell(saved_data.nets) || isempty(saved_data.nets)
        continue;
    end

    nets = saved_data.nets;
    K = numel(nets);
    first_size = size(nets{1});
    if numel(first_size) ~= 2 || first_size(1) ~= first_size(2)
        continue;
    end
    N = first_size(1);
    degree_matrix = zeros(N, K);
    matrix_types = strings(1, K);
    valid_file = true;
    for layer_idx = 1:K
        M = nets{layer_idx};
        if ~isequal(size(M), [N N])
            valid_file = false;
            break;
        end
        [degree_matrix(:, layer_idx), matrix_types(layer_idx)] = ...
            matrixToDegree(M, file_records(file_idx).name, layer_idx);
    end
    if ~valid_file, continue; end
    if N ~= file_records(file_idx).N_name
        error('文件名中的 N=%d 与矩阵实际尺寸 N=%d 不一致：%s', ...
            file_records(file_idx).N_name, N, file_records(file_idx).name);
    end

    pearson_matrix = eye(K);
    spearman_matrix = eye(K);
    rank_matrix = zeros(N, K);
    index_vector = (1:N)';
    rho_index_degree = NaN(K, 1);
    for layer_idx = 1:K
        rank_matrix(:, layer_idx) = averageRank(degree_matrix(:, layer_idx));
        rho_index_degree(layer_idx) = correlationValue( ...
            averageRank(index_vector), rank_matrix(:, layer_idx));
    end
    for k = 1:K
        for ell = k+1:K
            pearson_matrix(k, ell) = correlationValue(degree_matrix(:, k), degree_matrix(:, ell));
            pearson_matrix(ell, k) = pearson_matrix(k, ell);
            spearman_matrix(k, ell) = correlationValue(rank_matrix(:, k), rank_matrix(:, ell));
            spearman_matrix(ell, k) = spearman_matrix(k, ell);
        end
    end
    upper_mask = triu(true(K), 1);
    pearson_values = pearson_matrix(upper_mask);
    spearman_values = spearman_matrix(upper_mask);
    mean_pearson = meanFinite(pearson_values);
    mean_spearman = meanFinite(spearman_values);
    median_spearman = medianFinite(spearman_values);

    hub_masks = cell(1, numel(HUB_FRACTIONS));
    hub_jaccard = cell(1, numel(HUB_FRACTIONS));
    early_hub_fraction = NaN(K, numel(HUB_FRACTIONS));
    for frac_idx = 1:numel(HUB_FRACTIONS)
        hub_count = max(1, round(HUB_FRACTIONS(frac_idx) * N));
        masks = false(N, K);
        for layer_idx = 1:K
            [~, order] = sort(degree_matrix(:, layer_idx), 'descend');
            masks(order(1:hub_count), layer_idx) = true;
            early_count = max(1, round(HUB_FRACTIONS(frac_idx) * N));
            early_hub_fraction(layer_idx, frac_idx) = ...
                sum(masks(1:early_count, layer_idx)) / hub_count;
        end
        hub_masks{frac_idx} = masks;
        hub_jaccard{frac_idx} = jaccardMatrix(masks);
    end
    jaccard5_values = hub_jaccard{1}(upper_mask);
    jaccard10_values = hub_jaccard{2}(upper_mask);
    jaccard20_values = hub_jaccard{3}(upper_mask);

    % 随机重标号零假设：每层独立置换节点标签，保留各层度数序列。
    null_mean_spearman = NaN(N_PERMUTATIONS, 1);
    null_mean_jaccard10 = NaN(N_PERMUTATIONS, 1);
    base_hub_masks10 = hub_masks{2};
    for perm_idx = 1:N_PERMUTATIONS
        permuted_ranks = zeros(N, K);
        permuted_hubs10 = false(N, K);
        for layer_idx = 1:K
            permutation = randperm(N);
            permuted_ranks(:, layer_idx) = rank_matrix(permutation, layer_idx);
            permuted_hubs10(:, layer_idx) = base_hub_masks10(permutation, layer_idx);
        end
        perm_spearman = NaN(nnz(upper_mask), 1);
        perm_hubs = NaN(nnz(upper_mask), 1);
        pair_idx = 0;
        for k = 1:K
            for ell = k+1:K
                pair_idx = pair_idx + 1;
                perm_spearman(pair_idx) = correlationValue( ...
                    permuted_ranks(:, k), permuted_ranks(:, ell));
                intersection_count = sum(permuted_hubs10(:, k) & permuted_hubs10(:, ell));
                union_count = sum(permuted_hubs10(:, k) | permuted_hubs10(:, ell));
                perm_hubs(pair_idx) = intersection_count / union_count;
            end
        end
        null_mean_spearman(perm_idx) = meanFinite(perm_spearman);
        null_mean_jaccard10(perm_idx) = meanFinite(perm_hubs);
    end

    % 对 top 5%、top 10%、top 20% 分别生成独立等大小集合的蒙特卡洛基线。
    mc_jaccard_baseline = cell(1, numel(HUB_FRACTIONS));
    mc_jaccard_interval = NaN(numel(HUB_FRACTIONS), 2);
    mc_jaccard_mean = NaN(1, numel(HUB_FRACTIONS));
    mc_jaccard_95 = NaN(1, numel(HUB_FRACTIONS));
    for frac_idx = 1:numel(HUB_FRACTIONS)
        hub_count = max(1, round(HUB_FRACTIONS(frac_idx) * N));
        mc_jaccard_baseline{frac_idx} = ...
            randomSetJaccardBaseline(N, hub_count, N_HUB_BASELINE);
        mc_jaccard_mean(frac_idx) = mean(mc_jaccard_baseline{frac_idx});
        mc_jaccard_interval(frac_idx, :) = ...
            quantileLinear(mc_jaccard_baseline{frac_idx}, [0.025 0.975]);
        mc_jaccard_95(frac_idx) = quantileLinear(mc_jaccard_baseline{frac_idx}, 0.95);
    end

    perm_spearman_mean = meanFinite(null_mean_spearman);
    perm_spearman_interval = quantileLinear(null_mean_spearman, [0.025 0.975]);
    perm_spearman_percentile = empiricalPercentile(null_mean_spearman, mean_spearman);
    perm_spearman_p = empiricalUpperP(null_mean_spearman, mean_spearman);
    perm_jaccard_mean = meanFinite(null_mean_jaccard10);
    perm_jaccard_interval = quantileLinear(null_mean_jaccard10, [0.025 0.975]);
    perm_jaccard_percentile = empiricalPercentile(null_mean_jaccard10, mean(jaccard10_values));
    perm_jaccard_p = empiricalUpperP(null_mean_jaccard10, mean(jaccard10_values));
    % 以单侧 5% 经验检验定义“高于随机重标号基线”。
    degree_alignment = isfinite(mean_spearman) && mean_spearman > 0 && perm_spearman_p <= 0.05;
    hub_alignment = isfinite(mean(jaccard10_values)) && perm_jaccard_p <= 0.05;
    mean_rho_index = meanFinite(rho_index_degree);
    age_order_effect = mean_rho_index < 0 && quantileLinear(rho_index_degree, 0.95) < 0;
    if degree_alignment && hub_alignment && age_order_effect
        conclusion = "支持跨层 hub 对齐，且节点编号/共同加入顺序很可能有贡献";
    elseif degree_alignment && hub_alignment
        conclusion = "支持跨层 hub 对齐，但节点编号年龄效应未在所有层稳定显示";
    elseif degree_alignment || hub_alignment
        conclusion = "存在部分跨层对齐证据，需结合各项统计量判断";
    else
        conclusion = "当前统计未显示显著的跨层 hub 对齐证据";
    end

    type_list = unique(matrix_types, 'stable');
    type_summary = strjoin(type_list, '、');
    summary_filename(file_idx) = string(file_records(file_idx).name);
    summary_matrix_type(file_idx) = type_summary;
    summary_N(file_idx) = N;
    summary_m(file_idx) = file_records(file_idx).m_name;
    summary_K(file_idx) = K;
    summary_mean_pearson(file_idx) = mean_pearson;
    summary_mean_spearman(file_idx) = mean_spearman;
    summary_median_spearman(file_idx) = median_spearman;
    summary_min_spearman(file_idx) = minFinite(spearman_values);
    summary_max_spearman(file_idx) = maxFinite(spearman_values);
    summary_jaccard5(file_idx) = meanFinite(jaccard5_values);
    summary_jaccard10(file_idx) = meanFinite(jaccard10_values);
    summary_jaccard20(file_idx) = meanFinite(jaccard20_values);
    summary_mean_rho_index_degree(file_idx) = mean_rho_index;
    summary_median_rho_index_degree(file_idx) = medianFinite(rho_index_degree);
    summary_min_rho_index_degree(file_idx) = minFinite(rho_index_degree);
    summary_max_rho_index_degree(file_idx) = maxFinite(rho_index_degree);
    summary_early_hub5(file_idx) = meanFinite(early_hub_fraction(:, 1));
    summary_early_hub10(file_idx) = meanFinite(early_hub_fraction(:, 2));
    summary_early_hub20(file_idx) = meanFinite(early_hub_fraction(:, 3));
    summary_perm_mean_spearman(file_idx) = perm_spearman_mean;
    summary_perm_spearman_lo(file_idx) = perm_spearman_interval(1);
    summary_perm_spearman_hi(file_idx) = perm_spearman_interval(2);
    summary_perm_spearman_percentile(file_idx) = perm_spearman_percentile;
    summary_perm_p_spearman(file_idx) = perm_spearman_p;
    summary_perm_mean_jaccard10(file_idx) = perm_jaccard_mean;
    summary_perm_jaccard10_lo(file_idx) = perm_jaccard_interval(1);
    summary_perm_jaccard10_hi(file_idx) = perm_jaccard_interval(2);
    summary_perm_jaccard10_percentile(file_idx) = perm_jaccard_percentile;
    summary_perm_p_jaccard10(file_idx) = perm_jaccard_p;
    summary_mc_jaccard5_mean(file_idx) = mc_jaccard_mean(1);
    summary_mc_jaccard5_lo(file_idx) = mc_jaccard_interval(1, 1);
    summary_mc_jaccard5_hi(file_idx) = mc_jaccard_interval(1, 2);
    summary_mc_jaccard5_95(file_idx) = mc_jaccard_95(1);
    summary_mc_jaccard10_mean(file_idx) = mc_jaccard_mean(2);
    summary_mc_jaccard10_lo(file_idx) = mc_jaccard_interval(2, 1);
    summary_mc_jaccard10_hi(file_idx) = mc_jaccard_interval(2, 2);
    summary_mc_jaccard10_95(file_idx) = mc_jaccard_95(2);
    summary_mc_jaccard20_mean(file_idx) = mc_jaccard_mean(3);
    summary_mc_jaccard20_lo(file_idx) = mc_jaccard_interval(3, 1);
    summary_mc_jaccard20_hi(file_idx) = mc_jaccard_interval(3, 2);
    summary_mc_jaccard20_95(file_idx) = mc_jaccard_95(3);
    summary_degree_alignment(file_idx) = degree_alignment;
    summary_hub_alignment(file_idx) = hub_alignment;
    summary_age_order_effect(file_idx) = age_order_effect;
    summary_conclusion(file_idx) = conclusion;

    file_diagnostics{file_idx} = struct( ...
        'filename', file_records(file_idx).name, 'N', N, 'm', file_records(file_idx).m_name, ...
        'K', K, 'matrix_types', matrix_types, 'degree_matrix', degree_matrix, ...
        'pearson_matrix', pearson_matrix, 'spearman_matrix', spearman_matrix, ...
        'hub_jaccard_5', hub_jaccard{1}, 'hub_jaccard_10', hub_jaccard{2}, ...
        'hub_jaccard_20', hub_jaccard{3}, 'rho_index_degree', rho_index_degree, ...
        'early_hub_fraction', early_hub_fraction, ...
        'null_mean_spearman', null_mean_spearman, 'null_mean_jaccard10', null_mean_jaccard10, ...
        'mc_jaccard5', mc_jaccard_baseline{1}, ...
        'mc_jaccard10', mc_jaccard_baseline{2}, ...
        'mc_jaccard20', mc_jaccard_baseline{3}, 'degree_alignment', degree_alignment, ...
        'hub_alignment', hub_alignment, 'age_order_effect', age_order_effect, ...
        'conclusion', conclusion);

end

%% 3. 检查分析结果
valid_rows = summary_N > 0;
if ~any(valid_rows)
    error('没有找到可分析的 BA 网络 MAT 文件。');
end

%% 4. 选择代表网络并绘图
% 优先使用 N=200、m=20；若文件不存在，则选择 N 最接近 200 的网络。
valid_record_indices = find(valid_rows);
rep_idx = find([file_records.N_name] == 200 & [file_records.m_name] == 20, 1);
if isempty(rep_idx) || isempty(file_diagnostics{rep_idx})
    valid_N = [file_records(valid_record_indices).N_name];
    [~, near_idx] = min(abs(valid_N - 200));
    rep_idx = valid_record_indices(near_idx);
end
rep = file_diagnostics{rep_idx};
rep_tag = sprintf('N%d_m%d', rep.N, rep.m);
rep_layer_count = rep.K;
rep_layers = unique(round(linspace(1, rep_layer_count, min(3, rep_layer_count))));

% 1. 跨层 Spearman 节点度相关矩阵
fig1 = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 720 620]);
imagesc(rep.spearman_matrix, [-1 1]); axis image; colorbar;
colormap(fig1, blueWhiteRed(256));
set(gca, 'XTick', 1:rep.K, 'YTick', 1:rep.K);
xlabel('层'); ylabel('层');
title(sprintf('%s：跨层 degree Spearman 相关', rep.filename), 'Interpreter', 'none');
fig1_path = fullfile(fig_dir, ['ba_cross_layer_spearman_', rep_tag, '.png']);
exportgraphics(fig1, fig1_path, 'Resolution', 300);

% 2. top 10% 枢纽节点的跨层 Jaccard 矩阵
fig2 = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 720 620]);
imagesc(rep.hub_jaccard_10, [0 1]); axis image; colorbar;
colormap(fig2, parula(256));
set(gca, 'XTick', 1:rep.K, 'YTick', 1:rep.K);
xlabel('层'); ylabel('层');
title(sprintf('%s：top 10%% hub 的 Jaccard 重叠', rep.filename), 'Interpreter', 'none');
fig2_path = fullfile(fig_dir, ['ba_cross_layer_hub_jaccard10_', rep_tag, '.png']);
exportgraphics(fig2, fig2_path, 'Resolution', 300);

% 3. 代表层的节点编号与度数关系
fig3 = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 1000 560]);
hold on;
line_colors = lines(numel(rep_layers));
legend_text = cell(1, numel(rep_layers));
for layer_plot_idx = 1:numel(rep_layers)
    layer_idx = rep_layers(layer_plot_idx);
    rho_value = rep.rho_index_degree(layer_idx);
    plot(1:rep.N, rep.degree_matrix(:, layer_idx), '-', ...
        'Color', line_colors(layer_plot_idx, :), 'LineWidth', 1.1);
    legend_text{layer_plot_idx} = sprintf('第 %d 层，rho=%.3f', layer_idx, rho_value);
end
xlabel('节点编号 i'); ylabel('degree');
title(sprintf('%s：节点编号与度数（代表层）', rep.filename), 'Interpreter', 'none');
legend(legend_text, 'Location', 'best'); grid on; box on;
fig3_path = fullfile(fig_dir, ['ba_degree_by_node_', rep_tag, '.png']);
exportgraphics(fig3, fig3_path, 'Resolution', 300);

% 4. 原始统计量与随机重标号零假设分布
rep_summary_row = find(summary_filename == string(rep.filename), 1);
rep_diag = rep;
fig4 = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 1100 500]);
tl = tiledlayout(fig4, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
histogram(ax, rep_diag.null_mean_spearman, 'FaceColor', [0.35 0.65 0.85]);
hold(ax, 'on');
xline(ax, summary_mean_spearman(rep_summary_row), 'r-', '原始值', 'LineWidth', 2);
xline(ax, summary_perm_spearman_hi(rep_summary_row), 'k--', '置换 95% 上界', 'LineWidth', 1.4);
xlabel(ax, '跨层平均 Spearman 相关'); ylabel(ax, '次数');
title(ax, '随机重标号：degree 相关零假设'); grid(ax, 'on'); box(ax, 'on');
ax = nexttile(tl);
histogram(ax, rep_diag.null_mean_jaccard10, 'FaceColor', [0.50 0.70 0.40]);
hold(ax, 'on');
xline(ax, summary_jaccard10(rep_summary_row), 'r-', '原始值', 'LineWidth', 2);
xline(ax, summary_perm_jaccard10_hi(rep_summary_row), 'k--', '置换 95% 上界', 'LineWidth', 1.4);
xline(ax, summary_mc_jaccard10_mean(rep_summary_row), 'b:', '独立随机集合均值', 'LineWidth', 1.6);
xlabel(ax, '跨层平均 top 10% hub Jaccard'); ylabel(ax, '次数');
title(ax, '随机重标号：hub 重叠零假设'); grid(ax, 'on'); box(ax, 'on');
fig4_path = fullfile(fig_dir, ['ba_permutation_baseline_', rep_tag, '.png']);
exportgraphics(fig4, fig4_path, 'Resolution', 300);

%% 5. 关闭图窗并清理工作区
close([fig1 fig2 fig3 fig4]);
clear;

%% 6. 局部函数
function [degree, matrix_type] = matrixToDegree(M, filename, layer_idx)
% 根据对称性、行和及非对角元符号识别邻接矩阵或拉普拉斯矩阵。
if size(M, 1) ~= size(M, 2)
    error('%s 第 %d 层不是方阵。', filename, layer_idx);
end
scale = max(1, norm(M, 'fro'));
tol = 1e-10 * scale;
is_symmetric = norm(M - M', 'fro') <= tol;
row_sum_error = norm(full(sum(M, 2)), inf);
diag_values = full(diag(M));
offdiag = M - spdiags(diag(M), 0, size(M, 1), size(M, 2));
offdiag_values = nonzeros(offdiag);
if isempty(offdiag_values), max_offdiag = 0; else, max_offdiag = max(offdiag_values); end
if is_symmetric && row_sum_error <= tol && all(diag_values >= -tol) && max_offdiag <= tol
    degree = diag_values;
    matrix_type = "拉普拉斯矩阵";
    return;
end
if is_symmetric && max(abs(diag_values)) <= tol && ...
        (isempty(nonzeros(M)) || min(nonzeros(M)) >= -tol)
    degree = full(sum(M, 2));
    matrix_type = "邻接矩阵";
    return;
end
error(['无法自动识别 %s 第 %d 层是邻接矩阵还是拉普拉斯矩阵。' ...
    '请检查矩阵是否对称、拉普拉斯行和是否为零。'], filename, layer_idx);
end

function ranks = averageRank(values)
% 计算升序平均秩，避免依赖统计工具箱函数 tiedrank。
values = values(:);
[sorted_values, order] = sort(values);
ranks = zeros(size(values));
first = 1;
while first <= numel(values)
    last = first;
    while last < numel(values) && sorted_values(last + 1) == sorted_values(first)
        last = last + 1;
    end
    ranks(order(first:last)) = (first + last) / 2;
    first = last + 1;
end
end

function value = correlationValue(x, y)
% 计算皮尔逊相关系数；常数向量的相关系数记为 NaN。
x = double(x(:));
y = double(y(:));
valid = isfinite(x) & isfinite(y);
x = x(valid); y = y(valid);
if numel(x) < 2
    value = NaN;
    return;
end
x = x - mean(x); y = y - mean(y);
denominator = sqrt(sum(x.^2) * sum(y.^2));
if denominator == 0
    value = NaN;
else
    value = sum(x .* y) / denominator;
end
end

function matrix = jaccardMatrix(masks)
% 对每两层的枢纽节点集合计算 Jaccard 系数。
K = size(masks, 2);
matrix = eye(K);
for k = 1:K
    for ell = k+1:K
        intersection_count = sum(masks(:, k) & masks(:, ell));
        union_count = sum(masks(:, k) | masks(:, ell));
        matrix(k, ell) = intersection_count / union_count;
        matrix(ell, k) = matrix(k, ell);
    end
end
end

function values = randomSetJaccardBaseline(N, set_size, n_samples)
% 从独立、等大小随机集合的超几何交集分布抽样 Jaccard 基线。
overlap_min = max(0, 2*set_size - N);
overlap_values = (overlap_min:set_size)';
log_choose = @(n, k) gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1);
log_probability = log_choose(set_size, overlap_values) + ...
    log_choose(N - set_size, set_size - overlap_values) - log_choose(N, set_size);
probability = exp(log_probability - max(log_probability));
probability = probability / sum(probability);
cumulative = cumsum(probability);
cumulative(end) = 1;
uniform_draws = rand(n_samples, 1);
lower_index = ones(n_samples, 1);
upper_index = numel(cumulative) * ones(n_samples, 1);
while any(lower_index < upper_index)
    active = lower_index < upper_index;
    mid = floor((lower_index(active) + upper_index(active)) / 2);
    draw = uniform_draws(active);
    go_left = draw <= cumulative(mid);
    active_indices = find(active);
    left_indices = active_indices(go_left);
    right_indices = active_indices(~go_left);
    upper_index(left_indices) = mid(go_left);
    lower_index(right_indices) = mid(~go_left) + 1;
end
overlap = overlap_values(lower_index);
values = overlap ./ (2*set_size - overlap);
end

function q = quantileLinear(values, probabilities)
% 通过线性插值计算分位数，避免依赖统计工具箱函数 prctile。
values = sort(values(isfinite(values)));
probabilities = probabilities(:)';
q = NaN(size(probabilities));
if isempty(values), return; end
for idx = 1:numel(probabilities)
    position = 1 + (numel(values) - 1) * probabilities(idx);
    lower = floor(position);
    upper = ceil(position);
    if lower == upper
        q(idx) = values(lower);
    else
        q(idx) = values(lower) + (position - lower) * (values(upper) - values(lower));
    end
end
end

function value = empiricalPercentile(null_values, observed)
% 计算观测值在置换分布中的经验百分位。
null_values = null_values(isfinite(null_values));
if isempty(null_values) || ~isfinite(observed)
    value = NaN;
else
    value = 100 * (sum(null_values < observed) + 0.5 * sum(null_values == observed)) / numel(null_values);
end
end

function value = empiricalUpperP(null_values, observed)
% 单侧经验 p 值，检验原始统计量是否高于随机重标号基线。
null_values = null_values(isfinite(null_values));
if isempty(null_values) || ~isfinite(observed)
    value = NaN;
else
    value = (1 + sum(null_values >= observed)) / (numel(null_values) + 1);
end
end

function value = meanFinite(values)
values = values(isfinite(values));
if isempty(values), value = NaN; else, value = mean(values); end
end

function value = medianFinite(values)
values = values(isfinite(values));
if isempty(values), value = NaN; else, value = median(values); end
end

function value = minFinite(values)
values = values(isfinite(values));
if isempty(values), value = NaN; else, value = min(values); end
end

function value = maxFinite(values)
values = values(isfinite(values));
if isempty(values), value = NaN; else, value = max(values); end
end

function cmap = blueWhiteRed(n)
% 构造以白色为中点的蓝-白-红发散色图。
anchors = [0.10 0.25 0.75; 1 1 1; 0.75 0.10 0.10];
cmap = interp1([0 0.5 1], anchors, linspace(0, 1, n), 'linear');
end
