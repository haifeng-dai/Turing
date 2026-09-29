%% analyze_dispersion_and_dominant_modes_alpha010.m: 色散关系与主导模态空间结构
% 只读取已有谱分析结果和网络实例，不生成网络，也不运行动力学仿真。
clear; clc; close all;

%% 1. 路径、参数与已有谱结果
script_dir = fileparts(mfilename('fullpath'));
res_dir = fullfile(script_dir, 'results');
fig_dir = fullfile(script_dir, 'fig');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

input_file = fullfile(res_dir, ...
    'compare_forward_threshold_spectral_er_N200_K5_p0.030_a0.100.mat');
if ~exist(input_file, 'file')
    error('找不到已有谱分析结果文件：%s', input_file);
end
S = load(input_file);

RATIO_LIST = getRequiredField(S, {'RATIO_LIST', 'ratio_list'});
lambda_all = getRequiredField(S, {'lambda_all', 'LAMBDA_ALL'});
lambda_star_all = getRequiredField(S, {'lambda_star', 'LAMBDA_STAR'});
q_star_all = getRequiredField(S, {'q_star', 'Q_STAR'});
phi_star_all = getRequiredField(S, {'phi_star', 'PHI_STAR'});
sigma_c_spec_all = getRequiredField(S, {'sigma_c_spec', 'SIGMA_C_SPEC'});

RATIO_LIST = RATIO_LIST(:)';
lambda_star_all = lambda_star_all(:)';
q_star_all = q_star_all(:)';
sigma_c_spec_all = sigma_c_spec_all(:)';
if isfield(S, 'N'), N = double(S.N); else, N = 200; end
if isfield(S, 'K'), K = double(S.K); else, K = 5; end
if isfield(S, 'ALPHA_FIXED')
    alpha = double(S.ALPHA_FIXED);
elseif isfield(S, 'alpha')
    alpha = double(S.alpha);
else
    alpha = 0.1;
end
if abs(alpha - 0.1) > 1e-12
    error('谱结果中的 alpha=%.12g，与本脚本要求的 alpha=0.1 不一致。', alpha);
end
if any([numel(lambda_star_all), numel(q_star_all), numel(sigma_c_spec_all), ...
        numel(lambda_all)] ~= numel(RATIO_LIST))
    error('谱分析 MAT 文件中的 ratio 与谱结果字段长度不一致。');
end
if N ~= 200 || K ~= 5
    warning('当前结果为 N=%d、K=%d；预期诊断配置为 N=200、K=5。', N, K);
end

selected_ratios = [4.4 4.5];
n_selected = numel(selected_ratios);
selected_indices = zeros(1, n_selected);
for j = 1:n_selected
    [ratio_error, selected_indices(j)] = min(abs(RATIO_LIST - selected_ratios(j)));
    if ratio_error > 1e-7
        error('已有 RATIO_LIST 中找不到 r=%.1f。', selected_ratios(j));
    end
end

lambda_star = lambda_star_all(selected_indices);
q_star = q_star_all(selected_indices);
sigma_c_spec = sigma_c_spec_all(selected_indices);
lambda_all_selected = cell(1, n_selected);
mu_all_selected = cell(1, n_selected);
mu_cont = cell(1, n_selected);
dominant_basis = cell(1, n_selected);
critical_modes = cell(1, n_selected);
critical_projector_diag = cell(1, n_selected);
IPR_full = NaN(1, n_selected);
% 对完整 NK 维单位归一化 supra-eigenvector，完全离域时 IPR=1/(NK)。
IPR_reference = 1/(N*K);
% 仅作数值一致性检查，不作为论文正式指标。
sync_projection_norm = NaN(1, n_selected);
q_danger = NaN(1, n_selected);
mu_max = NaN(1, n_selected);
multiplicity = NaN(1, n_selected);
gap_left = NaN(1, n_selected);
gap_right = NaN(1, n_selected);

for j = 1:n_selected
    idx = selected_indices(j);
    lambda_all_selected{j} = real(getSpectrumAt(lambda_all, idx));
    lambda_all_selected{j} = lambda_all_selected{j}(:);
    if numel(lambda_all_selected{j}) ~= N*K
        error('r=%.1f 的超拉普拉斯谱长度为 %d，应为 N*K=%d。', ...
            selected_ratios(j), numel(lambda_all_selected{j}), N*K);
    end
    if any(diff(lambda_all_selected{j}) < -1e-10)
        error('r=%.1f 的 lambda_all 不是升序排列。', selected_ratios(j));
    end
    if ~isfinite(sigma_c_spec(j)) || sigma_c_spec(j) <= 0
        error('r=%.1f 的 sigma_c_spec 无效。', selected_ratios(j));
    end
end

%% 2. 数值计算连续色散关系的最优谱位置
lambda_upper = (10/3) / alpha;
sigma_of_lambda = @(lambda) (110/3 + 4*alpha.*lambda) ./ ...
    ((alpha.*lambda) .* (10/3 - alpha.*lambda));
lambda_grid_opt = linspace(lambda_upper*1e-8, lambda_upper*(1-1e-8), 100001);
sigma_grid_opt = sigma_of_lambda(lambda_grid_opt);
[~, opt_grid_idx] = min(sigma_grid_opt);
opt_lower = lambda_grid_opt(max(1, opt_grid_idx-1));
opt_upper = lambda_grid_opt(min(numel(lambda_grid_opt), opt_grid_idx+1));
lambda_opt = fminbnd(sigma_of_lambda, opt_lower, opt_upper);
sigma_min = sigma_of_lambda(lambda_opt);

%% 3. 读取已有网络并计算离散与连续增长率
topology_type = getOptionalField(S, {'TOPO_TYPE', 'topology_type'}, 'ER');
topology_type = upper(char(topology_type));
topology_param = getOptionalField(S, ...
    {'P_VAL_FIXED', 'TOPO_PARAM', 'topology_param'}, 0.030);
if strcmp(topology_type, 'BA')
    topology_param = getOptionalField(S, ...
        {'M_VAL_FIXED', 'M_FIXED', 'TOPO_PARAM'}, topology_param);
end
topology_file = resolveTopologyFile(S, script_dir, topology_type, N, topology_param);
if ~exist(topology_file, 'file')
    error('找不到已有网络实例：%s', topology_file);
end
topology_data = load(topology_file, 'nets');
if ~isfield(topology_data, 'nets') || numel(topology_data.nets) < K
    error('已有网络文件层数不足，需要 %d 层：%s', K, topology_file);
end

L_intra = cell(K, 1);
for layer_idx = 1:K
    Lk = sparse(topology_data.nets{layer_idx});
    if ~isequal(size(Lk), [N N])
        error('网络第 %d 层尺寸不是 N×N：%s', layer_idx, topology_file);
    end
    if norm(full(sum(Lk, 2)), inf) > 1e-8 || ...
            norm(full(Lk-Lk'), 'fro') > 1e-8 || any(diag(Lk) < -1e-12)
        error('网络第 %d 层不是当前分析要求的对称拉普拉斯矩阵。', layer_idx);
    end
    L_intra{layer_idx} = Lk;
end
L_intra_all = blkdiag(L_intra{:});
adj_inter = ones(K) - eye(K);
L_inter = diag(sum(adj_inter, 2)) - adj_inter;
L_inter_big = kron(L_inter, speye(N));

J = [10/3, -5; 10, -4];
degeneracy_tolerance = zeros(1, n_selected);
% 层同步子空间投影仅用于内部数值一致性检查。
P_sync = kron(ones(K)/K, speye(N));
for j = 1:n_selected
    idx = selected_indices(j);
    spectrum = lambda_all_selected{j};
    if abs(spectrum(1)) > 1e-8
        warning('r=%.1f 的最小特征值并非数值零：lambda_1=%.6g。', ...
            selected_ratios(j), spectrum(1));
    end

    q_saved = q_star(j);
    if ~isfinite(q_saved) || abs(q_saved-round(q_saved)) > 1e-8 || ...
            q_saved < 1 || q_saved > numel(spectrum)
        error('r=%.1f 的 q_star 不是有效谱索引。', selected_ratios(j));
    end
    q_saved = round(q_saved);
    if abs(spectrum(q_saved)-lambda_star(j)) > 1e-6*max(1,abs(lambda_star(j)))
        warning('r=%.1f 的 lambda_all(q_star) 与 lambda_star 不一致。', selected_ratios(j));
    end
    if q_saved > 1, gap_left(j) = spectrum(q_saved)-spectrum(q_saved-1); end
    if q_saved < numel(spectrum), gap_right(j) = spectrum(q_saved+1)-spectrum(q_saved); end

    degeneracy_tolerance(j) = max(1e-8, 1e-6*max(1,abs(lambda_star(j))));
    degenerate_indices = find(abs(spectrum-lambda_star(j)) < degeneracy_tolerance(j));
    multiplicity(j) = numel(degenerate_indices);
    if isempty(degenerate_indices)
        error('r=%.1f 的 lambda_star 在 lambda_all 中没有对应特征值。', ...
            selected_ratios(j));
    end
    if ~ismember(q_saved, degenerate_indices)
        warning('r=%.1f 的 q_star 不在 lambda_star 对应的特征值簇中。', ...
            selected_ratios(j));
    end

    mu_all_selected{j} = growthRates(spectrum, alpha, sigma_c_spec(j), J);
    [mu_max(j), q_local] = max(mu_all_selected{j}(2:end));
    q_danger(j) = q_local + 1;
    lambda_danger = spectrum(q_danger(j));
    if ~ismember(q_danger(j), degenerate_indices)
        warning(['r=%.1f 的最危险非均匀模态 q=%d 未落在 lambda_star 的简并簇中；' ...
            '请检查谱结果和网络 realization。'], selected_ratios(j), q_danger(j));
    end
    if abs(lambda_danger-lambda_star(j)) > degeneracy_tolerance(j)
        warning('r=%.1f 的最危险特征值与已保存 lambda_star 偏差较大。', selected_ratios(j));
    end

    L_super = L_intra_all + RATIO_LIST(idx)*L_inter_big;
    if multiplicity(j) > 1
        % 简并时从已有网络重算特征子空间，避免解释任意的单个特征向量。
        [V_full, D_full] = eig(full(L_super));
        [lambda_recomputed, order] = sort(real(diag(D_full)), 'ascend');
        V_full = V_full(:, order);
        spectrum_error = max(abs(lambda_recomputed-spectrum));
        if spectrum_error > 1e-6*max(1,max(abs(spectrum)))
            error('r=%.1f 的已有网络谱与 MAT 文件不一致，无法可靠重建简并子空间。', ...
                selected_ratios(j));
        end
        cluster = abs(lambda_recomputed-lambda_star(j)) < degeneracy_tolerance(j);
        V_mode = orth(V_full(:, cluster));
        multiplicity(j) = size(V_mode, 2);
    else
        V_mode = getEigenvectorAt(phi_star_all, idx);
        if numel(V_mode) ~= N*K || norm(V_mode) == 0
            error('r=%.1f 的 phi_star 长度无效。', selected_ratios(j));
        end
        V_mode = V_mode(:)/norm(V_mode);
        residual = norm(L_super*V_mode-lambda_star(j)*V_mode) / max(norm(V_mode), eps);
        if residual > 1e-6*max(1,abs(lambda_star(j)))
            warning('r=%.1f 的 phi_star 与已有网络特征向量残差较大：%.4g。', ...
                selected_ratios(j), residual);
        end
    end
    dominant_basis{j} = V_mode;

    % 临界模态或临界特征子空间的空间表示。
    if multiplicity(j) > 1
        % 由于临界特征值简并，不存在唯一临界特征向量。
        % 显示临界特征子空间正交投影矩阵的对角线；该量对
        % 简并子空间内部的正交基变换不敏感。
        projector_diag = sum(abs(V_mode).^2, 2);
        critical_projector_diag{j} = reshape(projector_diag, [N, K])';
        critical_modes{j} = [];
        IPR_full(j) = NaN;
        sync_projection_norm(j) = ...
            norm(P_sync*V_mode, 'fro') / norm(V_mode, 'fro');
    else
        % 非简并时使用前面已经校验过的真实主导特征向量。
        phi = real(V_mode(:));
        phi = phi / norm(phi, 2);

        % 固定整体符号以便图形可重复；整体正负号不改变特征向量含义。
        [~, imax] = max(abs(phi));
        if phi(imax) < 0
            phi = -phi;
        end

        critical_modes{j} = reshape(phi, [N, K])';
        critical_projector_diag{j} = [];
        IPR_full(j) = sum(abs(phi).^4) / (sum(abs(phi).^2)^2);
        sync_projection_norm(j) = ...
            norm((speye(N*K)-P_sync)*phi) / norm(phi);
    end

    if abs(selected_ratios(j)-4.4) < 1e-8 && sync_projection_norm(j) > 1e-6
        warning('r=4.4 的临界子空间并非纯 transverse，请检查。');
    end
end

lambda_cont_max = max(cellfun(@max, lambda_all_selected));
lambda_plot_max = 1.03 * max(lambda_cont_max, lambda_opt);
lambda_cont = linspace(0, lambda_plot_max, 2000);
for j = 1:n_selected
    mu_cont{j} = growthRates(lambda_cont, alpha, sigma_c_spec(j), J);
end

%% 4. 绘制色散关系与临界模态表示
all_mu_values = [mu_all_selected{1}(:); mu_all_selected{2}(:); ...
    mu_cont{1}(:); mu_cont{2}(:); 0];
mu_min_plot = min(all_mu_values);
mu_max_plot = max(all_mu_values);
mu_padding = max(0.04*(mu_max_plot-mu_min_plot), 1e-3);
common_ylim = [mu_min_plot-mu_padding, mu_max_plot+mu_padding];

out_img = fullfile(fig_dir, 'com2_3.png');
out_mat = fullfile(res_dir, 'com2_3.mat');
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [80 80 1350 900]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for j = 1:n_selected
    ax = nexttile(layout, j);
    hold(ax, 'on');
    plot(ax, lambda_cont, mu_cont{j}, 'k-', 'LineWidth', 1.8, ...
        'DisplayName', 'Continuous dispersion');
    spectrum = lambda_all_selected{j};
    mu_values = mu_all_selected{j};
    scatter(ax, spectrum(2:end), mu_values(2:end), ...
        13, [0.20 0.48 0.78], 'filled', 'MarkerFaceAlpha', 0.45, ...
        'DisplayName', 'Discrete modes');
    scatter(ax, spectrum(1), mu_values(1), ...
        45, [0.65 0.65 0.65], 'o', 'filled', 'DisplayName', 'Uniform mode');
    scatter(ax, lambda_star(j), mu_max(j), 110, ...
        [0.85 0.20 0.16], 'p', 'filled', 'MarkerEdgeColor', 'k', ...
        'DisplayName', 'Dominant mode');
    yline(ax, 0, 'k:', 'HandleVisibility', 'off');
    xline(ax, lambda_star(j), '--', 'Color', [0.85 0.20 0.16], ...
        'HandleVisibility', 'off');
    xline(ax, lambda_opt, '-.', 'Color', [0.35 0.60 0.35], ...
        'DisplayName', '\lambda_{opt}');
    xlim(ax, [0 lambda_plot_max]);
    ylim(ax, common_ylim);
    xlabel(ax, 'Supra-Laplacian eigenvalue \lambda');
    ylabel(ax, 'Maximum growth rate \mu');
    title(ax, sprintf('(%c) r=%.1f, \\sigma_c=%.4f', ...
        char('a'+j-1), selected_ratios(j), sigma_c_spec(j)), 'Interpreter', 'tex');
    legend(ax, 'show', 'Location', 'best', 'FontSize', 8);
    grid(ax, 'on'); box(ax, 'on');
end

signed_map = blueWhiteRedMap(257);
signed_color_limit = 0.4;
for j = 1:n_selected
    ax = nexttile(layout, j+2);
    if multiplicity(j) > 1
        Mproj = critical_projector_diag{j};
        imagesc(ax, 1:N, 1:K, Mproj);
        set(ax, 'YDir', 'normal', 'YTick', 1:K);
        % 与右下图共用同一套蓝—白—红色图。
        colormap(ax, signed_map);

        % 与右下特征向量图统一使用固定的对称色标范围。
        clim(ax, [-signed_color_limit signed_color_limit]);

        xlabel(ax, 'Node index i');
        ylabel(ax, 'Layer k');
        title(ax, sprintf(['(%c) Critical eigenspace projector\n' ...
            'r=%.1f, \\lambda_*=%.4f, multiplicity=%d'], ...
            char('a'+j+1), selected_ratios(j), ...
            lambda_star(j), multiplicity(j)), 'Interpreter', 'tex');
        cb = colorbar(ax, 'Location', 'eastoutside');
        cb.Label.String = 'Diagonal of critical eigenspace projector';
    else
        M = critical_modes{j};
        imagesc(ax, 1:N, 1:K, M);
        set(ax, 'YDir', 'normal', 'YTick', 1:K);

        clim(ax, [-signed_color_limit signed_color_limit]);
        colormap(ax, signed_map);
        xlabel(ax, 'Node index i');
        ylabel(ax, 'Layer k');
        title(ax, sprintf(['(%c) Normalized critical eigenvector\n' ...
            'r=%.1f, \\lambda_*=%.4f, IPR=%.4g'], ...
            char('a'+j+1), selected_ratios(j), ...
            lambda_star(j), IPR_full(j)), 'Interpreter', 'tex');
        cb = colorbar(ax, 'Location', 'eastoutside');
        cb.Label.String = 'Normalized eigenvector component';
    end
    box(ax, 'on');
end
exportgraphics(fig, out_img, 'Resolution', 300);

%% 5. 检查参考值并打印指定摘要
expected_lambda = [22 9.0105];
expected_sigma = [18.2353 18.3751];
expected_multiplicity = [4 1];
for j = 1:n_selected
    if abs(lambda_star(j)-expected_lambda(j)) > 0.5 || ...
            abs(sigma_c_spec(j)-expected_sigma(j)) > 0.2
        warning(['r=%.1f 的 lambda_star 或 sigma_c_spec 与预期参考值偏差较大；' ...
            '请检查输入谱文件与网络 realization。'], selected_ratios(j));
    end
    if multiplicity(j) ~= expected_multiplicity(j)
        warning('r=%.1f 的主导特征值重数为 %d，预期约为 %d。', ...
            selected_ratios(j), multiplicity(j), expected_multiplicity(j));
    end
    if abs(mu_max(j)) > 1e-5
        warning('r=%.1f 的临界增长率 |mu_max|=%.4g 未接近零。', ...
            selected_ratios(j), abs(mu_max(j)));
    end
    fprintf(['r=%.1f, lambda_star=%.8g, mu_max=%.4g, ' ...
        'multiplicity=%d'], selected_ratios(j), lambda_star(j), ...
        mu_max(j), multiplicity(j));
    if multiplicity(j) == 1
        fprintf(', IPR_full=%.8g, 1/(NK)=%.8g', ...
            IPR_full(j), IPR_reference);
    else
        fprintf(', IPR_full=not defined for a unique mode (degenerate)');
    end
    fprintf('\n');
end
fprintf('r=4.4 critical-subspace synchronous projection norm = %.4g\n', ...
    sync_projection_norm(1));
fprintf('lambda_opt=%.8g\n', lambda_opt);

save(out_mat, 'selected_ratios', 'lambda_all_selected', 'mu_all_selected', ...
    'lambda_cont', 'mu_cont', 'lambda_star', 'q_star', 'q_danger', ...
    'mu_max', 'sigma_c_spec', 'lambda_opt', 'sigma_min', 'critical_modes', ...
    'critical_projector_diag', 'IPR_full', 'IPR_reference', ...
    'sync_projection_norm', 'multiplicity', 'gap_left', 'gap_right', ...
    'degeneracy_tolerance', 'alpha', 'N', 'K', 'topology_file', 'input_file', ...
    'dominant_basis', '-v7');
fprintf('图片：%s\nMAT：%s\n', out_img, out_mat);

%% 局部函数
function map = blueWhiteRedMap(n_colors)
% 构造零值为浅奶白色的蓝白红发散色图。
half_count = floor((n_colors+1)/2);
red_count = n_colors-half_count+1;
neutral_color = [1.00 0.985 0.95];
blue_to_neutral = [linspace(0, neutral_color(1), half_count)', ...
    linspace(0, neutral_color(2), half_count)', ...
    linspace(1, neutral_color(3), half_count)'];
neutral_to_red = [linspace(neutral_color(1), 1, red_count)', ...
    linspace(neutral_color(2), 0, red_count)', ...
    linspace(neutral_color(3), 0, red_count)'];
map = [blue_to_neutral; neutral_to_red(2:end, :)];
end

function value = getRequiredField(S, candidates)
% 按候选名称读取必需字段。
for i = 1:numel(candidates)
    if isfield(S, candidates{i})
        value = S.(candidates{i});
        return;
    end
end
error('MAT 文件缺少必需字段，可接受名称：%s', strjoin(candidates, ', '));
end

function value = getOptionalField(S, candidates, default_value)
% 按候选名称读取可选字段；不存在时返回默认值。
value = default_value;
for i = 1:numel(candidates)
    if isfield(S, candidates{i})
        value = S.(candidates{i});
        return;
    end
end
end

function spectrum = getSpectrumAt(lambda_all, idx)
% 兼容 cell、列矩阵和行矩阵形式的逐 ratio 谱数据。
if iscell(lambda_all)
    spectrum = lambda_all{idx};
elseif size(lambda_all, 2) >= idx
    spectrum = lambda_all(:, idx);
elseif size(lambda_all, 1) >= idx
    spectrum = lambda_all(idx, :)';
else
    error('lambda_all 中找不到索引为 %d 的谱。', idx);
end
end

function vector = getEigenvectorAt(phi_all, idx)
% 兼容 cell 和按 ratio 分列保存的特征向量。
if iscell(phi_all)
    vector = phi_all{idx};
elseif size(phi_all, 2) >= idx
    vector = phi_all(:, idx);
else
    error('phi_star 中找不到索引为 %d 的特征向量。', idx);
end
end

function mu = growthRates(lambda_values, alpha, sigma, J)
% 计算各空间谱位置对应的最大实部线性增长率。
D = diag([1, sigma]);
lambda_values = lambda_values(:)';
mu = zeros(size(lambda_values));
for i = 1:numel(lambda_values)
    mu(i) = max(real(eig(J - alpha*lambda_values(i)*D)));
end
end

function topology_file = resolveTopologyFile(S, script_dir, topology_type, N, topology_param)
% 优先采用谱结果中记录的网络路径，否则按网络类型和参数定位已有实例。
topology_file = '';
if isfield(S, 'topology_file')
    saved_path = char(S.topology_file);
    if exist(saved_path, 'file')
        topology_file = saved_path;
    elseif exist(fullfile(script_dir, saved_path), 'file')
        topology_file = fullfile(script_dir, saved_path);
    end
end
if ~isempty(topology_file), return; end

topology_dir = fullfile(script_dir, 'results', 'topology', topology_type);
switch topology_type
    case 'ER'
        topology_file = fullfile(topology_dir, sprintf('N%d_p%.3f.mat', N, topology_param));
    case 'BA'
        topology_file = fullfile(topology_dir, sprintf('N%d_m%d.mat', N, topology_param));
    case 'WS'
        topology_file = fullfile(topology_dir, sprintf('N%d_pr%.2f.mat', N, topology_param));
    case 'SF'
        topology_file = fullfile(topology_dir, sprintf('N%d_g%.1f.mat', N, topology_param));
    otherwise
        error('不支持的网络类型：%s', topology_type);
end
end
