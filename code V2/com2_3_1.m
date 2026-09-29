%% com2_3_1.m：主导特征值变化与 r=4.5 临界特征向量
% 读取已有谱分析结果，不生成网络，也不运行动力学仿真。
clear; clc; close all;

%% 1. 读取已有谱分析结果
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
lambda_star = getRequiredField(S, {'lambda_star', 'LAMBDA_STAR'});
lambda_all = getRequiredField(S, {'lambda_all', 'LAMBDA_ALL'});
phi_star = getRequiredField(S, {'phi_star', 'PHI_STAR'});
RATIO_LIST = RATIO_LIST(:)';
lambda_star = lambda_star(:)';
if isfield(S, 'N'), N = double(S.N); else, N = 200; end
if isfield(S, 'K'), K = double(S.K); else, K = 5; end
if isfield(S, 'ALPHA_FIXED')
    alpha = double(S.ALPHA_FIXED);
elseif isfield(S, 'alpha')
    alpha = double(S.alpha);
else
    alpha = 0.1;
end
if abs(alpha-0.1) > 1e-12
    error('输入结果中的 alpha=%.12g，不符合本图要求的 alpha=0.1。', alpha);
end
if numel(lambda_star) ~= numel(RATIO_LIST)
    error('lambda_star 与 RATIO_LIST 长度不一致。');
end

%% 2. 计算连续色散关系的最优谱位置
lambda_upper = (10/3) / alpha;
sigma_of_lambda = @(lambda) (110/3 + 4*alpha.*lambda) ./ ...
    ((alpha.*lambda) .* (10/3 - alpha.*lambda));
lambda_grid = linspace(lambda_upper*1e-8, lambda_upper*(1-1e-8), 100001);
sigma_grid = sigma_of_lambda(lambda_grid);
[~, grid_min_idx] = min(sigma_grid);
bracket_lo = lambda_grid(max(1, grid_min_idx-1));
bracket_hi = lambda_grid(min(numel(lambda_grid), grid_min_idx+1));
lambda_opt = fminbnd(sigma_of_lambda, bracket_lo, bracket_hi);

%% 3. 读取 r=4.5 的非简并主导特征向量
target_ratio = 4.5;
[ratio_error, target_idx] = min(abs(RATIO_LIST-target_ratio));
if ratio_error > 1e-7
    error('RATIO_LIST 中找不到 r=4.5。');
end
target_ratio = RATIO_LIST(target_idx);
target_lambda = lambda_star(target_idx);

spectrum = real(getSpectrumAt(lambda_all, target_idx));
spectrum = spectrum(:);
degeneracy_tolerance = max(1e-8, 1e-6*max(1, abs(target_lambda)));
multiplicity = sum(abs(spectrum-target_lambda) < degeneracy_tolerance);
if multiplicity ~= 1
    warning('r=%.1f 的主导特征值重数为 %d，不是预期的非简并情形。', ...
        target_ratio, multiplicity);
end

phi = getEigenvectorAt(phi_star, target_idx);
if numel(phi) ~= N*K || norm(phi) == 0
    error('r=%.1f 的 phi_star 长度无效。', target_ratio);
end
phi = real(phi(:));
phi = phi / norm(phi, 2);
[~, imax] = max(abs(phi));
if phi(imax) < 0
    phi = -phi;
end
M = reshape(phi, [N, K])';
IPR_full = sum(abs(phi).^4) / (sum(abs(phi).^2)^2);

%% 4. 绘制组合图
selected_ratios = [4.3 4.4 4.5 4.6];
selected_indices = zeros(size(selected_ratios));
for j = 1:numel(selected_ratios)
    [ratio_error, selected_indices(j)] = min(abs(RATIO_LIST-selected_ratios(j)));
    if ratio_error > 1e-7
        warning('RATIO_LIST 中没有 r=%.1f，未绘制该标记。', selected_ratios(j));
        selected_indices(j) = NaN;
    end
end

out_img = fullfile(fig_dir, 'com2_3_1.png');
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 80 1350 1050]);
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% 上图：复用 dominant_mode_switch_alpha010 的 lambda_star 曲线。
ax1 = nexttile(layout, 1);
plot(ax1, RATIO_LIST, lambda_star, '-o', 'LineWidth', 1.5, ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.0 0.447 0.741]);
hold(ax1, 'on');
yline(ax1, lambda_opt, '--', sprintf('\\lambda_{opt}=%.4g', lambda_opt), ...
    'LineWidth', 1.4, 'LabelHorizontalAlignment', 'left');
for j = 1:numel(selected_ratios)
    if isfinite(selected_indices(j))
        idx = selected_indices(j);
        plot(ax1, RATIO_LIST(idx), lambda_star(idx), 'kp', ...
            'MarkerFaceColor', [1.0 0.75 0.0], 'MarkerSize', 10, ...
            'HandleVisibility', 'off');
    end
end
xlim(ax1, [min(RATIO_LIST) max(RATIO_LIST)]);
xlabel(ax1, '\beta/\alpha', 'Interpreter', 'tex');
ylabel(ax1, '\lambda_*', 'Interpreter', 'tex');
title(ax1, '最危险超拉普拉斯模态的特征值');
grid(ax1, 'on'); box(ax1, 'on');

% 下图：复用 com2_3.png 右下的 r=4.5 signed eigenvector heatmap。
ax2 = nexttile(layout, 2);
imagesc(ax2, 1:N, 1:K, M);
set(ax2, 'YDir', 'normal', 'YTick', 1:K);
signed_map = blueWhiteRedMap(257);
signed_color_limit = 0.4;
clim(ax2, [-signed_color_limit signed_color_limit]);
colormap(ax2, signed_map);
xlabel(ax2, 'Node index i');
ylabel(ax2, 'Layer k');
title(ax2, sprintf(['Normalized critical eigenvector\n' ...
    'r=%.1f, \\lambda_*=%.4f, IPR=%.4g'], ...
    target_ratio, target_lambda, IPR_full), 'Interpreter', 'tex');
cb = colorbar(ax2, 'Location', 'eastoutside');
cb.Label.String = 'Normalized eigenvector component';
box(ax2, 'on');

exportgraphics(fig, out_img, 'Resolution', 300);
fprintf('图片：%s\n', out_img);

%% 局部函数
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

function spectrum = getSpectrumAt(lambda_all, idx)
% 兼容 cell、按 ratio 分列和按 ratio 分行的谱数据。
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
