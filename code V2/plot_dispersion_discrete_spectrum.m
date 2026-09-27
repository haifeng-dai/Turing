%% plot_dispersion_discrete_spectrum.m: 色散关系与离散超拉普拉斯谱点
clear; clc; close all;

%% 1. 路径与参数配置
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(script_dir, 'simulations'), fullfile(script_dir, 'networks'));

TOPO_TYPE = 'ER';
TOPO_PARAM = 0.03;  % ER：连接概率 p；WS：重连概率；BA：参数 m；SF：幂律指数 gamma

cfg.N = 200;
cfg.K = 5;
cfg.alpha = 0.05;
SIGMA = 20;  % Fig. 2 参数范围内的一个固定谱截面，可按需要调整
NOISE_REFERENCE = 0.01;  % Fig. 2 的噪声参数；不进入确定性线性增长率
RATIO_LIST = [0.1, 6.0, 10.0];  % 与 Fig. 2 组合图选取的 beta / alpha 一致
RATIO_LABELS = {'弱耦合', '转折附近', '强耦合'};

%% 2. 加载层内网络并构造层间网络
topology_dir = fullfile(script_dir, 'results', 'topology', upper(TOPO_TYPE));
switch upper(TOPO_TYPE)
    case 'ER'
        topology_file = fullfile(topology_dir, ...
            sprintf('N%d_p%.3f.mat', cfg.N, TOPO_PARAM));
    case 'WS'
        topology_file = fullfile(topology_dir, ...
            sprintf('N%d_pr%.2f.mat', cfg.N, TOPO_PARAM));
    case 'BA'
        topology_file = fullfile(topology_dir, ...
            sprintf('N%d_m%d.mat', cfg.N, TOPO_PARAM));
    case 'SF'
        topology_file = fullfile(topology_dir, ...
            sprintf('N%d_g%.1f.mat', cfg.N, TOPO_PARAM));
    otherwise
        error('不支持的拓扑类型：%s', TOPO_TYPE);
end

if ~exist(topology_file, 'file')
    error('找不到拓扑文件：%s', topology_file);
end

topology_data = load(topology_file, 'nets');
if numel(topology_data.nets) < cfg.K
    error('拓扑文件中的层数不足：需要 %d 层，实际只有 %d 层。', ...
        cfg.K, numel(topology_data.nets));
end

L_intra = cell(cfg.K, 1);
for k = 1:cfg.K
    Lk = topology_data.nets{k};
    % 当前拓扑库保存的是拉普拉斯矩阵，而不是邻接矩阵。
    if max(abs(full(sum(Lk, 2)))) > 1e-8 || any(diag(Lk) < -1e-12)
        error('nets{%d} 不是当前脚本要求的拉普拉斯矩阵，请先从邻接矩阵构造 D-A。', k);
    end
    L_intra{k} = sparse(Lk);
end
adj_inter = ones(cfg.K) - eye(cfg.K);
L_inter = diag(sum(adj_inter, 2)) - adj_inter;

%% 3. 计算反应项线性化矩阵与离散增长率
% 稳态为 (u*, v*) = (5, 10)，反应项雅可比矩阵如下。
J = [10 / 3, -5; 10, -4];
L_intra_all = blkdiag(L_intra{:});
L_inter_big = kron(L_inter, speye(cfg.N));

lambda_all = cell(size(RATIO_LIST));
mu_all = cell(size(RATIO_LIST));
lambda_unstable = cell(size(RATIO_LIST));
mu_unstable = cell(size(RATIO_LIST));
lambda_max_growth = zeros(size(RATIO_LIST));
mu_max_growth = zeros(size(RATIO_LIST));

for i = 1:numel(RATIO_LIST)
    ratio = RATIO_LIST(i);
    % 采用论文记号：L = L^L + (beta/alpha) L^I。
    L_super = L_intra_all + ratio * L_inter_big;
    lambda = sort(real(eig(full(L_super))));
    mu = zeros(size(lambda));

    for q = 1:numel(lambda)
        mode_matrix = J - cfg.alpha * lambda(q) * diag([1, SIGMA]);
        mu(q) = max(real(eig(mode_matrix)));
    end

    lambda_all{i} = lambda;
    mu_all{i} = mu;
    unstable_mask = mu > 0;
    lambda_unstable{i} = lambda(unstable_mask);
    mu_unstable{i} = mu(unstable_mask);

    [mu_max_growth(i), max_idx] = max(mu);
    lambda_max_growth(i) = lambda(max_idx);
end

%% 4. 在同一个坐标轴中叠加三组色散关系和离散谱点
lambda_max = max(cellfun(@max, lambda_all));
lambda_grid = linspace(0, 1.03 * lambda_max, 1000);
mu_grid = zeros(size(lambda_grid));

for q = 1:numel(lambda_grid)
    mode_matrix = J - cfg.alpha * lambda_grid(q) * diag([1, SIGMA]);
    mu_grid(q) = max(real(eig(mode_matrix)));
end

mu_limit = max(mu_grid);
mu_limit = max(mu_limit, max(cellfun(@max, mu_all)));
mu_min = min(mu_grid);
mu_min = min(mu_min, min(cellfun(@min, mu_all)));
mu_padding = 0.08 * max(mu_limit - mu_min, eps);

h = figure('Visible', 'off', 'Color', 'w', ...
    'Name', '色散关系与离散超拉普拉斯谱点', ...
    'Units', 'normalized', 'Position', [0.12, 0.16, 0.7, 0.68]);
ax = axes(h);
hold(ax, 'on');
colors = lines(numel(RATIO_LIST));

% 三组参数共用一条确定性连续色散曲线。
plot(ax, lambda_grid, mu_grid, 'k-', 'LineWidth', 2.0, ...
    'DisplayName', '共同连续色散关系');

for i = 1:numel(RATIO_LIST)
    ratio_name = sprintf('\\beta/\\alpha=%.3g（%s）', ...
        RATIO_LIST(i), RATIO_LABELS{i});

    % 同色小点表示该 beta/alpha 下的全部离散谱点。
    scatter(ax, lambda_all{i}, mu_all{i}, 12, colors(i, :), ...
        'filled', 'HandleVisibility', 'off');

    % 失稳点加大并加黑边，便于在三组点叠加后识别。
    if ~isempty(lambda_unstable{i})
        scatter(ax, lambda_unstable{i}, mu_unstable{i}, 25, colors(i, :), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.35, ...
            'HandleVisibility', 'off');
    end

    % 每组的最危险谱点用五角星标出。
    scatter(ax, lambda_max_growth(i), mu_max_growth(i), 65, colors(i, :), ...
        'p', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');

    % 用图例代理对象说明三组离散谱点的颜色。
    plot(ax, nan, nan, 'o', 'Color', colors(i, :), ...
        'MarkerFaceColor', colors(i, :), 'LineWidth', 1.0, ...
        'DisplayName', [ratio_name, ' 离散谱点']);
end

yline(ax, 0, 'k-', 'LineWidth', 0.9, 'DisplayName', '\mu=0');
xlabel(ax, '\lambda_q', 'FontSize', 15, 'Interpreter', 'tex');
ylabel(ax, '\mu_q', 'FontSize', 15, 'Interpreter', 'tex');
title(ax, sprintf('共同色散关系与离散谱点（\\sigma=%.6g，ER，p=%.3f，N=%d，K=%d，\\eta=%.2g）', ...
    SIGMA, TOPO_PARAM, cfg.N, cfg.K, NOISE_REFERENCE), 'FontSize', 15);
xlim(ax, [0, 1.03 * lambda_max]);
ylim(ax, [mu_min - mu_padding, mu_limit + mu_padding]);
grid(ax, 'on');
box(ax, 'on');
set(ax, 'FontSize', 12);
legend(ax, 'Location', 'best', 'FontSize', 10);

% 在图的右上角集中列出三组失稳点数量和最大增长率。
summary_text = cell(numel(RATIO_LIST), 1);
for i = 1:numel(RATIO_LIST)
    summary_text{i} = sprintf('%.3g：%d/%d 个失稳，最大 \\mu=%.3g', ...
        RATIO_LIST(i), numel(lambda_unstable{i}), numel(lambda_all{i}), mu_max_growth(i));
end
text(ax, 0.98, 0.98, summary_text, 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'BackgroundColor', 'w', 'Margin', 4);

%% 5. 保存图像和谱数据
fig_dir = fullfile(script_dir, 'fig');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
out_img = fullfile(fig_dir, sprintf( ...
    'dispersion_discrete_spectrum_%s_N%d_K%d_p%.3f_s%.6g.png', ...
    lower(TOPO_TYPE), cfg.N, cfg.K, TOPO_PARAM, SIGMA));
exportgraphics(h, out_img, 'Resolution', 300);

data_dir = fullfile(script_dir, 'results');
out_data = fullfile(data_dir, sprintf( ...
    'dispersion_discrete_spectrum_%s_N%d_K%d_p%.3f_s%.6g.mat', ...
    lower(TOPO_TYPE), cfg.N, cfg.K, TOPO_PARAM, SIGMA));
save(out_data, 'SIGMA', 'NOISE_REFERENCE', 'RATIO_LIST', 'lambda_all', 'mu_all', ...
    'lambda_unstable', 'mu_unstable', 'lambda_max_growth', 'mu_max_growth');

fprintf('色散关系图已保存：%s\n', out_img);
fprintf('谱数据已保存：%s\n', out_data);
for i = 1:numel(RATIO_LIST)
    fprintf('beta/alpha=%.6g：失稳点 %d/%d，最危险点 lambda=%.8g，mu=%.8g\n', ...
        RATIO_LIST(i), numel(lambda_unstable{i}), numel(lambda_all{i}), ...
        lambda_max_growth(i), mu_max_growth(i));
end
