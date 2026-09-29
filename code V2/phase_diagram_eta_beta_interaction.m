%% phase_diagram_eta_beta_interaction.m
% 联合扫描乘性噪声强度 eta 与层间/层内耦合比 beta/alpha 对 Turing 阈值的影响
clear; clc; close all;

project_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(project_dir, 'simulations'), fullfile(project_dir, 'networks'));

% 仅用于快速验证调用链；正常实验保持 false。
SMOKE_TEST = strcmp(getenv('TURING_ETA_RATIO_SMOKE_TEST'), '1');

%% 1. 固定参数与扫描变量
TOPO_TYPE = 'ER';
TOPO_PARAM = 0.03;
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.T_END = 500;
DYNA.steps = 2;
DYNA.init_perturb = 0.1;
DYNA.sigma_min = 0;
DYNA.sigma_max = 30;

etaList = 0:0.1:2;
ratioList = 0:0.1:10.0;

if SMOKE_TEST
    % 微型网格和短积分仅用于验证调用链，不代表实验数据。
    etaList = [0 0.1];
    ratioList = [0 0.25];
    DYNA.T_END = 0.2;
    DYNA.sigma_max = 2;
end

%% 2. 并行环境与固定网络/反向扫描种子
if isempty(gcp('nocreate'))
    if SMOKE_TEST
        parpool('Threads', 2);
    else
        parpool('Threads', 64);
    end
else
    pool = gcp;
    fprintf('[POOL] 使用现有并行池 (%d 个 worker)\n', pool.NumWorkers);
end

net_path = fullfile(project_dir, 'results', 'topology', 'ER', ...
    sprintf('N%d_p%.3f.mat', DYNA.N, TOPO_PARAM));
seed_path = fullfile(project_dir, 'results', ...
    sprintf('evolution_%s_N200_K5_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', ...
    lower(TOPO_TYPE)));
if ~exist(net_path, 'file'), error('找不到网络文件: %s', net_path); end
if ~exist(seed_path, 'file'), error('找不到全局斑图种子: %s', seed_path); end

net_data = load(net_path, 'nets');
seed_data = load(seed_path, 'Y');
L_intra_lib = net_data.nets;
y_universal_seed = seed_data.Y(end, :)';

adj_inter = ones(DYNA.K) - eye(DYNA.K);
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

%% 3. 按 eta 分行并行扫描，并逐行保存独立结果
numEta = length(etaList);
numRatio = length(ratioList);
sigmaF = nan(numEta, numRatio);
sigmaB = nan(numEta, numRatio);

result_dir = fullfile(project_dir, 'results');
fig_dir = fullfile(project_dir, 'fig');
if ~exist(result_dir, 'dir'), mkdir(result_dir); end
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
file_suffix = '';
if SMOKE_TEST, file_suffix = '_smoke_test'; end
file_stem = sprintf('eta_beta_interaction_%s_N%d_K%d_p%.3f_a%.3f%s', ...
    lower(TOPO_TYPE), DYNA.N, DYNA.K, TOPO_PARAM, DYNA.alpha, file_suffix);
result_path = fullfile(result_dir, [file_stem, '.mat']);

fprintf('[SIM] 按 eta 分行扫描：%d 行，每行 %d 个参数点；每行独立保存。\n', ...
    numEta, numRatio);
tStart = tic;
for etaIdx = 1:numEta
    eta_result_path = fullfile(result_dir, sprintf('%s_eta%.2f.mat', file_stem, etaList(etaIdx)));
    if exist(eta_result_path, 'file')
        eta_data = load(eta_result_path, 'etaIdx', 'sigmaFRow', 'sigmaBRow');
        if isfield(eta_data, 'sigmaFRow') && isfield(eta_data, 'sigmaBRow') && ...
                isequal(size(eta_data.sigmaFRow), [1, numRatio]) && ...
                isequal(size(eta_data.sigmaBRow), [1, numRatio])
            sigmaF(etaIdx, :) = eta_data.sigmaFRow;
            sigmaB(etaIdx, :) = eta_data.sigmaBRow;
            fprintf('[LOAD] eta=%.3f 已读取，跳过扫描。\n', etaList(etaIdx));
            continue;
        end
    end

    fprintf('[SIM] 开始 eta=%.3f (%d/%d)。\n', etaList(etaIdx), etaIdx, numEta);
    sigmaFRow = zeros(1, numRatio);
    sigmaBRow = zeros(1, numRatio);
    parfor ratioIdx = 1:numRatio
        dyna_local = DYNA;
        dyna_local.noise = etaList(etaIdx);
        dyna_local.beta = ratioList(ratioIdx) * dyna_local.alpha;
        dyna_local.L_intra = L_intra_lib;
        dyna_local.y_seed = y_universal_seed;

        thresholds = find_thresholds(TOPO_TYPE, TOPO_PARAM, dyna_local);
        sigmaFRow(ratioIdx) = thresholds(1);
        sigmaBRow(ratioIdx) = thresholds(2);
    end

    sigmaF(etaIdx, :) = sigmaFRow;
    sigmaB(etaIdx, :) = sigmaBRow;
    save(eta_result_path, 'etaList', 'ratioList', 'etaIdx', 'sigmaFRow', 'sigmaBRow', ...
        'TOPO_TYPE', 'TOPO_PARAM', 'DYNA');
    fprintf('[SAVE] eta=%.3f 已保存，累计耗时 %.2f 秒。\n', ...
        etaList(etaIdx), toc(tStart));
end
if any(isnan(sigmaF(:))) || any(isnan(sigmaB(:)))
    error('仍有 eta 文件缺失或无效，暂不生成总结果和图。');
end
fprintf('[DONE] 所有 eta 文件已读取或生成，耗时 %.2f 秒。\n', toc(tStart));

% 总结果只在当前 etaList 的所有行齐全后生成。
% 每台机器应使用不重叠的 etaList 子集，避免同时写同一个 eta 文件。

deltaSigma = sigmaF - sigmaB;

etaZeroIdx = find(abs(etaList) < eps, 1);
ratioZeroIdx = find(abs(ratioList) < eps, 1);
if isempty(etaZeroIdx) || isempty(ratioZeroIdx)
    error('interactionIndex 计算要求 etaList 和 ratioList 均包含 0。');
end
interactionIndex = sigmaF ...
    - repmat(sigmaF(etaZeroIdx, :), numEta, 1) ...
    - repmat(sigmaF(:, ratioZeroIdx), 1, numRatio) ...
    + sigmaF(etaZeroIdx, ratioZeroIdx);

assert(isequal(size(sigmaF), [numEta, numRatio]));
assert(isequal(size(sigmaB), [numEta, numRatio]));
assert(isequal(size(deltaSigma), [numEta, numRatio]));
assert(isequal(size(interactionIndex), [numEta, numRatio]));

%% 4. 保存完整结果
save(result_path, 'etaList', 'ratioList', 'sigmaF', 'sigmaB', ...
    'deltaSigma', 'interactionIndex', 'TOPO_TYPE', 'TOPO_PARAM', 'DYNA');
fprintf('[SAVE] 数据已保存: %s\n', result_path);

%% 5. 2 × 2 参数平面图（仅显示实际扫描点，不进行插值或平滑）
[ratioGrid, etaGrid] = meshgrid(ratioList, etaList);
fig = figure('Color', 'w', 'Units', 'normalized', ...
    'Position', [0.12, 0.12, 0.76, 0.72]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'loose', 'Padding', 'loose');

ax1 = nexttile;
scatter(ax1, ratioGrid(:), etaGrid(:), 150, sigmaF(:), 's', 'filled');
title(ax1, '(a) Forward threshold \sigma_c^f');
xlabel(ax1, '\beta/\alpha'); ylabel(ax1, '\eta');
colorbar(ax1); grid(ax1, 'on'); box(ax1, 'on');

ax2 = nexttile;
scatter(ax2, ratioGrid(:), etaGrid(:), 150, sigmaB(:), 's', 'filled');
title(ax2, '(b) Backward threshold \sigma_c^b');
xlabel(ax2, '\beta/\alpha'); ylabel(ax2, '\eta');
colorbar(ax2); grid(ax2, 'on'); box(ax2, 'on');

ax3 = nexttile;
scatter(ax3, ratioGrid(:), etaGrid(:), 150, deltaSigma(:), 's', 'filled');
title(ax3, '(c) Observed hysteresis width \sigma_c^f - \sigma_c^b');
xlabel(ax3, '\beta/\alpha'); ylabel(ax3, '\eta');
colorbar(ax3); grid(ax3, 'on'); box(ax3, 'on');

ax4 = nexttile;
scatter(ax4, ratioGrid(:), etaGrid(:), 150, interactionIndex(:), 's', 'filled', ...
    'MarkerEdgeColor', [0.4, 0.4, 0.4], 'LineWidth', 0.5);
title(ax4, '(d) Non-additive interaction I_f');
xlabel(ax4, '\beta/\alpha'); ylabel(ax4, '\eta');
colorbar(ax4); grid(ax4, 'on'); box(ax4, 'on');
interactionLimit = max(abs(interactionIndex(:)));
if interactionLimit == 0, interactionLimit = 1; end
clim(ax4, [-interactionLimit, interactionLimit]);
halfMap = 128;
blueWhiteRed = [linspace(0.15, 1, halfMap)', linspace(0.30, 1, halfMap)', ones(halfMap, 1); ...
    ones(halfMap, 1), linspace(1, 0.20, halfMap)', linspace(1, 0.15, halfMap)'];
colormap(ax4, blueWhiteRed);

xPad = 0.02 * max(ratioList(end) - ratioList(1), 1);
yPad = 0.02 * max(etaList(end) - etaList(1), 1);
xLimits = [ratioList(1) - xPad, ratioList(end) + xPad];
yLimits = [etaList(1) - yPad, etaList(end) + yPad];
set([ax1, ax2, ax3, ax4], 'FontSize', 11, 'Layer', 'top', ...
    'XLim', xLimits, 'YLim', yLimits);
set([ax1, ax2, ax3, ax4], 'XTick', ratioList, 'YTick', etaList);
yLabels = [ax1.YLabel, ax2.YLabel, ax3.YLabel, ax4.YLabel];
set(yLabels, 'Units', 'normalized');
for labelIdx = 1:numel(yLabels)
    labelPos = yLabels(labelIdx).Position;
    labelPos(1) = -0.14;
    yLabels(labelIdx).Position = labelPos;
end
figure_path = fullfile(fig_dir, [file_stem, '.png']);
exportgraphics(fig, figure_path, 'Resolution', 300);
fprintf('[DONE] 参数平面图已保存: %s\n', figure_path);

%% 6. Smoke test：确认项目默认随机数设置可复现相同轨迹
if SMOKE_TEST
    testCfg = DYNA;
    testCfg.noise = etaList(end);
    testCfg.beta = ratioList(end) * testCfg.alpha;
    testCfg.L_intra = L_intra_lib;
    testCfg.sigma = 1;
    testCfg.y0 = [];
    testCfg.detect_convergence = false;

    rng(0, 'twister');
    [~, testY1] = solve_multiplex(testCfg);
    rng(0, 'twister');
    [~, testY2] = solve_multiplex(testCfg);
    assert(isequal(testY1, testY2), '项目默认随机数设置未能复现相同轨迹。');
    assert(isgraphics(fig), '2 × 2 figure 未能正常生成。');
    fprintf('[SMOKE] 并行、项目默认随机数可复现、矩阵尺寸及 2 × 2 绘图检查通过。\n');
end
