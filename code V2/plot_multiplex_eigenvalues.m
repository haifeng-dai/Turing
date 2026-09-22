%% plot_multiplex_eigenvalues.m: 分析多层网络特征值随层间耦合强度(beta/alpha)的变化
clear; clc; close all;

%% 1. 基本参数设置
N = 200;
K = 5;
p = 0.01; % 使用标准连接概率
ratio_list = 0:0.01:1; % beta/alpha 的比率范围，扫描至 25 以观察饱和趋势

%% 2. 加载拓扑结构
topo_file = fullfile('results', 'topology', 'ER', sprintf('N%d_p%.3f.mat', N, p));
if ~exist(topo_file, 'file')
    error('未找到拓扑文件: %s', topo_file);
end
data = load(topo_file);
nets = data.nets(1:K); % 【核心修复】显式截取前 K 层，确保与 L_inter 的维度 NK 匹配

% 构造全系统层内拉普拉斯 (Block Diagonal)
L_intra_all = blkdiag(nets{:});

% 构造全连通层间拉普拉斯 (K-1 同步分量)
adj_inter = ones(K) - eye(K);
L_inter_base = diag(sum(adj_inter, 2)) - adj_inter;
L_inter_all = kron(L_inter_base, speye(N));

%% 3. 循环计算特征值与 Gamma 演化
num_r = length(ratio_list);
lam2 = zeros(num_r, 1);
lamNK = zeros(num_r, 1);
gamma_vec = zeros(num_r, 1); % 稳定性参数 Gamma (Stochastic Stability)

% 动力学分析参数
SIGMA_FIXED = 16.0; 
ETA_FIXED   = 0.3;   % 【关键】乘性噪声强度，对 Stability 有显著偏置作用
J     = [10/3, -5; 10, -4]; % 反应 Jacobian
D_eff = diag([1, SIGMA_FIXED]);

fprintf('[CALC] 正在扫描 %d 组耦合强度比率 (beta/alpha)...\n', num_r);
simStart = tic;

parfor i = 1:num_r
    beta_val = ratio_list(i);
    % 合成全系统超拉普拉斯算子
    L_total = L_intra_all + beta_val * L_inter_all;

    % 1. 特征谱计算
    ev = sort(real(eig(L_total)));
    lam2(i)  = ev(2);
    lamNK(i) = ev(end);

    % 2. 随机稳定性分析 (Stochastic Gamma)
    % 公式修正：M = J - mu*D + (eta^2/2)*(mu*D)^2
    NK_total = length(ev);
    max_re_eig_total = -inf;
    for k = 1:NK_total
        mu = ev(k);
        % 每个模式下的等效矩阵 (计入二阶噪声修正项)
        Mi = J - mu * D_eff + (ETA_FIXED^2 / 2) * (mu^2 * (D_eff^2));
        max_re_eig_total = max(max_re_eig_total, max(real(eig(Mi))));
    end
    gamma_vec(i) = -max_re_eig_total; 
end

fprintf('[DONE] 特征值计算完成，耗时 %.2f 秒。\n', toc(simStart));

%% 4. 数据可视化 (双轴展示)
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.45, 0.45]);
ax = axes('Position', [0.15, 0.15, 0.70, 0.75]); hold on;

% --- A. 左轴：绘制 lambda_2 (反映结构融合) ---
yyaxis left
p1 = plot(ratio_list, lam2, '-', 'Color', [0 0.447 0.741], 'LineWidth', 2.5);
ylabel('\lambda_2 (Algebraic Connectivity)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'YColor', [0 0.447 0.741]);

% --- B. 右轴：绘制 lambda_max (反映谱能量上限) ---
yyaxis right
p2 = plot(ratio_list, lamNK, '-', 'Color', [0.85 0.325 0.098], 'LineWidth', 2.5);
ylabel('\lambda_{NK} (Spectral Radius)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'YColor', [0.85 0.325 0.098]);

% --- C. 属性精修 ---
xlabel('Inter-layer Coupling Ratio \beta / \alpha', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Eigenvalue Evolution of Multiplex ER Network (p=%.3f)', p), 'FontSize', 14);
grid on; box on;
set(ax, 'FontSize', 13, 'LineWidth', 1.2);

% 标注图例
legend([p1, p2], {'\lambda_2 (left axis)', '\lambda_{max} (right axis)'}, ...
    'Location', 'northwest', 'FontSize', 12, 'EdgeColor', 'none');

% 导出图像
plots_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
out_img = fullfile(plots_dir, 'multiplex_eigenvalue_evolution.png');
exportgraphics(h, out_img, 'Resolution', 300);

fprintf('[DONE] 第一张演化图已生成并存至: %s\n', out_img);

%% 5. 数据可视化 (图 2: 特征值差值与比值)
% --- A. 准备数据并过滤奇异值 ---
lam_diff  = lamNK - lam2;
lam_ratio = lamNK ./ (lam2 + eps);
% 过滤掉 lambda2=0 的点 (beta=0)，以防对数坐标轴崩溃或缩放失败
valid_idx = (ratio_list > 1e-4);

h2 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.1, 0.45, 0.45]);
ax2 = axes('Position', [0.15, 0.15, 0.70, 0.75]); hold on;

% --- B. 左轴：绘制特征值差值 ---
yyaxis left
p3 = plot(ratio_list, lam_diff, '-', 'Color', [0.49, 0.18, 0.56], 'LineWidth', 2.5);
ylabel('\lambda_{NK} - \lambda_2 (Difference)', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0, max(lam_diff)*1.1]); % 自动适配范围，保留顶部空间
set(gca, 'YColor', [0.49, 0.18, 0.56]);

% --- C. 右轴：使用对数坐标绘制比值 (Log Scale) ---
yyaxis right
p4 = plot(ratio_list(valid_idx), lam_ratio(valid_idx), '-', 'Color', [0.47, 0.67, 0.19], 'LineWidth', 2.5);
ylabel('\lambda_{NK} / \lambda_2 (Ratio, Log_1_0 Scale)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'YScale', 'log', 'YColor', [0.47, 0.67, 0.19]); % 【关键修复】对数轴

% 设置一个合理的比值范围：从 1 到最大值的 10 倍，便于观察
max_r = max(lam_ratio(valid_idx));
ylim([1, 10^(ceil(log10(max_r))+1)]);

% --- D. 属性精修 ---
xlabel('Inter-layer Coupling Ratio \beta / \alpha', 'FontSize', 14, 'FontWeight', 'bold');
title('Eigenvalue Span and Relative Connectivity Ratio', 'FontSize', 14);
grid on; box on;
set(ax2, 'FontSize', 13, 'LineWidth', 1.2);

% 标注图例
legend([p3, p4], {'Difference (Linear)', 'Ratio (Log Scale)'}, ...
    'Location', 'northeast', 'FontSize', 12, 'EdgeColor', 'none', 'FontWeight', 'bold');

% 导出第二张图像
out_img2 = fullfile(plots_dir, 'multiplex_eigenvalue_stats.png');
exportgraphics(h2, out_img2, 'Resolution', 300);

fprintf('[DONE] 优化后的统计图已生成并存至: %s\n', out_img2);

%% 6. 数据可视化 (图 3: Gamma 稳定性分析)
h3 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.05, 0.45, 0.45]);
ax3 = axes('Position', [0.15, 0.15, 0.75, 0.75]); hold on;

plot(ratio_list, gamma_vec, 'LineWidth', 2.5, 'Color', [0.64, 0.08, 0.18]);
yline(0, '--k', 'LineWidth', 1.5); % 临界稳定性线 (\gamma=0)

xlabel('Inter-layer Coupling Ratio \beta / \alpha', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Stability Parameter \gamma', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Stochastic Stability \\gamma vs. Coupling (\\sigma=%.1f, \\eta=%.1f)', SIGMA_FIXED, ETA_FIXED), 'FontSize', 14);
grid on; box on;
set(ax3, 'FontSize', 13, 'LineWidth', 1.2);

% 快速说明：gamma < 0 时系统不稳定，支持斑图形成
text(0.6*max(ratio_list), 0.8*min(gamma_vec), 'Pattern Region (\gamma < 0)', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.64, 0.08, 0.18]);

% 导出第三张图像
out_img3 = fullfile(plots_dir, 'multiplex_gamma_evolution.png');
exportgraphics(h3, out_img3, 'Resolution', 300);

fprintf('[DONE] 稳定性参数演化图已生成并存至: %s\n', out_img3);
