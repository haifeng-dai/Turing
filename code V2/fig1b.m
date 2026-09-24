%% fig1b.m: 专门展示 p=0.01 情况下的噪声(eta)对滞回阈值(sigma)的影响
clear; clc; close all;

% 1. 参数设置与数据加载
TOPO_TYPE = 'ER';
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta  = 0.1 * DYNA.alpha;
DYNA.noise = 0.01;
DYNA.T_END = 500;
DYNA.steps = 2;
DYNA.init_perturb = 0.1;
DYNA.sigma_min = 0;
DYNA.sigma_max = 35;
% 0.01 - 35 0.02 - 25

P_VAL = 0.01;           % 指定要测试和绘制的单个 p
ETA_LIST = linspace(0, 2, 201);

res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
data_path = fullfile(res_dir, sprintf('scan_eta_p%.3f_%s_N%d_K%d_a%.3f_b%.3f.mat', ...
    P_VAL, lower(TOPO_TYPE), DYNA.N, DYNA.K, DYNA.alpha, DYNA.beta));

FORCE_RERUN = false;   % 控制是否强制重跑: true = 强制重新计算覆盖旧数据, false = 存在则直接读取

% 缺失或强制重算时执行仿真
if FORCE_RERUN || ~exist(data_path, 'file')
    sigma_eta_connectivity_er(DYNA, TOPO_TYPE, P_VAL, ETA_LIST);
end

data = load(data_path);
eta_range = data.ETA_LIST;
sf_vec = data.sf_vec;
sb_vec = data.sb_vec;

% ========== 标注文字坐标配置 (易改参数) ==========
% Forward (蓝色) 标签坐标 - 文字、箭头起点、箭头终点完全独立
label_forward.txt_x = 0.47;          % 文字 x 坐标
label_forward.txt_y = 30;          % 文字 y 坐标
label_forward.arrow_start_x = 0.47;  % 箭头起点 x 坐标
label_forward.arrow_start_y = 29.5;  % 箭头起点 y 坐标
label_forward.arrow_end_x = 0.40;    % 箭头终点 x 坐标
label_forward.arrow_end_y = 25.7;    % 箭头终点 y 坐标

% Backward (橙色) 标签坐标 - 文字、箭头起点、箭头终点完全独立
label_backward.txt_x = 0.80;         % 文字 x 坐标
label_backward.txt_y = 8;          % 文字 y 坐标
label_backward.arrow_start_x = 0.8; % 箭头起点 x 坐标
label_backward.arrow_start_y = 9.0;  % 箭头起点 y 坐标
label_backward.arrow_end_x = 1;   % 箭头终点 x 坐标
label_backward.arrow_end_y = 12;   % 箭头终点 y 坐标

% 2. 绘图展示 (三相图风格: 区域填充 + 虚线边界)
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.43]);
ax = axes('Position', [0.18, 0.18, 0.75, 0.75]);
hold on;

% 动态确定 Y 轴边界 (保证填充美感)
max_val = max(sf_vec) * 1.15;
min_val = min(sb_vec) * 0.85;
fill_max_vec = ones(size(eta_range)) * max_val;
fill_min_vec = ones(size(eta_range)) * min_val;

% --- A. 区域填充 (按层级叠加) ---
% 1. 下层区域: 稳态/同步区 (Uniform State) - 采用浅青色
fill([eta_range, fliplr(eta_range)], [sb_vec, fill_min_vec], ...
    [0.9 0.95 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');

% 2. 中层区域: 双稳态/滞回区 (Bistable Region) - 采用浅灰色
fill([eta_range, fliplr(eta_range)], [sf_vec, fliplr(sb_vec)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.8, 'HandleVisibility', 'off');

% 3. 上层区域: 斑图/失稳区 (Pattern Region) - 采用浅桃色
fill([eta_range, fliplr(eta_range)], [fill_max_vec, fliplr(sf_vec)], ...
    [1.0 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');

% --- C. 边界曲线 (不画点，仅虚线) ---
p1 = plot(eta_range, sf_vec, '--', 'Color', [0 0.447 0.741], 'LineWidth', 2.5, ...
    'DisplayName', 'Forward Critical \sigma_c^f');

p2 = plot(eta_range, sb_vec, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 2.5, ...
    'DisplayName', 'Backward Critical \sigma_c^b');

% --- D. 在区域中添加文字标签 (采用用户手动调整的位置) ---
t1 = text(0.3, 12, 'Uniform', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

t2 = text(0.23, 24, 'Bistable', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

t3 = text(1.5, 26, 'Pattern', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

% 3. 属性修饰 (!!! 必须在此处提前设定范围，以供后续箭头坐标计算准确)
xlabel('$\eta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex');
grid off; box on;
xlim([0, 2]); ylim([min_val, 32]);
set(ax, 'FontSize', 18, 'Layer', 'top', 'TickLabelInterpreter', 'latex');

% --- F. 用 annotation 绘制真实的曲线箭头标注 (采用修正后的映射) ---
% 定义坐标转化公式 (基于最终的 xlim/ylim 运行)
ax_pos = get(ax, 'Position');
xl = xlim(ax); yl = ylim(ax);
d2f_x = @(xd) ax_pos(1) + (xd - xl(1))/(xl(2)-xl(1)) * ax_pos(3);
d2f_y = @(yd) ax_pos(2) + (yd - yl(1))/(yl(2)-yl(1)) * ax_pos(4);

% 1. Forward 蓝色曲线：完全解耦的文字、箭头起点、箭头终点
txt_f_x = label_forward.txt_x;
txt_f_y = label_forward.txt_y;
arrow_f_start_x = label_forward.arrow_start_x;
arrow_f_start_y = label_forward.arrow_start_y;
arrow_f_end_x = label_forward.arrow_end_x;
arrow_f_end_y = label_forward.arrow_end_y;
annotation(h, 'arrow', [d2f_x(arrow_f_start_x), d2f_x(arrow_f_end_x)], [d2f_y(arrow_f_start_y), d2f_y(arrow_f_end_y)], ...
    'Color', [0 0.447 0.741], 'LineWidth', 2, 'HeadStyle', 'vback2', 'HeadLength', 8);
tf = text(txt_f_x, txt_f_y, '$\sigma_c^f$', ...
    'HorizontalAlignment', 'left', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0 0.447 0.741]);

% 2. Backward 橙色曲线：完全解耦的文字、箭头起点、箭头终点
txt_b_x = label_backward.txt_x;
txt_b_y = label_backward.txt_y;
arrow_b_start_x = label_backward.arrow_start_x;
arrow_b_start_y = label_backward.arrow_start_y;
arrow_b_end_x = label_backward.arrow_end_x;
arrow_b_end_y = label_backward.arrow_end_y;
annotation(h, 'arrow', [d2f_x(arrow_b_start_x), d2f_x(arrow_b_end_x)], [d2f_y(arrow_b_start_y), d2f_y(arrow_b_end_y)], ...
    'Color', [0.85 0.325 0.098], 'LineWidth', 2, 'HeadStyle', 'vback2', 'HeadLength', 8);
tb = text(txt_b_x, txt_b_y, '$\sigma_c^b$', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.85 0.325 0.098]);

% 强制置顶
uistack([t1, t2, t3, tf, tb], 'top');

% 添加面板标签 (b) 在左上角外部
annotation(h, 'textbox', [0.05, 0.86, 0.08, 0.08], ...
    'String', '(b)', 'FontSize', 18, ...
    'BackgroundColor', 'none', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% 4. 导出图像
% plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'V2', 'manuscirpt', 'figures');
% if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
% out_img = fullfile(plots_dir, 'fig1b.eps');
% exportgraphics(h, out_img);
% fprintf('[DONE] 多区域精确标注的三相图已导出: %s\n', out_img);

% 测试环境 PNG 导出
test_plots_dir = fullfile(fileparts(mfilename('fullpath')), 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end
test_out_img = fullfile(test_plots_dir, 'fig1b.png');
exportgraphics(h, test_out_img, 'Resolution', 300);
fprintf('[DONE] 测试图片已保存: %s\n', test_out_img);
