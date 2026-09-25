%% fig_hetero_comparison.m: 分别独立导出 BA、WS 与 BA-WS-ER 异构网络的 (beta, sigma) 三相图
% 完全对齐图 2c (fig2abc.m) 风格: 独立的单张图窗、相对坐标箭头与文字解耦映射、单张高清导出
clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'simulations'), fullfile(script_dir, 'networks'));

RUN_SIMULATION = false;  % 控制是否运行仿真 (true: 自动计算并保存, false: 直接读取已有数据绘图)
FORCE_RERUN    = false;  % 是否强制重新计算已存在的数据

%% 1. 全局物理与动力学配置 (完全对齐图 2c)
N = 200;                % 节点数
K = 3;                  % 层数设为 3 层 (以支持 BA-WS-ER 异构对称对比)
DYNA.N = N;
DYNA.K = K;
DYNA.alpha = 0.10;      % 对齐图 2c 的 alpha = 0.10
DYNA.noise = 0.01;      % 微弱背景扰动
DYNA.T_END = 500;
DYNA.steps = 2;
DYNA.init_perturb = 0.1;
DYNA.sigma_min = 0;
DYNA.sigma_max = 80;

% 横轴参数：层间耦合比 ratio = beta / alpha
RATIO_LIST = 0:0.25:10.0; % 41 个采样点
num_r = length(RATIO_LIST);

res_dir = fullfile(script_dir, 'results');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end
test_plots_dir = fullfile(script_dir, 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end

%% 2. 定义三种对比网络案例与各自独立的绘图标注坐标配置
% ========== 为每个网络独立配置标注文字和箭头坐标 (0~1 归一化) ==========
% 案例 1: Pure BA (3层全 BA, m=4)
CASES(1).name = 'BA';
CASES(1).label = 'Pure BA Multiplex';
CASES(1).tag = 'pure_ba';
CASES(1).suffix = 'ba';
CONFIG(1).forward.txt_x_norm = 0.05;
CONFIG(1).forward.txt_y_norm = 0.65;
CONFIG(1).forward.arrow_start_x_norm = 0.08;
CONFIG(1).forward.arrow_start_y_norm = 0.60;
CONFIG(1).forward.arrow_end_x_norm = 0.12;
CONFIG(1).forward.arrow_end_y_norm = 0.45;
CONFIG(1).backward.txt_x_norm = 0.08;
CONFIG(1).backward.txt_y_norm = 0.05;
CONFIG(1).backward.arrow_start_x_norm = 0.13;
CONFIG(1).backward.arrow_start_y_norm = 0.08;
CONFIG(1).backward.arrow_end_x_norm = 0.20;
CONFIG(1).backward.arrow_end_y_norm = 0.20;
CONFIG(1).region_x_norms = [0.55, 0.70, 0.55]; % Uniform, Bistable, Pattern
CONFIG(1).region_y_norms = [0.10, 0.50, 0.92];

% 案例 2: Pure WS (3层全 WS, K=6, pr=0.10)
CASES(2).name = 'WS';
CASES(2).label = 'Pure WS Multiplex';
CASES(2).tag = 'pure_ws';
CASES(2).suffix = 'ws';
CONFIG(2).forward.txt_x_norm = 0.05;
CONFIG(2).forward.txt_y_norm = 0.65;
CONFIG(2).forward.arrow_start_x_norm = 0.08;
CONFIG(2).forward.arrow_start_y_norm = 0.60;
CONFIG(2).forward.arrow_end_x_norm = 0.12;
CONFIG(2).forward.arrow_end_y_norm = 0.45;
CONFIG(2).backward.txt_x_norm = 0.08;
CONFIG(2).backward.txt_y_norm = 0.05;
CONFIG(2).backward.arrow_start_x_norm = 0.13;
CONFIG(2).backward.arrow_start_y_norm = 0.08;
CONFIG(2).backward.arrow_end_x_norm = 0.20;
CONFIG(2).backward.arrow_end_y_norm = 0.20;
CONFIG(2).region_x_norms = [0.55, 0.70, 0.55];
CONFIG(2).region_y_norms = [0.10, 0.50, 0.92];

% 案例 3: BA-WS-ER Heterogeneous (Layer 1=BA, Layer 2=WS, Layer 3=ER)
CASES(3).name = 'HETERO';
CASES(3).label = 'BA-WS-ER Heterogeneous Multiplex';
CASES(3).tag = 'hetero_ba_ws_er';
CASES(3).suffix = 'hetero';
CONFIG(3).forward.txt_x_norm = 0.05;
CONFIG(3).forward.txt_y_norm = 0.65;
CONFIG(3).forward.arrow_start_x_norm = 0.08;
CONFIG(3).forward.arrow_start_y_norm = 0.60;
CONFIG(3).forward.arrow_end_x_norm = 0.12;
CONFIG(3).forward.arrow_end_y_norm = 0.45;
CONFIG(3).backward.txt_x_norm = 0.08;
CONFIG(3).backward.txt_y_norm = 0.05;
CONFIG(3).backward.arrow_start_x_norm = 0.13;
CONFIG(3).backward.arrow_start_y_norm = 0.08;
CONFIG(3).backward.arrow_end_x_norm = 0.20;
CONFIG(3).backward.arrow_end_y_norm = 0.20;
CONFIG(3).region_x_norms = [0.55, 0.70, 0.55];
CONFIG(3).region_y_norms = [0.10, 0.50, 0.92];

%% 3. 并行仿真区 (独立计算与落盘)
seed_file = fullfile(res_dir, 'evolution_er_N200_K5_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat');
if ~exist(seed_file, 'file')
    error('未找到基准种子文件 %s，请先运行 pattern_evolution.m 生成参考种子。', seed_file);
end
tmp_seed = load(seed_file);
y_universal_seed = tmp_seed.Y(end, :)';

% 预加载三种底层网络
tmp_ba = load(fullfile(res_dir, 'topology', 'BA', 'N200_m4.mat'));
tmp_ws = load(fullfile(res_dir, 'topology', 'WS', 'N200_K6_pr0.10.mat'));
tmp_er = load(fullfile(res_dir, 'topology', 'ER', 'N200_p0.030.mat'));

NET_SETS(1).L = tmp_ba.nets(1:K); % Pure BA
NET_SETS(2).L = tmp_ws.nets(1:K); % Pure WS
NET_SETS(3).L = {tmp_ba.nets{1}, tmp_ws.nets{1}, tmp_er.nets{1}}; % BA-WS-ER

% 准备层间拉普拉斯算子
adj_inter = (ones(K) - eye(K));
DYNA.L_inter = diag(sum(adj_inter, 2)) - adj_inter;

for c_idx = 1:3
    case_info = CASES(c_idx);
    mat_name = sprintf('scan_beta_phase_%s_N%d_K%d_a%.3f.mat', case_info.tag, N, K, DYNA.alpha);
    data_path = fullfile(res_dir, mat_name);

    if RUN_SIMULATION && (FORCE_RERUN || ~exist(data_path, 'file'))
        fprintf('\n============================================================\n');
        fprintf('[SIM] 正在计算: %s (任务 %d/3)...\n', case_info.label, c_idx);
        fprintf('============================================================\n');

        if isempty(gcp('nocreate')), parpool('Threads', 64); end

        L_current = NET_SETS(c_idx).L;
        sf_vec = zeros(1, num_r);
        sb_vec = zeros(1, num_r);
        simStart = tic;

        parfor r_i = 1:num_r
            cur_ratio = RATIO_LIST(r_i);
            dyna_local = DYNA;
            dyna_local.beta = dyna_local.alpha * cur_ratio;
            dyna_local.L_intra = L_current;
            dyna_local.y_seed = y_universal_seed;

            res = find_thresholds('ER', 0.030, dyna_local);
            sf_vec(r_i) = res(1);
            sb_vec(r_i) = res(2);
            fprintf('  [%s] ratio=%.2f 完成 -> sf=%.2f, sb=%.2f\n', ...
                case_info.name, cur_ratio, res(1), res(2));
        end

        fprintf('[DONE] %s 计算完成，耗时: %.2f 秒。\n', case_info.label, toc(simStart));
        save(data_path, 'sf_vec', 'sb_vec', 'RATIO_LIST', 'DYNA', 'case_info');
    end
end

%% 4. 【独立绘图与导出区】(与 fig2abc.m 保持 100% 一致的单独图窗与导出)
fprintf('\n[VIS] 开始逐个独立绘制并导出相图...\n');

for c_idx = 1:3
    case_info = CASES(c_idx);
    cfg = CONFIG(c_idx);
    mat_name = sprintf('scan_beta_phase_%s_N%d_K%d_a%.3f.mat', case_info.tag, N, K, DYNA.alpha);
    data_path = fullfile(res_dir, mat_name);

    if ~exist(data_path, 'file')
        warning('未找到数据文件: %s，跳过该图绘制。', mat_name);
        continue;
    end

    data = load(data_path);
    x_data = data.RATIO_LIST;
    sf_vec = data.sf_vec;
    sb_vec = data.sb_vec;

    % 避免填充时出现负宽度
    swap_mask = sf_vec < sb_vec;
    if any(swap_mask)
        tmp = sf_vec(swap_mask);
        sf_vec(swap_mask) = sb_vec(swap_mask);
        sb_vec(swap_mask) = tmp;
    end

    % 为每个网络创建独立的单一图窗 (尺寸完全对齐 fig2abc.m)
    h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.46, 0.52]);
    ax = axes('Position', [0.14, 0.20, 0.80, 0.74]); 
    hold(ax, 'on');

    max_val = max(sf_vec) * 1.12;
    min_val = max(0, min(sb_vec) * 0.88);
    fill_max_vec = ones(size(x_data)) * max_val;
    fill_min_vec = ones(size(x_data)) * min_val;

    % 1. 下层: Uniform (浅青)
    fill([x_data, fliplr(x_data)], [sb_vec, fill_min_vec], ...
        [0.90 0.95 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');

    % 2. 中层: Bistable (浅灰)
    fill([x_data, fliplr(x_data)], [sf_vec, fliplr(sb_vec)], ...
        [0.90 0.90 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');

    % 3. 上层: Pattern (浅桃粉)
    fill([x_data, fliplr(x_data)], [fill_max_vec, fliplr(sf_vec)], ...
        [1.00 0.95 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');

    % 4. 临界边界虚线
    plot(x_data, sf_vec, '--', 'Color', [0.00 0.447 0.741], 'LineWidth', 2.5, ...
        'DisplayName', 'Forward Critical \sigma_c^f');
    plot(x_data, sb_vec, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 2.5, ...
        'DisplayName', 'Backward Critical \sigma_c^b');

    % 5. 坐标轴修饰
    xlabel('$\beta / \alpha$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$\sigma$', 'FontSize', 18, 'Interpreter', 'latex');
    box on; grid off;
    xlim([0, 10]);
    ylim([min_val, max_val]);
    set(ax, 'FontSize', 18, 'Layer', 'top', 'TickLabelInterpreter', 'latex');

    % 6. 坐标映射 (供区域文字和箭头精确定位)
    ax_pos = get(ax, 'Position');
    xl = xlim(ax);
    yl = ylim(ax);
    norm2data_x = @(xn) xl(1) + xn * (xl(2) - xl(1));
    norm2data_y = @(yn) yl(1) + yn * (yl(2) - yl(1));
    data2fig_x = @(xd) ax_pos(1) + (xd - xl(1)) / (xl(2) - xl(1)) * ax_pos(3);
    data2fig_y = @(yd) ax_pos(2) + (yd - yl(1)) / (yl(2) - yl(1)) * ax_pos(4);
    clamp01 = @(v) max(0, min(1, v));

    % 区域文字 - 使用归一化映射
    t1 = text(norm2data_x(cfg.region_x_norms(1)), min_val + cfg.region_y_norms(1)*(max_val-min_val), 'Uniform', ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
    t2 = text(norm2data_x(cfg.region_x_norms(2)), min_val + cfg.region_y_norms(2)*(max_val-min_val), 'Bistable', ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
    t3 = text(norm2data_x(cfg.region_x_norms(3)), min_val + cfg.region_y_norms(3)*(max_val-min_val), 'Pattern', ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

    % Forward 箭头与标注 (蓝色)
    txt_f_data_x = norm2data_x(cfg.forward.txt_x_norm);
    txt_f_data_y = norm2data_y(cfg.forward.txt_y_norm);
    arr_f_s_x = clamp01(data2fig_x(norm2data_x(cfg.forward.arrow_start_x_norm)));
    arr_f_s_y = clamp01(data2fig_y(norm2data_y(cfg.forward.arrow_start_y_norm)));
    arr_f_e_x = clamp01(data2fig_x(norm2data_x(cfg.forward.arrow_end_x_norm)));
    arr_f_e_y = clamp01(data2fig_y(norm2data_y(cfg.forward.arrow_end_y_norm)));
    annotation(h, 'arrow', [arr_f_s_x, arr_f_e_x], [arr_f_s_y, arr_f_e_y], ...
        'Color', [0.00 0.447 0.741], 'LineWidth', 2, 'HeadStyle', 'vback2', 'HeadLength', 8);
    tf = text(txt_f_data_x, txt_f_data_y, '$\sigma_c^f$', ...
        'HorizontalAlignment', 'left', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.00 0.447 0.741]);

    % Backward 箭头与标注 (橙色)
    txt_b_data_x = norm2data_x(cfg.backward.txt_x_norm);
    txt_b_data_y = norm2data_y(cfg.backward.txt_y_norm);
    arr_b_s_x = clamp01(data2fig_x(norm2data_x(cfg.backward.arrow_start_x_norm)));
    arr_b_s_y = clamp01(data2fig_y(norm2data_y(cfg.backward.arrow_start_y_norm)));
    arr_b_e_x = clamp01(data2fig_x(norm2data_x(cfg.backward.arrow_end_x_norm)));
    arr_b_e_y = clamp01(data2fig_y(norm2data_y(cfg.backward.arrow_end_y_norm)));
    annotation(h, 'arrow', [arr_b_s_x, arr_b_e_x], [arr_b_s_y, arr_b_e_y], ...
        'Color', [0.85 0.325 0.098], 'LineWidth', 2, 'HeadStyle', 'vback2', 'HeadLength', 8);
    tb = text(txt_b_data_x, txt_b_data_y, '$\sigma_c^b$', ...
        'HorizontalAlignment', 'left', 'FontSize', 18, 'Interpreter', 'latex', 'Color', [0.85 0.325 0.098]);

    uistack([t1, t2, t3, tf, tb], 'top');

    % 单图导出 (以独立图片输出)
    test_out_img = fullfile(test_plots_dir, sprintf('phase_beta_sigma_%s.png', case_info.suffix));
    exportgraphics(h, test_out_img, 'Resolution', 300);
    fprintf('[DONE] 单个测试相图已保存 (%s): %s\n', case_info.label, test_out_img);
end

fprintf('\n所有相图处理与导出完毕！\n');
