%% fig2.m
% Combine fig2abc and fig2def into a single 2x3 layout.
% Row 1: (beta/alpha, sigma) phase diagrams for selected alpha values.
% Row 2: (alpha, sigma) phase diagrams for selected beta/alpha ratios.
clear; clc; close all;

%% 1. Config
TOPO_TYPE = 'ER';
N = 200;
K = 5;
P_VAL_FIXED = 0.030;
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

% --- Row 1 Config (alpha values) ---
SELECTED_ALPHAS = [0.03, 0.05, 0.10];
ALPHA_TOL = 1e-4;
data_path_abc = fullfile(res_dir, sprintf('scan_beta_alpha_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), N, K, P_VAL_FIXED));

% alpha_config for (a, b, c)
alpha_config(1).alpha = 0.03;
alpha_config(1).forward.txt_x_norm = 0.25; alpha_config(1).forward.txt_y_norm = 0.65;
alpha_config(1).forward.arrow_start_x_norm = 0.26; alpha_config(1).forward.arrow_start_y_norm = 0.6;
alpha_config(1).forward.arrow_end_x_norm = 0.23; alpha_config(1).forward.arrow_end_y_norm = 0.5;
alpha_config(1).backward.txt_x_norm = 0.06; alpha_config(1).backward.txt_y_norm = 0.08;
alpha_config(1).backward.arrow_start_x_norm = 0.1; alpha_config(1).backward.arrow_start_y_norm = 0.13;
alpha_config(1).backward.arrow_end_x_norm = 0.15; alpha_config(1).backward.arrow_end_y_norm = 0.23;
alpha_config(1).region_x_norms = [0.7, 0.35, 0.75]; alpha_config(1).region_y_norms = [0.10, 0.3, 0.55];

alpha_config(2).alpha = 0.05;
alpha_config(2).forward.txt_x_norm = 0.20; alpha_config(2).forward.txt_y_norm = 0.62;
alpha_config(2).forward.arrow_start_x_norm = 0.2; alpha_config(2).forward.arrow_start_y_norm = 0.6;
alpha_config(2).forward.arrow_end_x_norm = 0.1; alpha_config(2).forward.arrow_end_y_norm = 0.48;
alpha_config(2).backward.txt_x_norm = 0.05; alpha_config(2).backward.txt_y_norm = 0.33;
alpha_config(2).backward.arrow_start_x_norm = 0.08; alpha_config(2).backward.arrow_start_y_norm = 0.27;
alpha_config(2).backward.arrow_end_x_norm = 0.15; alpha_config(2).backward.arrow_end_y_norm = 0.14;
alpha_config(2).region_x_norms = [0.7, 0.35, 0.6]; alpha_config(2).region_y_norms = [0.10, 0.2, 0.70];

alpha_config(3).alpha = 0.10;
alpha_config(3).forward.txt_x_norm = 0.05; alpha_config(3).forward.txt_y_norm = 0.65;
alpha_config(3).forward.arrow_start_x_norm = 0.08; alpha_config(3).forward.arrow_start_y_norm = 0.6;
alpha_config(3).forward.arrow_end_x_norm = 0.12; alpha_config(3).forward.arrow_end_y_norm = 0.45;
alpha_config(3).backward.txt_x_norm = 0.08; alpha_config(3).backward.txt_y_norm = 0.05;
alpha_config(3).backward.arrow_start_x_norm = 0.13; alpha_config(3).backward.arrow_start_y_norm = 0.08;
alpha_config(3).backward.arrow_end_x_norm = 0.2; alpha_config(3).backward.arrow_end_y_norm = 0.2;
alpha_config(3).region_x_norms = [0.55, 0.7, 0.55]; alpha_config(3).region_y_norms = [0.1, 0.5, 0.92];

% --- Row 2 Config (ratio values) ---
SELECTED_RATIOS = [0.1, 6.0, 10.0];
RATIO_TOL = 1e-4;
data_path_def = fullfile(res_dir, sprintf('scan_alpha_ratio_%s_N%d_K%d_p%.3f.mat', ...
    lower(TOPO_TYPE), N, K, P_VAL_FIXED));

% ratio_config for (d, e, f)
ratio_config(1).ratio = 0.1;
ratio_config(1).forward.txt_x_norm = 0.15; ratio_config(1).forward.txt_y_norm = 0.87;
ratio_config(1).forward.arrow_start_x_norm = 0.15; ratio_config(1).forward.arrow_start_y_norm = 0.85;
ratio_config(1).forward.arrow_end_x_norm = 0.05; ratio_config(1).forward.arrow_end_y_norm = 0.74;
ratio_config(1).backward.txt_x_norm = 0.08; ratio_config(1).backward.txt_y_norm = 0.10;
ratio_config(1).backward.arrow_start_x_norm = 0.11; ratio_config(1).backward.arrow_start_y_norm = 0.13;
ratio_config(1).backward.arrow_end_x_norm = 0.21; ratio_config(1).backward.arrow_end_y_norm = 0.23;
ratio_config(1).region_x_norms = [0.7, 0.4, 0.6]; ratio_config(1).region_y_norms = [0.08, 0.2, 0.72];

ratio_config(2).ratio = 1.0;
ratio_config(2).forward.txt_x_norm = 0.15; ratio_config(2).forward.txt_y_norm = 0.87;
ratio_config(2).forward.arrow_start_x_norm = 0.15; ratio_config(2).forward.arrow_start_y_norm = 0.85;
ratio_config(2).forward.arrow_end_x_norm = 0.05; ratio_config(2).forward.arrow_end_y_norm = 0.74;
ratio_config(2).backward.txt_x_norm = 0.02; ratio_config(2).backward.txt_y_norm = 0.10;
ratio_config(2).backward.arrow_start_x_norm = 0.07; ratio_config(2).backward.arrow_start_y_norm = 0.15;
ratio_config(2).backward.arrow_end_x_norm = 0.17; ratio_config(2).backward.arrow_end_y_norm = 0.27;
ratio_config(2).region_x_norms = [0.5, 0.83, 0.35]; ratio_config(2).region_y_norms = [0.2, 0.47, 0.75];

ratio_config(3).ratio = 10.0;
ratio_config(3).forward.txt_x_norm = 0.17; ratio_config(3).forward.txt_y_norm = 0.49;
ratio_config(3).forward.arrow_start_x_norm = 0.2; ratio_config(3).forward.arrow_start_y_norm = 0.43;
ratio_config(3).forward.arrow_end_x_norm = 0.20; ratio_config(3).forward.arrow_end_y_norm = 0.3;
ratio_config(3).backward.txt_x_norm = 0.35; ratio_config(3).backward.txt_y_norm = 0.1;
ratio_config(3).backward.arrow_start_x_norm = 0.35; ratio_config(3).backward.arrow_start_y_norm = 0.12;
ratio_config(3).backward.arrow_end_x_norm = 0.25; ratio_config(3).backward.arrow_end_y_norm = 0.2;
ratio_config(3).region_x_norms = [0.8, 0.65, 0.8]; ratio_config(3).region_y_norms = [0.08, 0.5, 0.9];

%% 2. Load Data
if ~exist(data_path_abc, 'file') || ~exist(data_path_def, 'file')
    error('Missing data files. Ensure sigma_beta_alpha.m and sigma_alpha_beta.m have been run.');
end
data_abc = load(data_path_abc);
data_def = load(data_path_def);

%% 3. Plotting
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.75, 0.55]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Helpers
clamp01 = @(v) max(0, min(1, v));

% --- Combined 2x3 Layout ---
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
axs = gobjects(6, 1);

% Row 1
for i = 1:3
    axs(i) = nexttile(t); hold on;
    target_alpha = SELECTED_ALPHAS(i);
    [~, idx_a] = min(abs(data_abc.ALPHA_LIST - target_alpha));
    actual_alpha = data_abc.ALPHA_LIST(idx_a);
    x_data = data_abc.RATIO_LIST(:)';
    sf_vec = data_abc.SF_Matrix(idx_a, :);
    sb_vec = data_abc.SB_Matrix(idx_a, :);
    mask = sf_vec < sb_vec;
    if any(mask), [sf_vec(mask), sb_vec(mask)] = deal(sb_vec(mask), sf_vec(mask)); end
    max_val = max(sf_vec) * 1.12; min_val = max(0, min(sb_vec) * 0.88);
    if i == 1, min_val = 8; end
    fill_max_vec = ones(size(x_data)) * max_val; fill_min_vec = ones(size(x_data)) * min_val;
    fill([x_data, fliplr(x_data)], [sb_vec, fill_min_vec], [0.90 0.95 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
    fill([x_data, fliplr(x_data)], [sf_vec, fliplr(sb_vec)], [0.90 0.90 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
    fill([x_data, fliplr(x_data)], [fill_max_vec, fliplr(sf_vec)], [1.00 0.95 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
    plot(x_data, sf_vec, '--', 'Color', [0.00 0.447 0.741], 'LineWidth', 2, 'DisplayName', 'Forward $\sigma_c^f$');
    plot(x_data, sb_vec, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 2, 'DisplayName', 'Backward $\sigma_c^b$');
    xlabel(sprintf('$\\beta/\\alpha$ ($\\alpha = %.2f$)', actual_alpha), 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    if i == 1, ylabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman'); end
    xlim([0, 10]); ylim([min_val, max_val]);
    set(axs(i), 'FontSize', 14, 'Layer', 'top', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
    % title(sprintf('$\\alpha = %.2f$', actual_alpha), 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
    if i == 1, legend('Location', 'NorthEast', 'FontSize', 12, 'Interpreter', 'latex', 'FontName', 'Times New Roman'); end
    grid off; box on;
    text(0.05, 0.95, sprintf('\\textbf{(%s)}', char('a'+i-1)), 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'FontName', 'Times New Roman');

    % Map coordinates for annotations
    drawnow;
    ax_pos = get(axs(i), 'Position');
    xl = xlim(axs(i)); yl = ylim(axs(i));
    norm2data_x = @(xn) xl(1) + xn * (xl(2) - xl(1));
    norm2data_y = @(yn) yl(1) + yn * (yl(2) - yl(1));
    data2fig_x = @(xd) ax_pos(1) + (xd - xl(1)) / (xl(2) - xl(1)) * ax_pos(3);
    data2fig_y = @(yd) ax_pos(2) + (yd - yl(1)) / (yl(2) - yl(1)) * ax_pos(4);

    config = alpha_config(i);
    % Region text
    text(norm2data_x(config.region_x_norms(1)), min_val + config.region_y_norms(1)*(max_val-min_val), 'Stable', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    text(norm2data_x(config.region_x_norms(2)), min_val + config.region_y_norms(2)*(max_val-min_val), 'Hysteresis', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    text(norm2data_x(config.region_x_norms(3)), min_val + config.region_y_norms(3)*(max_val-min_val), 'Turing', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    % Forward Arrow + Label
    af_sx = clamp01(data2fig_x(norm2data_x(config.forward.arrow_start_x_norm))); af_sy = clamp01(data2fig_y(norm2data_y(config.forward.arrow_start_y_norm)));
    af_ex = clamp01(data2fig_x(norm2data_x(config.forward.arrow_end_x_norm))); af_ey = clamp01(data2fig_y(norm2data_y(config.forward.arrow_end_y_norm)));
    % annotation(h, 'arrow', [af_sx, af_ex], [af_sy, af_ey], 'Color', [0.00 0.447 0.741], 'LineWidth', 1.5, 'HeadStyle', 'vback2', 'HeadLength', 6);
    % text(norm2data_x(config.forward.txt_x_norm), norm2data_y(config.forward.txt_y_norm), '$\sigma_c^f$', 'Color', [0.00 0.447 0.741], 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    % Backward Arrow + Label
    ab_sx = clamp01(data2fig_x(norm2data_x(config.backward.arrow_start_x_norm))); ab_sy = clamp01(data2fig_y(norm2data_y(config.backward.arrow_start_y_norm)));
    ab_ex = clamp01(data2fig_x(norm2data_x(config.backward.arrow_end_x_norm))); ab_ey = clamp01(data2fig_y(norm2data_y(config.backward.arrow_end_y_norm)));
    % annotation(h, 'arrow', [ab_sx, ab_ex], [ab_sy, ab_ey], 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5, 'HeadStyle', 'vback2', 'HeadLength', 6);
    % text(norm2data_x(config.backward.txt_x_norm), norm2data_y(config.backward.txt_y_norm), '$\sigma_c^b$', 'Color', [0.85 0.325 0.098], 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
end

% Row 2
for i = 1:3
    axs(i+3) = nexttile(t); hold on;
    target_ratio = SELECTED_RATIOS(i);
    [~, idx_r] = min(abs(data_def.RATIO_LIST - target_ratio));
    actual_ratio = data_def.RATIO_LIST(idx_r);
    x_data = data_def.ALPHA_LIST(:)';
    sf_vec = data_def.SF_Matrix(:, idx_r)'; sb_vec = data_def.SB_Matrix(:, idx_r)';
    max_val = max(sf_vec) * 1.15; min_val = 0;
    if i == 3, min_val = 10; end
    fill_max_vec = ones(size(x_data)) * max_val; fill_min_vec = ones(size(x_data)) * min_val;
    % Fills
    fill([x_data, fliplr(x_data)], [sb_vec, fill_min_vec], [0.90 0.95 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
    fill([x_data, fliplr(x_data)], [sf_vec, fliplr(sb_vec)], [0.90 0.90 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
    fill([x_data, fliplr(x_data)], [fill_max_vec, fliplr(sf_vec)], [1.00 0.95 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
    plot(x_data, sf_vec, '--', 'Color', [0.00 0.447 0.741], 'LineWidth', 2, 'DisplayName', 'Forward $\sigma_c^f$');
    plot(x_data, sb_vec, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 2, 'DisplayName', 'Backward $\sigma_c^b$');
    xlabel(sprintf('$\\alpha$ ($\\beta/\\alpha = %.1f$)', actual_ratio), 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    if i == 1, ylabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman'); end
    xlim([0.01, 0.1]); ylim([min_val, max_val]);
    set(axs(i+3), 'FontSize', 14, 'Layer', 'top', 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
    % title(sprintf('$\\beta/\\alpha = %.1f$', actual_ratio), 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
    grid off; box on;
    text(0.05, 0.95, sprintf('\\textbf{(%s)}', char('d'+i-1)), 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'FontName', 'Times New Roman');

    % Map coordinates for annotations
    drawnow;
    ax_pos = get(axs(i+3), 'Position');
    xl = xlim(axs(i+3)); yl = ylim(axs(i+3));
    norm2data_x = @(xn) xl(1) + xn * (xl(2) - xl(1));
    norm2data_y = @(yn) yl(1) + yn * (yl(2) - yl(1));
    data2fig_x = @(xd) ax_pos(1) + (xd - xl(1)) / (xl(2) - xl(1)) * ax_pos(3);
    data2fig_y = @(yd) ax_pos(2) + (yd - yl(1)) / (yl(2) - yl(1)) * ax_pos(4);

    config = ratio_config(i);
    % Region text
    text(norm2data_x(config.region_x_norms(1)), min_val + config.region_y_norms(1)*(max_val-min_val), 'Stable', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    text(norm2data_x(config.region_x_norms(2)), min_val + config.region_y_norms(2)*(max_val-min_val), 'Hysteresis', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    text(norm2data_x(config.region_x_norms(3)), min_val + config.region_y_norms(3)*(max_val-min_val), 'Turing', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    % Forward Arrow + Label
    af_sx = clamp01(data2fig_x(norm2data_x(config.forward.arrow_start_x_norm))); af_sy = clamp01(data2fig_y(norm2data_y(config.forward.arrow_start_y_norm)));
    af_ex = clamp01(data2fig_x(norm2data_x(config.forward.arrow_end_x_norm))); af_ey = clamp01(data2fig_y(norm2data_y(config.forward.arrow_end_y_norm)));
    % annotation(h, 'arrow', [af_sx, af_ex], [af_sy, af_ey], 'Color', [0.00 0.447 0.741], 'LineWidth', 1.5, 'HeadStyle', 'vback2', 'HeadLength', 6);
    % text(norm2data_x(config.forward.txt_x_norm), norm2data_y(config.forward.txt_y_norm), '$\sigma_c^f$', 'Color', [0.00 0.447 0.741], 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    % Backward Arrow + Label
    ab_sx = clamp01(data2fig_x(norm2data_x(config.backward.arrow_start_x_norm))); ab_sy = clamp01(data2fig_y(norm2data_y(config.backward.arrow_start_y_norm)));
    ab_ex = clamp01(data2fig_x(norm2data_x(config.backward.arrow_end_x_norm))); ab_ey = clamp01(data2fig_y(norm2data_y(config.backward.arrow_end_y_norm)));
    % annotation(h, 'arrow', [ab_sx, ab_ex], [ab_sy, ab_ey], 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5, 'HeadStyle', 'vback2', 'HeadLength', 6);
    % text(norm2data_x(config.backward.txt_x_norm), norm2data_y(config.backward.txt_y_norm), '$\sigma_c^b$', 'Color', [0.85 0.325 0.098], 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
end

% Export
out_path = fullfile(plots_dir, 'fig2.eps');
exportgraphics(h, out_path, 'ContentType', 'vector');
fprintf('[DONE] Combined Figure 2 saved to: %s\n', out_path);
