%% fig5.m
% Combine fig5a, b, c into a 1x3 tiledlayout.
% Panel (a): ER network hysteresis vs K
% Panel (b): WS network hysteresis vs K
% Panel (c): BA network hysteresis vs K

clear; clc; close all;

%% 1. Config
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'V2', 'manuscirpt', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

% Common base parameters
DYNA_BASE.N = 200;
DYNA_BASE.alpha = 0.05;
DYNA_BASE.beta  = 0.1 * DYNA_BASE.alpha;
DYNA_BASE.noise = 0.1;

K_VALS = 2:4:10; % [2, 6, 10]

%% 2. Plotting
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.3, 0.90, 0.33]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel (a): ER Network ---
ax1 = nexttile(t); hold on;
TOPO_PARAM_ER = 0.03;
colors = lines(length(K_VALS));
for i = 1:length(K_VALS)
    mat_name = sprintf('hysteresis_er_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        DYNA_BASE.N, K_VALS(i), TOPO_PARAM_ER, DYNA_BASE.alpha, DYNA_BASE.beta, DYNA_BASE.noise);
    data_path = fullfile(res_dir, mat_name);
    if exist(data_path, 'file')
        res = load(data_path);
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $K=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:10:length(res.sigma_range));
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $K=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:10:length(res.sigma_range));
    end
end
xlabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$A(\sigma)$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 30]); ylim([0, 140]);
legend('Location', 'NorthEast', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'NumColumns', 2);
grid on; box on;
set(ax1, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
text(0.05, 0.5, '\textbf{(a)}', 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman');

% --- Panel (b): WS Network ---
ax2 = nexttile(t); hold on;
K_AVG_WS = 6; P_REWIRE_WS = 0.1;
colors = lines(length(K_VALS));
for i = 1:length(K_VALS)
    mat_name = sprintf('hysteresis_ws_N%d_K%d_pr%.2f_mK%d_a%.3f_b%.3f_n%.2f_results.mat', ...
        DYNA_BASE.N, K_VALS(i), P_REWIRE_WS, K_AVG_WS, DYNA_BASE.alpha, DYNA_BASE.beta, DYNA_BASE.noise);
    data_path = fullfile(res_dir, mat_name);
    if exist(data_path, 'file')
        res = load(data_path);
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $K=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:5:length(res.sigma_range));
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $K=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:5:length(res.sigma_range));
    end
end
xlabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 40]); ylim([0, 140]);
legend('Location', 'NorthEast', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'NumColumns', 2);
grid on; box on;
set(ax2, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
text(0.05, 0.5, '\textbf{(b)}', 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman');

% --- Panel (c): BA Network ---
ax3 = nexttile(t); hold on;
M_BA = 3;
colors = lines(length(K_VALS));
for i = 1:length(K_VALS)
    mat_name = sprintf('hysteresis_%s_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        'ba', DYNA_BASE.N, K_VALS(i), M_BA, DYNA_BASE.alpha, DYNA_BASE.beta, DYNA_BASE.noise);
    data_path = fullfile(res_dir, mat_name);
    if exist(data_path, 'file')
        res = load(data_path);
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $K=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:5:length(res.sigma_range));
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $K=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:5:length(res.sigma_range));
    end
end
xlabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 25]); ylim([0, 100]);
legend('Location', 'NorthEast', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'NumColumns', 2);
grid on; box on;
set(ax3, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
text(0.05, 0.5, '\textbf{(c)}', 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman');

% Export
% out_path = fullfile(plots_dir, 'fig5.eps');
% exportgraphics(h, out_path, 'ContentType', 'vector');
test_plots_dir = fullfile(fileparts(mfilename('fullpath')), 'fig');
if ~exist(test_plots_dir, 'dir'), mkdir(test_plots_dir); end
test_out_path = fullfile(test_plots_dir, 'fig5.png');
exportgraphics(h, test_out_path, 'Resolution', 300);
fprintf('[DONE] Test Figure 5 saved to: %s\n', test_out_path);
