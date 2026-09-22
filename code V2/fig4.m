%% fig4.m
% Combine fig4a, b, c into a 1x3 tiledlayout.
% Panel (a): ER network hysteresis vs N
% Panel (b): WS network hysteresis vs N
% Panel (c): BA network hysteresis vs N

clear; clc; close all;

%% 1. Config
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'V2', 'manuscirpt', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

% Common base parameters
DYNA_BASE.K = 5;
DYNA_BASE.alpha = 0.05;
DYNA_BASE.beta  = 0.1 * DYNA_BASE.alpha;
DYNA_BASE.noise = 0.1;
N_VALS = [100, 500, 1000];

%% 2. Plotting
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.3, 0.90, 0.33]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel (a): ER Network ---
ax1 = nexttile(t); hold on;
TOPO_PARAM_ER = 0.03;
colors = lines(length(N_VALS));
for i = 1:length(N_VALS)
    mat_name = sprintf('hysteresis_er_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        N_VALS(i), DYNA_BASE.K, TOPO_PARAM_ER, DYNA_BASE.alpha, DYNA_BASE.beta, DYNA_BASE.noise);
    data_path = fullfile(res_dir, mat_name);
    if exist(data_path, 'file')
        res = load(data_path);
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $N=%d$', N_VALS(i)), ...
            'MarkerIndices', 1:5:length(res.sigma_range));
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $N=%d$', N_VALS(i)), ...
            'MarkerIndices', 1:5:length(res.sigma_range));
    end
end
xlabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$A(\sigma)$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 25]); ylim([0, 140]);
legend('Location', 'NorthEast', 'FontSize', 13, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'NumColumns', 2);
grid on; box on;
set(ax1, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
text(0.05, 0.5, '\textbf{(a)}', 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman');

% --- Panel (b): WS Network ---
ax2 = nexttile(t); hold on;
K_AVG_WS = 6; P_REWIRE_WS = 0.1;
colors = lines(length(N_VALS));
for i = 1:length(N_VALS)
    mat_name = sprintf('hysteresis_ws_N%d_K%d_pr%.2f_mK%d_a%.3f_b%.3f_n%.2f_results.mat', ...
        N_VALS(i), DYNA_BASE.K, P_REWIRE_WS, K_AVG_WS, DYNA_BASE.alpha, DYNA_BASE.beta, DYNA_BASE.noise);
    data_path = fullfile(res_dir, mat_name);
    if exist(data_path, 'file')
        res = load(data_path);
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $N=%d$', N_VALS(i)), ...
            'MarkerIndices', 1:5:length(res.sigma_range));
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $N=%d$', N_VALS(i)), ...
            'MarkerIndices', 1:5:length(res.sigma_range));
    end
end
xlabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 35]); ylim([0, 140]);
legend('Location', 'NorthEast', 'FontSize', 13, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'NumColumns', 2);
grid on; box on;
set(ax2, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
text(0.05, 0.5, '\textbf{(b)}', 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman');

% --- Panel (c): BA Network ---
ax3 = nexttile(t); hold on;
M_BA = 3;
colors = lines(length(N_VALS));
for i = 1:length(N_VALS)
    mat_name = sprintf('hysteresis_ba_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        N_VALS(i), DYNA_BASE.K, M_BA, DYNA_BASE.alpha, DYNA_BASE.beta, DYNA_BASE.noise);
    data_path = fullfile(res_dir, mat_name);
    if exist(data_path, 'file')
        res = load(data_path);
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $N=%d$', N_VALS(i)), ...
            'MarkerIndices', 1:5:length(res.sigma_range));
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $N=%d$', N_VALS(i)), ...
            'MarkerIndices', 1:5:length(res.sigma_range));
    end
end
xlabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 25]); ylim([0, 100]);
legend('Location', 'NorthEast', 'FontSize', 13, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'NumColumns', 2);
grid on; box on;
set(ax3, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
text(0.05, 0.5, '\textbf{(c)}', 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman');

% Export
out_path = fullfile(plots_dir, 'fig4.eps');
exportgraphics(h, out_path, 'ContentType', 'vector');
fprintf('[DONE] Combined Figure 4 saved to: %s\n', out_path);
