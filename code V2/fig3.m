%% fig3.m
% Combine fig3a, b, c into a 1x3 tiledlayout.
% Panel (a): ER network hysteresis (p)
% Panel (b): WS network hysteresis (k)
% Panel (c): BA network hysteresis (m)

clear; clc; close all;

%% 1. Config
res_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
plots_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'manuscript', 'V2', 'manuscirpt', 'figures');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

% DYNA settings (common across a, b, c)
DYNA.N = 200;
DYNA.K = 5;
DYNA.alpha = 0.05;
DYNA.beta  = 0.1 * DYNA.alpha;
DYNA.noise = 0.1;

%% 2. Plotting
h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.3, 0.90, 0.33]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel (a): ER Network ---
ax1 = nexttile(t); hold on;
P_VALS = [0.03, 0.04, 0.10];
colors = lines(length(P_VALS));
for i = 1:length(P_VALS)
    mat_name = sprintf('hysteresis_er_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        DYNA.N, DYNA.K, P_VALS(i), DYNA.alpha, DYNA.beta, DYNA.noise);
    data_path = fullfile(res_dir, mat_name);
    if exist(data_path, 'file')
        res = load(data_path);
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $p=%.2f$', P_VALS(i)), ...
            'MarkerIndices', 1:4:length(res.sigma_range));
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $p=%.2f$', P_VALS(i)), ...
            'MarkerIndices', 1:4:length(res.sigma_range));
    end
end
xlabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$A(\sigma)$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 25]); ylim([0, 140]);
legend('Location', 'NorthEast', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'NumColumns', 2);
grid on; box on;
set(ax1, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
text(0.05, 0.5, '\textbf{(a)}', 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman');

% --- Panel (b): WS Network ---
ax2 = nexttile(t); hold on;
K_VALS = [6, 8, 20];
P_REWIRE = 0.01;
colors = lines(length(K_VALS));
for i = 1:length(K_VALS)
    mat_name = sprintf('hysteresis_ws_N%d_K%d_pr%.2f_mK%d_a%.3f_b%.3f_n%.2f_results.mat', ...
        DYNA.N, DYNA.K, P_REWIRE, K_VALS(i), DYNA.alpha, DYNA.beta, DYNA.noise);
    data_path = fullfile(res_dir, mat_name);
    if exist(data_path, 'file')
        res = load(data_path);
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $k=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:4:length(res.sigma_range));
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $k=%d$', K_VALS(i)), ...
            'MarkerIndices', 1:4:length(res.sigma_range));
    end
end
xlabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 30]); ylim([0, 140]);
legend('Location', 'NorthEast', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'NumColumns', 2);
grid on; box on;
set(ax2, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
text(0.05, 0.5, '\textbf{(b)}', 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman');

% --- Panel (c): BA Network ---
ax3 = nexttile(t); hold on;
M_VALS = [3, 4, 10];
colors = lines(length(M_VALS));
for i = 1:length(M_VALS)
    mat_name = sprintf('hysteresis_ba_N%d_K%d_p%.3f_a%.3f_b%.3f_n%.2f_results.mat', ...
        DYNA.N, DYNA.K, M_VALS(i), DYNA.alpha, DYNA.beta, DYNA.noise);
    data_path = fullfile(res_dir, mat_name);
    if exist(data_path, 'file')
        res = load(data_path);
        plot(res.sigma_range, res.A_fwd, '-o', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Fwd, $m=%d$', M_VALS(i)), ...
            'MarkerIndices', 1:4:length(res.sigma_range));
        plot(res.sigma_range, res.A_bwd, '--^', 'Color', colors(i,:), 'LineWidth', 2, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Bwd, $m=%d$', M_VALS(i)), ...
            'MarkerIndices', 1:4:length(res.sigma_range));
    end
end
xlabel('$\sigma$', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
xlim([10, 18]); ylim([0, 80]);
legend('Location', 'NorthEast', 'FontSize', 14, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'NumColumns', 2);
grid on; box on;
set(ax3, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
text(0.05, 0.5, '\textbf{(c)}', 'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman');

% Export
out_path = fullfile(plots_dir, 'fig3.eps');
exportgraphics(h, out_path, 'ContentType', 'vector');
fprintf('[DONE] Combined Figure 3 saved to: %s\n', out_path);
