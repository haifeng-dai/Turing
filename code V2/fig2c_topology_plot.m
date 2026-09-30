function image_files = fig2c_topology_plot(state, cases, output_dir, output_tag)
% 三张填色相图和一张无区域填色的汇总图。
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
ratio = state.config.ratios;
forward_color = [0 0.447 0.741];
backward_color = [0.85 0.325 0.098];
image_files = cell(1, 4);
for c = 1:3
    h = figure('Color', 'w', 'Units', 'pixels', 'Position', [100 100 900 800]);
    layout = tiledlayout(h, 1, 1, 'Padding', 'loose');
    ax = nexttile(layout);
    hold(ax, 'on');
    sf = state.sf_mean(c, :); sb = state.sb_mean(c, :);
    sf_std = state.sf_std(c, :); sb_std = state.sb_std(c, :);
    % 只对填色边界排序，原始正向/反向统计及散点不交换。
    lower = min(sf, sb); upper = max(sf, sb);
    [ymin, ymax] = plotLimits(sf, sb, sf_std, sb_std, state.sigma_spec(c, :));
    fill(ax, [ratio fliplr(ratio)], [lower ymin*ones(size(ratio))], ...
        [0.90 0.95 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
    fill(ax, [ratio fliplr(ratio)], [upper fliplr(lower)], ...
        [0.90 0.90 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
    fill(ax, [ratio fliplr(ratio)], [ymax*ones(size(ratio)) fliplr(upper)], ...
        [1.00 0.95 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
    addBand(ax, ratio, sf, sf_std, forward_color);
    addBand(ax, ratio, sb, sb_std, backward_color);
    hs = plot(ax, ratio, state.sigma_spec(c, :), '-', 'Color', forward_color, 'LineWidth', 2);
    hf = plot(ax, ratio, sf, 'o', 'LineStyle', 'none', 'Color', forward_color, 'MarkerSize', 4);
    hb = plot(ax, ratio, sb, 's', 'LineStyle', 'none', 'Color', backward_color, 'MarkerSize', 4);
    % 原相图区的三种区域标签，纵向位置由当前边界确定。
    idx = max(1, round(0.7*numel(ratio)));
    text(ax, ratio(idx), (ymin+lower(idx))/2, 'Uniform', 'HorizontalAlignment', 'center', 'FontSize', 16);
    if upper(idx)-lower(idx) > 0.02
        text(ax, ratio(idx), (lower(idx)+upper(idx))/2, 'Bistable', 'HorizontalAlignment', 'center', 'FontSize', 16);
    end
    text(ax, ratio(idx), (upper(idx)+ymax)/2, 'Pattern', 'HorizontalAlignment', 'center', 'FontSize', 16);
    decorateAxes(ax, ratio, ymin, ymax);
    title(ax, cases(c).label, 'Interpreter', 'none', 'FontSize', 18, 'FontName', 'Arial');
    subtitle(ax, sprintf('Layers: %s | alpha=%.2f, eta=%.2f | mean +/- 1 SD, n=%d', ...
        strjoin(cases(c).layer_types, '-'), state.config.alpha, state.config.noise, ...
        numel(state.config.noise_seeds)), 'Interpreter', 'none', 'FontSize', 10, 'FontName', 'Arial');
    legend(ax, [hs hf hb], {'Forward spectrum', 'Forward mean', ...
        'Backward mean'}, 'Location', 'northwest', 'FontSize', 10, 'FontName', 'Arial');
    image_files{c} = fullfile(output_dir, sprintf('fig2c_topology_%s_%s.png', cases(c).tag, output_tag));
    savePNG(h, image_files{c});
end

h = figure('Color', 'w', 'Units', 'pixels', 'Position', [100 100 1100 800]);
layout = tiledlayout(h, 1, 1, 'Padding', 'loose');
ax = nexttile(layout); hold(ax, 'on');
colors = lines(3);
handles = gobjects(1, 9);
labels = cell(1, 9);
for c = 1:3
    handles(3*c-2) = plot(ax, ratio, state.sigma_spec(c, :), '-', ...
        'Color', colors(c, :), 'LineWidth', 2);
    handles(3*c-1) = errorbar(ax, ratio, state.sf_mean(c, :), state.sf_std(c, :), ...
        'o', 'LineStyle', 'none', 'Color', colors(c, :), 'MarkerSize', 4, 'CapSize', 3);
    handles(3*c) = errorbar(ax, ratio, state.sb_mean(c, :), state.sb_std(c, :), ...
        's', 'LineStyle', 'none', 'Color', colors(c, :), 'MarkerSize', 4, 'CapSize', 3);
    labels(3*c-2:3*c) = {sprintf('%s: forward spectrum', cases(c).label), ...
        sprintf('%s: forward mean +/- SD', cases(c).label), ...
        sprintf('%s: backward mean +/- SD', cases(c).label)};
end
[ymin, ymax] = plotLimits(state.sf_mean, state.sb_mean, state.sf_std, ...
    state.sb_std, state.sigma_spec);
decorateAxes(ax, ratio, ymin, ymax);
title(ax, 'Topology comparison at fixed networks', 'FontSize', 18, 'FontName', 'Arial');
subtitle(ax, sprintf('N=%d, K=%d, alpha=%.2f, eta=%.2f, n=%d noise realizations', ...
    state.config.N, state.config.K, state.config.alpha, state.config.noise, ...
    numel(state.config.noise_seeds)), 'Interpreter', 'none', 'FontSize', 10, 'FontName', 'Arial');
legend(ax, handles, labels, 'Location', 'southoutside', 'NumColumns', 3, 'FontSize', 10, 'FontName', 'Arial');
image_files{4} = fullfile(output_dir, sprintf('fig2c_topology_comparison_%s.png', output_tag));
savePNG(h, image_files{4});
end

function savePNG(h, path)
% 先低分辨率渲染，刷新无桌面模式下的字体度量和图例自动布局。
drawnow;
set(h, 'PaperPositionMode', 'auto');
preview = print(h, '-RGBImage', '-r72'); %#ok<NASGU>
drawnow;
print(h, path, '-dpng', '-r300');
end

function addBand(ax, x, average, deviation, color)
fill(ax, [x fliplr(x)], [average+deviation fliplr(average-deviation)], ...
    color, 'EdgeColor', 'none', 'FaceAlpha', 0.18, 'HandleVisibility', 'off');
end

function decorateAxes(ax, ratio, ymin, ymax)
xlabel(ax, '$\beta/\alpha$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel(ax, '$\sigma$', 'Interpreter', 'latex', 'FontSize', 18);
xlim(ax, [ratio(1) ratio(end)]); ylim(ax, [ymin ymax]);
box(ax, 'on'); grid(ax, 'off');
set(ax, 'FontSize', 16, 'FontName', 'Arial', 'Layer', 'top', 'TickLabelInterpreter', 'latex');
end

function [ymin, ymax] = plotLimits(sf, sb, sf_std, sb_std, spec)
values = [sf(:)+sf_std(:); sb(:)+sb_std(:); spec(:)];
values = values(isfinite(values));
if isempty(values), error('没有可绘制的有限阈值。'); end
ymin = max(0, min(sb(:)-sb_std(:))*0.88);
ymax = max(values)*1.15;
if ymax <= ymin, ymax = ymin+1; end
end
