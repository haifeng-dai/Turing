function image_file = fig2c_noise_species_plot(state,image_file)
% 绘制单个 ER 拓扑下的已有双组分基准及 u-only/v-only 两条消融曲线。
ratio = state.config.ratios;
colors = [0.18 0.18 0.18; 0.00 0.447 0.741; 0.85 0.325 0.098];
condition_labels = {'Both noisy (existing Fig. 2c)', ...
    'u-only noise','v-only noise'};

figure_handle = figure('Color','w','Units','pixels','Position',[100 100 1050 760]);
ax = axes(figure_handle);
hold(ax,'on');
all_values = [state.baseline_sf(:);state.baseline_sb(:);state.sf(:);state.sb(:)];

for condition = 1:3
    if condition == 1
        sf = state.baseline_sf;
        sb = state.baseline_sb;
    else
        sf = state.sf(condition-1,:);
        sb = state.sb(condition-1,:);
    end
    plot(ax,ratio,sf,'-o','Color',colors(condition,:), ...
        'MarkerFaceColor',colors(condition,:),'LineWidth',2, ...
        'MarkerSize',4,'DisplayName',[condition_labels{condition} ': forward']);
    plot(ax,ratio,sb,'--s','Color',colors(condition,:), ...
        'MarkerFaceColor','w','LineWidth',2, ...
        'MarkerSize',4,'DisplayName',[condition_labels{condition} ': backward']);
end

finite_values = all_values(isfinite(all_values));
if isempty(finite_values)
    error('没有可绘制的有限阈值。');
end
ymin = min(finite_values);
ymax = max(finite_values);
padding = max(0.5,0.08*(ymax-ymin));
xlabel(ax,'$\beta/\alpha$','Interpreter','latex','FontSize',17);
ylabel(ax,'$\sigma$','Interpreter','latex','FontSize',17);
xlim(ax,[ratio(1) ratio(end)]);
ylim(ax,[max(0,ymin-padding) ymax+padding]);
title(ax,'Species-specific noise ablation at the original Fig. 2(c) ER setting', ...
    'Interpreter','none','FontName','Arial','FontSize',16);
subtitle(ax,sprintf(['N=%d, K=%d, p=%.3f, alpha=%.2f, eta=%.2f; ' ...
    'new conditions use fixed noise seed=%d'], ...
    state.config.N,state.config.K,state.config.p,state.config.alpha, ...
    state.config.noise,state.config.noise_seed), ...
    'Interpreter','none','FontName','Arial','FontSize',11);
legend(ax,'Location','best','NumColumns',2,'FontName','Arial','FontSize',10);
grid(ax,'on');
box(ax,'on');
set(ax,'FontName','Arial','FontSize',13,'Layer','top', ...
    'TickLabelInterpreter','latex');
exportgraphics(figure_handle,image_file,'Resolution',300);
close(figure_handle);
end
