%% Plot Fig 07

% Load data
load('../res/res_08','rmse*');


figure('Position',[0 0 900 600], 'PaperUnits', 'inches', 'PaperSize', 6.92*[1,2/3]);
% All subjects

for i = 1:6
    plt_fn(2, 3, i, rmse_H0(i,:,:), rmse_H2(i,:,:), rmse_H1H2(i,:,:), sprintf('Subject %d',i));
end
subplot(2,3,1);
legend({'VC','PRIMO','TxLR'},'Interpreter','latex','Location','NorthWest');

saveas(gcf,'fig_07');
saveas(gcf,'fig_07.svg');
saveas(gcf,'fig_07','epsc');

%% Helper function
function plt_fn(m, n, idx, a, b, c, lbl)
    subplot(m, n, idx);
    bar([a(1,:,10);b(1,:,10);c(1,:,10)]');
    grid on;
    xlabel('Acceleration Factor (R)','Interpreter','latex');
    xticklabels(2:8);
    ylabel('RMSE','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    ylim([0 0.75]);
    xlim([0 8]);
    yticks(0:0.1:1);
    title(lbl,'Interpreter','latex');
    set(gca,'FontSize',16);
    hold on;
    plot((1:7)-0.225,min(a,[],3),'LineStyle','None','Marker','.','MarkerSize',16,'color','k');
    plot((1:7)-0.000,min(b,[],3),'LineStyle','None','Marker','.','MarkerSize',16,'color','k');
    plot((1:7)+0.225,min(c,[],3),'LineStyle','None','Marker','.','MarkerSize',16,'color','k');
end
