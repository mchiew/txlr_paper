% Plot Fig 3

load('../res/res_03','rmse*');

figure('Position',[0 0 400 400], 'PaperUnits', 'inches', 'PaperSize', 6.92*[0.5,0.5]);

clf();hold on;
ax = gca;
set(ax,'FontSize',16);

col = lines(3);
plt_fn(mean(rmse_H0,2), std(rmse_H0,[],2), col(1,:));
plt_fn(mean(rmse_H2,2), std(rmse_H2,[],2), col(2,:));
plt_fn(mean(rmse_H1H2,2), std(rmse_H1H2,[],2), col(3,:));

xlim([0.5 12.5]);
ylim([0 1.1]);
xticks(1:12);
yticks(0.1:0.1:1.1);
grid on;
axis square;

set(gca,'TickLabelInterpreter','latex');
ylabel('RMSE','Interpreter','latex');
xlabel('R','Interpreter','latex');
legend({'VC','PRIMO','TxLR'},'Interpreter','latex','Location','NorthWest');

saveas(gcf,'fig_03');
saveas(gcf,'fig_03.svg');
saveas(gcf,'fig_03','epsc');
%% Helper functions

function plt_fn(m, s, col)
    h = patch([1:length(m) length(m):-1:1], [(m-s)' fliplr((m+s)')],1,'facecolor',col,'edgecolor','none','facealpha',0.5);
    h.HandleVisibility = 'off';
    plot(mean(m,2), 'Color',col, 'LineWidth', 2);
end