%% Plot Fig S6
figure('Position',[0 0 800 800], 'PaperUnits', 'inches', 'PaperSize', 6.92*[1,1]);
load('../res/res_S1');
load('../res/res_S6');

cla();hold on;
plot(rmse_H1H2(1:50,50),'linewidth',2);
plot(rmse_txlr,'linewidth',2);
axis square;grid on;
set(gca,'FontSize',16);
xlabel('Iterations @ Rank 50','Interpreter','latex');
ylabel('RMSE','Interpreter','latex');
legend({'TxLR, R=8, uniform poisson disc sampling','TxLR, R=8, poisson disc with fully sampled 4x4 centre'},'Interpreter','latex');
xticks(10:10:50);
yticks(0:0.2:1);
set(gca,'TickLabelInterpreter','latex');
c=colorbar();
c.Visible = false;


saveas(gcf,'fig_S6');
saveas(gcf,'fig_S6.svg');
saveas(gcf,'fig_S6','epsc');