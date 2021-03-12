%% Plot Fig S6

figure('Position',[0 0 800 800], 'PaperUnits', 'inches', 'PaperSize', 6.92*[1,1]);

load('../res/res_S6');
ids = 'abcd';
for i = 1:4
    subplot(4,1,i)
    plot(rmse(:,2:end,i));
    grid on;hold on
    scatter(N_Chi2(2:end,i),diag(rmse(N_Chi2(2:end,i),2:end,i)),150,lines(3),'filled')
    ylabel('RMSE','Interpreter','latex');
    text(138,0.025,sprintf('R=%d',2*i),'FontSize',20, 'Interpreter', 'latex');
    text(2,0.29,ids(i),'FontSize',32, 'Interpreter', 'latex');
    line([50 50],[0 1],'linewidth',2,'color','k','linestyle','--');
    yticks(0:0.05:0.25);
    xticks(25:25:150);
end

linkaxes(findobj(gcf,'Type','axes'));
xlim([0 150]);ylim([0 0.26]);
xlabel('Iterations','Interpreter','latex');
subplot(4,1,1);legend({'PSNR = 60 dB','PSNR = 50 dB','PSNR = 40 dB'},'Interpreter','latex','location','NorthEast','Orientation','horizontal');

set(findobj(gcf,'Type','axes'),'FontSize',16);
set(findobj(gcf,'Type','Line'),'LineWidth',2);
set(findobj(gcf,'Type','axes'),'TickLabelInterpreter','latex');

saveas(gcf,'fig_S6');
saveas(gcf,'fig_S6.svg');
saveas(gcf,'fig_S6','epsc');
