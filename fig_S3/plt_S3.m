% Plot Fig S3

load('../res/res_S3','rmse*');

figure('Position',[0 0 400 400], 'PaperUnits', 'inches', 'PaperSize', 6.92*[0.5,0.5]);
clf();hold on;
ax = gca;
set(ax,'FontSize',16);

for i = 1:4
    ax.ColorOrderIndex = 1;
    scatter(sum(abs(rmse_18(:,i)-rmse_18(:,i)')<0.005,2).*randn(48,1)/250 + (i-1)*4+1 + (i-1),rmse_18(:,i),50,'filled')
    scatter(sum(abs(rmse_24(:,i)-rmse_24(:,i)')<0.005,2).*randn(48,1)/250 + (i-1)*4+2 + (i-1),rmse_24(:,i),50,'filled')
    scatter(sum(abs(rmse_36(:,i)-rmse_36(:,i)')<0.005,2).*randn(48,1)/250 + (i-1)*4+3 + (i-1),rmse_36(:,i),50,'filled')
    scatter(sum(abs(rmse_48(:,i)-rmse_48(:,i)')<0.005,2).*randn(48,1)/250 + (i-1)*4+4 + (i-1),rmse_48(:,i),50,'filled')
end
alpha(0.5);

boxplot([rmse_18(:,1) rmse_24(:,1) rmse_36(:,1) rmse_48(:,1),...
         rmse_18(:,2) rmse_24(:,2) rmse_36(:,2) rmse_48(:,2),...
         rmse_18(:,3) rmse_24(:,3) rmse_36(:,3) rmse_48(:,3),...
         rmse_18(:,4) rmse_24(:,4) rmse_36(:,4) rmse_48(:,4)],'positions',[1:4 6:9 11:14 16:19]);
set(findobj(gca,'Type','Line'),'Linewidth',1);
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2 0.75])
set(findobj(gca,'Tag','Outliers'),'Marker','None');

axis square;
grid on;
xlim([0.5 19.5]);
ylim([0 0.25]);
yticks(0:0.025:0.25);
xticks([1:4 6:9 11:14 16:19])

set(gca,'TickLabelInterpreter','latex');
xticklabels({'',' R=2','','','',' R=4','','','',' R=6','','','',' R=8','','',''})
ylabel('RMSE','Interpreter','latex');
legend({'matrix size = 18$\times$18','matrix size = 24$\times$24','matrix size = 36$\times$36','matrix size = 48$\times$48'},'Interpreter','latex', 'NumColumns',2, 'FontSize',10)

saveas(gcf,'fig_S3');
saveas(gcf,'fig_S3.svg');
saveas(gcf,'fig_S3','epsc');
