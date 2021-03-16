%% Plot Fig S2
figure('Position',[0 0 800 800], 'PaperUnits', 'inches', 'PaperSize', 6.92*[1,1]);
load('../res/res_S2');

subplot(2,2,3);
pltfn(rmse_H1H2,'TxLR', '(c)');


subplot(2,2,2);
pltfn(rmse_H2,'PRIMO', '(b)');

subplot(2,2,1);
pltfn(rmse_H0,'VC', '(a)');

subplot(2,2,4);
cla();hold on;
plot(rmse_H0(:,50),'linewidth',2);
plot(rmse_H2(:,50),'linewidth',2);
plot(rmse_H1H2(:,50),'linewidth',2);
text(0,1.1,'(d)','Interpreter','latex', 'FontSize', 24);
axis square;grid on;
set(gca,'FontSize',16);
xlabel('Iterations @ Rank 50','Interpreter','latex');
ylabel('RMSE','Interpreter','latex');
legend({'VC','PRIMO','TxLR'},'Interpreter','latex');
xticks(20:20:100);
yticks(0:0.2:1);
set(gca,'TickLabelInterpreter','latex');
c=colorbar();
c.Visible = false;


saveas(gcf,'fig_S2');
saveas(gcf,'fig_S2.svg');
saveas(gcf,'fig_S2','epsc');

%% Helper function
function pltfn(data, label, id)
    imagesc(data);
    axis xy square;
    xticks(20:20:100);
    xlabel('Rank Threshold','Interpreter','latex');
    ylabel('Iterations','Interpreter','latex')
    yticks(20:20:100);
    caxis([0 1]);
    c=colorbar();
    ylabel(c,'RMSE');
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter='latex';
    set(gca,'FontSize',16);
    title(label,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    text(0,110,id,'Interpreter','latex', 'FontSize', 24);
end
