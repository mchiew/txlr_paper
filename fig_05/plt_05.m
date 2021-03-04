%% Plot Fig 5

figure('Position',[0 0 825 625], 'PaperUnits', 'inches', 'PaperSize', 6.92*[1,625/825]);

load('../res/res_05a');
plt_fn1(1, rmse, 'Rank Thresholds', 'Kernel Size', 'Fixed Iterations (50)',...
        {'10','20','30','40','50','60','70','80'},...
        {'[3,3]','[4,4]','[5,5]','[6,6]','[7,7]','[8,8]','[9,9]','[10,10]'});
line([3.5 8.5;3.5 8.5;3.5 3.5;8.5 8.5],[2.5 2.5;5.5 5.5;2.5 5.5;2.5 5.5],'linewidth',2,'color','r');
text(0,0,'a','Interpreter','latex', 'FontSize', 24);
    
load('../res/res_05b');
plt_fn1(2, rmse, 'Rank Thresholds', 'Kernel Size', 'Optimal Iterations',...
        {'10','20','30','40','50','60','70','80'},...
        {'[3,3]','[4,4]','[5,5]','[6,6]','[7,7]','[8,8]','[9,9]','[10,10]'});
line([3.5 8.5;3.5 8.5;3.5 3.5;8.5 8.5],[2.5 2.5;5.5 5.5;2.5 5.5;2.5 5.5],'linewidth',2,'color','r');    
text(0,0,'b','Interpreter','latex', 'FontSize', 24);

load('../res/res_05c');
plt_fn2(3, rmse_a_zoom, rmse_b_zoom, 'Rank Thresholds', 'RMSE',...
        {'40','44','48','52','56','60','64','68','72','76'},...
        {'kernel = [5,5], 50 iters','kernel = [6,6], 50 iters','kernel = [7,7], 50 iters','kernel = [5,5], optimal','kernel = [6,6], optimal','kernel = [7,7], optimal'});
text(0,0.53,'c','Interpreter','latex', 'FontSize', 24);

load('../res/res_05d');
plt_fn1(4, rmse, '$\mathcal{T_C}$ Rank Threshold', '$\mathcal{R_C}$ Rank Threshold', '',...
        {'10','20','30','40','50','60','70','80'},...
        {'10','20','30','40','50','60','70','80'});
text(0,0,'d','Interpreter','latex', 'FontSize', 24);

saveas(gcf, 'fig_05');
saveas(gcf, 'fig_05.svg');
saveas(gcf, 'fig_05','epsc');
%% Helper function
function plt_fn1(idx, data, xlbl, ylbl, zlbl, xtik, ytik)
    subplot(2,2,idx);
    
    imagesc(data, [0,1]);
    
    c   =   colorbar();
    c.Label.String          = 'RMSE';
    c.Label.Interpreter     = 'latex';
    c.TickLabelInterpreter  = 'latex';
    
    title(zlbl,'Interpreter','latex');
    
    xlabel(xlbl,'Interpreter','latex')
    xticks(1:8);
    xticklabels(xtik)
    
    ylabel(ylbl,'Interpreter','latex')
    yticks(1:8);
    yticklabels(ytik)
    
    axis square;
    set(gca,'FontSize',16)
    set(gca,'TickLabelInterpreter','latex');
end

function plt_fn2(idx, data1, data2, xlbl, ylbl, xtik, leg)
    subplot(2,2,idx);
    ax  =   gca;
    hold on;
    
    ax.ColorOrderIndex  =   1;
    plot(data1','linewidth',2);
    
    ax.ColorOrderIndex  =   1;
    plot(data2','linewidth',2,'linestyle','--');
    
    colorbar('Visible',false);
    
    xlabel(xlbl,'Interpreter','latex');
    xticks(2:2:20);
    xticklabels(xtik);
    xlim([0 21]);
    
    ylabel(ylbl,'Interpreter','latex');
    yticks(0:0.05:0.5);
    ylim([0 0.5]);
    
    legend(leg,'location','northwest','Interpreter','latex', 'FontSize', 10);
    legend('boxoff');
    set(ax,'TickLabelInterpreter','latex');    
    set(ax,'FontSize',16);
    grid on;
    axis square;
end
