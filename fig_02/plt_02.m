%% Plot Fig 2

load('../res/res_02');

figure('Position',[0 0 800 400], 'PaperUnits', 'inches', 'PaperSize', 6.92*[1,0.5]);

plt_fn(1, s1, r1, '$\mathcal{T_C}$-Unfolding', 'a');
plt_fn(2, s2, r2, '$\mathcal{R_C}$-Unfolding', 'b');

saveas(gcf,'fig_02');
saveas(gcf,'fig_02.svg');
saveas(gcf,'fig_02','epsc');

%% Helper functions
function plt_fn(idx, data1, data2, label, id)
    subplot(1,2,idx);
    hold on;
    
    col = lines(2);
    
    plot(mean(data1,2),'linewidth',2)
    plot(mean(data2,2),'linewidth',2)
    
    patch([1:200 200:-1:1], [min(data1,[],2)' fliplr(max(data1,[],2)')],1,'facecolor',col(1,:),'edgecolor','none');
    patch([1:200 200:-1:1], [min(data2,[],2)' fliplr(max(data2,[],2)')],1,'facecolor',col(2,:),'edgecolor','none');
    alpha(0.5);
    
    grid on;
    axis square;
    set(gca,'FontSize',16);
    set(gca,'TickLabelInterpreter','latex');
    
    xticks(50:50:200);
    yticks(0:0.1:1);
    
    xlabel('Singular Value Index','Interpreter','latex')
    ylabel('Normalized Singular Value','Interpreter','latex');
    legend({label,'Random Matrix'},'Interpreter','latex','Location','SouthEast')
    
    text(10,1.05,id,'Interpreter','latex', 'FontSize', 24);
end