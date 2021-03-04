%% Plot Fig 4

load('../res/res_04');

figure('Position',[0 0 800 400], 'PaperUnits', 'inches', 'PaperSize', 6.92*[1,0.5]);
ids = {'a', 'b'};

for i = 1:2
    plt_fn(i, cat(3,rmse_H0(:,:,i),rmse_H2(:,:,i),rmse_H1H2(:,:,i)), ids{i});
end

saveas(gcf,'fig_04');
saveas(gcf,'fig_04.svg');
saveas(gcf,'fig_04','epsc');

%% Helper function
function plt_fn(idx, data, id)
    subplot(1,2,idx);
    hold on;
    
    col = lines(3);

    title(sprintf('R=%d',4*idx),'Interpreter','latex');
    for i = 1:3
        for j = 1:3
            h((i-1)*3+j) =  scatter(sum(abs(data(:,j,i)-data(:,j,i)')<0.01,2).*randn(48,1)/200 + (j-1)*3 + i + j - 1,data(:,j,i),100,col(i,:),'filled');
        end
    end
    alpha(0.5);
    
    boxplot(reshape(permute(data,[1,3,2]),size(data,1),[]),'Positions',[1:3 5:7 9:11]);
    
    set(findobj(gca,'Type','Line'),'Linewidth',1);
    set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2 0.75])
    set(findobj(gca,'Tag','Outliers'),'Marker','None');
    set(gca,'TickLabelInterpreter','latex');
    
    xticks([2,6,10]);
    xticklabels(70:-10:50);
    xlabel('PSNR (dB)', 'Interpreter', 'latex');
    
    ylim([0 1]);
    ylabel('RMSE','Interpreter','latex');
    
    grid on;
    axis square;
    set(gca,'FontSize',16);
    set(gca,'TickLabelInterpreter','latex');
    
    if idx == 1
        %legend(h([1 4 7]),'PSNR = 70 dB ($\sigma=3.2\times 10^{-4}$)','PSNR = 60 dB ($\sigma=1.0\times 10^{-3}$)','PSNR = 50 dB ($\sigma=3.2\times 10^{-3}$)','Interpreter','latex')
        legend(h([1 4 7]),'VC','PRIMO','TxLR','Interpreter','latex');
        legend('boxoff');
    end   
    
    text(1,1.05,id,'Interpreter','latex', 'FontSize', 24);
end