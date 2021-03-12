% Plot Fig 6

figure('Position',[0 0 800 400], 'PaperUnits', 'inches', 'PaperSize', 6.92*[1,0.5]);
load('../res/res_06');

%% (a)

subplot('Position',[0.05 0.1 0.45 .9]); hold on;
set(gca,'FontSize',16);

scatter(sum(abs(rmse_H1-rmse_H1')<0.01,2).*randn(48,1)/250 + 1,rmse_H1,100,'filled','MarkerFaceColor',[0.5,0,0.5]);
scatter(sum(abs(rmse_H2-rmse_H2')<0.01,2).*randn(48,1)/250 + 2,rmse_H2,100,'filled')
scatter(sum(abs(rmse_H1H2-rmse_H1H2')<0.01,2).*randn(48,1)/250 + 3,rmse_H1H2,100,'filled')
alpha(0.5);

boxplot([rmse_H1 rmse_H2 rmse_H1H2]);
set(findobj(gca,'Type','Line'),'Linewidth',1);
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2 0.75])
set(findobj(gca,'Tag','Outliers'),'Marker','None');

axis square;
grid on;
xlim([0.5 3.5]);
ylim([0 0.4]);

set(gca,'TickLabelInterpreter','latex');
xticklabels({'\parbox{5em}{$\mathcal{T_C}$ only\\ }','\parbox{5em}{$\mathcal{R_C}$ only\\(PRIMO)}','\parbox{5em}{$\mathcal{T_C}+\mathcal{R_C}$\\(TxLR)}'});
ylabel('RMSE','Interpreter','latex');
yticks(0:0.05:0.4);

text(0,0.4,'a','Interpreter','latex', 'FontSize', 24);

%% (b)
addpath('../lib');

q       =   matfile('../data/syn_data');
data    =   permute(reshape(double(crop_k(q.syn,[24,24])),24,24,[],8,8),[1,2,4,5,3]);
data    =   data/max(abs(data(:)));

%clf(); hold on;
set(gca,'FontSize',16);

Rx  = 3;
Tx  = 4;
z   = 24;

% Plot magnitude images
plt_fn([0.5500 0.725 0.1 0.2], fftshift(ifft2(data(:,:,Rx,Tx,z))), [0 5E-3], 'Truth');
ylabel('Image','Interpreter','latex','Visible','on');
text(-5,-2.5,'b','Interpreter','latex', 'FontSize', 24);
plt_fn([0.6625 0.725 0.1 0.2], fftshift(ifft2(out_H1(:,:,Rx,Tx,z))), [0 5E-3], '$\mathcal{T_C}$ only');
plt_fn([0.7750 0.725 0.1 0.2], fftshift(ifft2(out_H2(:,:,Rx,Tx,z))), [0 5E-3], '$\mathcal{R_C}$ only');
plt_fn([0.8875 0.725 0.1 0.2], fftshift(ifft2(out_H1H2(:,:,Rx,Tx,z))), [0 5E-3], '$\mathcal{T_C}+\mathcal{R_C}$');

% Plot absolute difference from ground truth
plt_fn([0.6625 0.5 0.1 0.2], fftshift(ifft2(data(:,:,Rx,Tx,z)-out_H1(:,:,Rx,Tx,z))), [0 2E-3], '');
ylabel('diff','Interpreter','latex','Visible','on');
text(12,22,sprintf('%3.3f',rmse_H1(z)),'Color',0.99*[1,1,1],'FontSize', 16,'Interpreter','latex');
plt_fn([0.7750 0.5 0.1 0.2], fftshift(ifft2(data(:,:,Rx,Tx,z)-out_H2(:,:,Rx,Tx,z))), [0 2E-3], '');
text(12,22,sprintf('%3.3f',rmse_H2(z)),'Color',0.99*[1,1,1],'FontSize', 16,'Interpreter','latex');
plt_fn([0.8875 0.5 0.1 0.2], fftshift(ifft2(data(:,:,Rx,Tx,z)-out_H1H2(:,:,Rx,Tx,z))), [0 2E-3], '');
text(12,22,sprintf('%3.3f',rmse_H1H2(z)),'Color',0.99*[1,1,1],'FontSize', 16,'Interpreter','latex');

% Get Tx sensitivities
kernel      =   [5,5];
imsize      =   [64,64];
thresh      =   50;

mask        =   sum(sum(abs(ifft2(padarray(data(:,:,:,:,z),0.5*[imsize(1)-size(data,1),imsize(2)-size(data,2),0,0]))).^2,3),4).^0.5;
mask        =   fftshift(mask > 0.05*max(mask(:)));
sens_ref    =   tx_espirit(data(:,:,:,:,z), imsize, kernel, thresh).*mask;
sens_H1     =   tx_espirit(out_H1(:,:,:,:,z), imsize, kernel, thresh).*mask;
sens_H2     =   tx_espirit(out_H2(:,:,:,:,z), imsize, kernel, thresh).*mask;
sens_H1H2   =   tx_espirit(out_H1H2(:,:,:,:,z), imsize, kernel, thresh).*mask;

% Phase align for RMSE comparison (since absolute phase is irrelevant)
sens_H1     =   phs_align(sens_ref, sens_H1);
sens_H2     =   phs_align(sens_ref, sens_H2);
sens_H1H2   =   phs_align(sens_ref, sens_H1H2);

% Plot magnitude Tx sens
plt_fn([0.5500 0.225 0.1 0.2], sens_ref(:,:,Tx), [0 1], '');
ylabel('Tx sens','Interpreter','latex','Visible','on');
text(-13.33,-6.66,'c','Interpreter','latex', 'FontSize', 24);
plt_fn([0.6625 0.225 0.1 0.2], sens_H1(:,:,Tx), [0 1], '');
plt_fn([0.7750 0.225 0.1 0.2], sens_H2(:,:,Tx), [0 1], '');
plt_fn([0.8875 0.225 0.1 0.2], sens_H1H2(:,:,Tx), [0 1], '');

% Plot absolute difference Tx sens from ground truth
plt_fn([0.6625 0 0.1 0.2], sens_ref(:,:,Tx)-sens_H1(:,:,Tx), [0 1], '');
ylabel('diff','Interpreter','latex','Visible','on');
text(34,58,sprintf('%3.3f',rmse(sens_H1(:,:,Tx),sens_ref(:,:,Tx))),'Color',0.99*[1,1,1],'FontSize', 16,'Interpreter','latex');
plt_fn([0.7750 0 0.1 0.2], sens_ref(:,:,Tx)-sens_H2(:,:,Tx), [0 1], '');
text(34,58,sprintf('%3.3f',rmse(sens_H2(:,:,Tx),sens_ref(:,:,Tx))),'Color',0.99*[1,1,1],'FontSize', 16,'Interpreter','latex');
plt_fn([0.8875 0 0.1 0.2], sens_ref(:,:,Tx)-sens_H1H2(:,:,Tx), [0 1], '');
text(34,58,sprintf('%3.3f',rmse(sens_H1H2(:,:,Tx),sens_ref(:,:,Tx))),'Color',0.99*[1,1,1],'FontSize', 16,'Interpreter','latex');

saveas(gcf,'fig_06');
saveas(gcf,'fig_06.svg');
saveas(gcf,'fig_06','epsc');

%% Helper functions
function plt_fn(pos, img, clim, label)
    axes('Position', pos);
    imagesc(abs(img),clim);
    axis square; axis off
    title(label,'Interpreter','latex');
    xticklabels '';
    yticklabels '';
    set(gca,'FontSize',16);
end

function b = phs_align(a,b)
    b   =   b.*exp(-1j*angle(sum(conj(a).*b,3)));
end

function err = rmse(a,b)
    err =   norm(a(:)-b(:))/norm(b(:));
end
