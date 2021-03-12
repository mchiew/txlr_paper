%% Plot Fig 09

% Load data
load('../res/res_09');

% Add paths
addpath('../lib');

% Create figure
figure('Position',[0 0 800 400], 'PaperUnits', 'inches', 'PaperSize', 6.92*[1,0.5]);
ids = 'abc';

%% Plot Fig 09a
subplot('Position',[0.05 0.1 0.45 .9]); hold on;
ax = gca;
set(ax,'FontSize',16);
R = [1,3,5,7];

for i = 1:4
    ax.ColorOrderIndex = 1;
    scatter(sum(abs(rmse_vc(:,R(i))-rmse_vc(:,R(i))')<0.005,2).*randn(20,1)/250 + (i-1)*3+1 + (i-1),rmse_vc(:,R(i)),100,'filled')
    scatter(sum(abs(rmse_primo(:,R(i))-rmse_primo(:,R(i))')<0.005,2).*randn(20,1)/250 + (i-1)*3+2 + (i-1),rmse_primo(:,R(i)),100,'filled')
    scatter(sum(abs(rmse_txlr(:,R(i))-rmse_txlr(:,R(i))')<0.005,2).*randn(20,1)/250 + (i-1)*3+3 + (i-1),rmse_txlr(:,R(i)),100,'filled')
end
alpha(0.5);

boxplot([rmse_vc(:,1) rmse_primo(:,1) rmse_txlr(:,1) rmse_vc(:,3) rmse_primo(:,3) rmse_txlr(:,3) rmse_vc(:,5) rmse_primo(:,5) rmse_txlr(:,5) rmse_vc(:,7) rmse_primo(:,7) rmse_txlr(:,7)],'positions',[1:3 5:7 9:11 13:15]);
set(findobj(gca,'Type','Line'),'Linewidth',1);
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2 0.75])
set(findobj(gca,'Tag','Outliers'),'Marker','None');

axis square;
grid on;
xlim([0.5 15.5]);
ylim([0 1]);
yticks(0:0.1:1);
xticks([1:3 5:7 9:11 13:15])

set(gca,'TickLabelInterpreter','latex');
xticklabels({'','R=2','','','R=4','','','R=6','','','R=8',''})
ylabel('RMSE','Interpreter','latex');
legend({'VC','PRIMO','TxLR'},'Interpreter','latex','Location','NorthWest');
text(-1.5,1,ids(1),'Interpreter','latex', 'FontSize', 32);

%% Plot Fig 09b

% Load data
q       =   matfile('../data/brain_data');
data    =   permute(reshape(double(q.kdata(:,:,7:26,:,:)),48,48,20,32,8),[1,2,4,5,3]);


i = 7;

% ESPIRiT Parameters
kernel      =   [5,5];
imsize      =   [64,64];
thresh      =   50;

% Compute masks
for z = 1:20
    mask(:,:,z)         =   sum(sum(abs(ifft2(padarray(data(:,:,:,:,z),0.5*[imsize(1)-size(data,1),imsize(2)-size(data,2),0,0]))).^2,3),4).^0.5;
    mask(:,:,z)         =   fftshift(mask(:,:,z) > 0.1*max(reshape(mask(:,:,z),[],1)));
end

% Crop data down to 24x24 matrix 
data    =   reshape(crop_k(data, [24,24]),[24,24,32,8,20]);

% Compute Tx-sensitivities
for z = 1:20
    
    sens_ref(:,:,:,z)   =   tx_espirit(data(:,:,:,:,z), imsize, kernel, thresh).*mask(:,:,z);
    sens_txlr(:,:,:,z)  =   tx_espirit(out_txlr(:,:,:,:,z,i), imsize, kernel, thresh).*mask(:,:,z);
    sens_primo(:,:,:,z) =   tx_espirit(out_primo(:,:,:,:,z,i), imsize, kernel, thresh).*mask(:,:,z);
    sens_vc(:,:,:,z)    =   tx_espirit(out_vc(:,:,:,:,z,i), imsize, kernel, thresh).*mask(:,:,z);
    
    % Phase align for RMSE comparison (since absolute phase is irrelevant)
    sens_txlr(:,:,:,z)  =   phs_align(sens_ref(:,:,:,z), sens_txlr(:,:,:,z));
    sens_primo(:,:,:,z) =   phs_align(sens_ref(:,:,:,z), sens_primo(:,:,:,z));
    sens_vc(:,:,:,z)    =   phs_align(sens_ref(:,:,:,z), sens_vc(:,:,:,z));
end
%%
% Plot options
c = 8;
z = 16;
lbls        =   {'Truth','VC','PRIMO','TxLR';'','','',''};


% Plot magnitude sensitivities
plt_mag([0.5500 0.725 0.1 0.2], sens_ref(:,:,c,z),   lbls{1,1});
ylabel('Mag','Interpreter','latex','Visible','on');
text(-13.33,-6.66,ids(2),'Interpreter','latex', 'FontSize', 32);
plt_mag([0.6625 0.725 0.1 0.2], sens_vc(:,:,c,z),    lbls{1,2});
plt_mag([0.7750 0.725 0.1 0.2], sens_primo(:,:,c,z),    lbls{1,3});
plt_mag([0.8875 0.725 0.1 0.2], sens_txlr(:,:,c,z),  lbls{1,4});

cbar = colorbar();
cbar.Position = [0.585 0.5 0.025 0.2];

% Plot absolute differences
plt_mag([0.6625 0.5 0.1 0.2], abs(sens_ref(:,:,c,z))-abs(sens_vc(:,:,c,z)),     '');
text(34,58,sprintf('%3.3f', rmse(sens_vc(:,:,c,z),sens_ref(:,:,c,z))),'Color',0.99*[1,1,1],'FontSize', 16,'Interpreter','latex');
ylabel('diff','Interpreter','latex','Visible','on');
plt_mag([0.7750 0.5 0.1 0.2], abs(sens_ref(:,:,c,z))-abs(sens_primo(:,:,c,z)),     '');
text(34,58,sprintf('%3.3f', rmse(sens_primo(:,:,c,z),sens_ref(:,:,c,z))),'Color',0.99*[1,1,1],'FontSize', 16,'Interpreter','latex');
plt_mag([0.8875 0.5 0.1 0.2], abs(sens_ref(:,:,c,z))-abs(sens_txlr(:,:,c,z)),   '');
text(34,58,sprintf('%3.3f', rmse(sens_txlr(:,:,c,z),sens_ref(:,:,c,z))),'Color',0.99*[1,1,1],'FontSize', 16,'Interpreter','latex');

% Plot phase sensitivities
plt_phs([0.5500 0.225 0.1 0.2], angle(sens_ref(:,:,c,z)),   lbls{2,1});
ylabel('Phs','Interpreter','latex','Visible','on');
text(-13.33,-6.66,ids(3),'Interpreter','latex', 'FontSize', 32);
plt_phs([0.6625 0.225 0.1 0.2], angle(sens_vc(:,:,c,z)),    lbls{2,2});
plt_phs([0.7750 0.225 0.1 0.2], angle(sens_primo(:,:,c,z)),    lbls{2,3});
plt_phs([0.8875 0.225 0.1 0.2], angle(sens_txlr(:,:,c,z)),  lbls{2,4});

cbar = colorbar();
cbar.Position = [0.585 0.0 0.025 0.2];

% Plot phase differences
plt_phs([0.6625 0.0 0.1 0.2], angle(sens_ref(:,:,c,z))-angle(sens_vc(:,:,c,z)),     '');
ylabel('diff','Interpreter','latex','Visible','on');
plt_phs([0.7750 0.0 0.1 0.2], angle(sens_ref(:,:,c,z))-angle(sens_primo(:,:,c,z)),     '');
plt_phs([0.8875 0.0 0.1 0.2], angle(sens_ref(:,:,c,z))-angle(sens_txlr(:,:,c,z)),   '');


saveas(gcf,'fig_09');
saveas(gcf,'fig_09.svg');
saveas(gcf,'fig_09','epsc');

%% Helper function
function plt_mag(pos, img, label)
    subplot('Position',pos);
    imagesc(flip(abs(img),1),[0 1]);colormap(gca,parula);
    axis square; axis off
    title(label,'Interpreter','latex')
    xticklabels '';
    yticklabels '';
    set(gca,'FontSize',16);
end

function plt_phs(pos, img, label)
    subplot('Position',pos);
    imagesc(flip(img,1),[-pi,pi]);colormap(gca,hsv);
    axis square; axis off
    title(label,'Interpreter','latex')
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
