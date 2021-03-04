%% Plot Fig 9

% Load data
load('../res/res_09');

figure('Position',[0 0 800 800], 'PaperUnits', 'inches', 'PaperSize', 6.92*[1,1]);
ids = 'ab';

%% Plot Fig 9
% Add paths
addpath('../lib');

% Load data

q       =   matfile('../data/body_data');
data    =   zeros(24,24,8,8,6);
data(:,:,:,:,1) =   double(crop_k(q.kdata1,[24,24],true));
data(:,:,:,:,2) =   double(crop_k(q.kdata2,[24,24],true));
data(:,:,:,:,3) =   double(crop_k(q.kdata3(1:72,:,:,:),[24,24],true));
data(:,:,:,:,4) =   double(crop_k(q.kdata4(1:72,:,:,:),[24,24],true));
data(:,:,:,:,5) =   double(crop_k(q.kdata5(1:72,:,:,:),[24,24],true));
data(:,:,:,:,6) =   double(crop_k(q.kdata6(1:72,:,:,:),[24,24],true));

% Subject, dataset and coil to plot
subj        =   4;
c           =   3;
rank        =   50;

% ESPIRiT Parameters
kernel      =   [5,5];
imsize      =   [64,64];
thresh      =   rank;

lbls        =   {'Truth','VC','PRIMO','TxLR';'','','',''};

R           =   [4,8];
for i = 1:2
    % Compute and mask Tx-sensitivities
    mask        =   sum(sum(abs(ifft2(padarray(data(:,:,:,:,subj),0.5*[imsize(1)-size(data,1),imsize(2)-size(data,2),0,0]))).^2,3),4).^0.5;
    mask        =   fftshift(mask > 0.05*max(mask(:)));
    sens_ref    =   tx_espirit(data(:,:,:,:,subj), imsize, kernel, thresh).*mask;
    sens_H1H2   =   tx_espirit(out_H1H2(:,:,:,:,subj,R(i)-1,rank/5), imsize, kernel, thresh).*mask;
    sens_H2     =   tx_espirit(out_H2(:,:,:,:,subj,R(i)-1,rank/5), imsize, kernel, thresh).*mask;
    sens_H0     =   tx_espirit(out_H0(:,:,:,:,subj,R(i)-1,rank/5), imsize, kernel, thresh).*mask;

    % Phase align for RMSE comparison (since absolute phase is irrelevant)
    sens_H1H2   =   phs_align(sens_ref, sens_H1H2);
    sens_H2     =   phs_align(sens_ref, sens_H2);
    sens_H0     =   phs_align(sens_ref, sens_H0);

    % Plot magnitude sensitivities
    plt_fn([0.0500 0.735-0.5*(i-1) 0.2 0.2], abs(sens_ref(:,:,c)),   lbls{i,1});
    ylabel(sprintf('R=%d',R(i)),'Interpreter','latex','Visible','on');
    text(-13.33,-6.66,ids(i),'Interpreter','latex', 'FontSize', 24);
    plt_fn([0.2750 0.735-0.5*(i-1) 0.2 0.2], abs(sens_H0(:,:,c)),    lbls{i,2});
    plt_fn([0.5000 0.735-0.5*(i-1) 0.2 0.2], abs(sens_H2(:,:,c)),    lbls{i,3});
    plt_fn([0.7250 0.735-0.5*(i-1) 0.2 0.2], abs(sens_H1H2(:,:,c)),  lbls{i,4});

    cbar = colorbar();
    cbar.Position = [0.935 1.01-i*0.5 0.025 0.2];
    
    % Plot absolute differences
    plt_fn([0.2750 0.51-0.5*(i-1) 0.2 0.2], sens_ref(:,:,c)-sens_H0(:,:,c),     '');
    text(43,60,sprintf('%3.3f', rmse(sens_H0(:,:,c),sens_ref(:,:,c))),'Color',0.99*[1,1,1],'FontSize', 20,'Interpreter','latex');
    ylabel(sprintf('diff'),'Interpreter','latex','Visible','on');
    plt_fn([0.5000 0.51-0.5*(i-1) 0.2 0.2], sens_ref(:,:,c)-sens_H2(:,:,c),     '');
    text(43,60,sprintf('%3.3f', rmse(sens_H2(:,:,c),sens_ref(:,:,c))),'Color',0.99*[1,1,1],'FontSize', 20,'Interpreter','latex');
    plt_fn([0.7250 0.51-0.5*(i-1) 0.2 0.2], sens_ref(:,:,c)-sens_H1H2(:,:,c),   '');
    text(43,60,sprintf('%3.3f', rmse(sens_H1H2(:,:,c),sens_ref(:,:,c))),'Color',0.99*[1,1,1],'FontSize', 20,'Interpreter','latex');
end

saveas(gcf,'fig_09');
saveas(gcf,'fig_09.svg');
saveas(gcf,'fig_09','epsc');

%% Helper function
function plt_fn(pos, img, label)
    axes('Position', pos);
    imagesc(abs(rot90(img)), [0,1]);
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