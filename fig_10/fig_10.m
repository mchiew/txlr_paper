%% Fig 10 - Real Brain Data
% Add paths
addpath('../lib');

% Load bottom 20 slices of brain dataset and crop to 24x24, with noise covariance matrices 
q       =   matfile('../data/brain_data');
data    =   permute(reshape(double(crop_k(q.kdata(:,:,7:26,:,:),[24,24])),24,24,20,32,8),[1,2,4,5,3]);
noise   =   cov(q.noise);

% Load 24x24 masks
q       =   matfile('../data/masks');
mask    =   reshape(q.mask24,24,24,1,[],12);

%% Recon

% Parameters
kernel  =   [5,5];
R       =   2:8;
r       =   [50 50];

% Loop over datasets and accleration factor
out_txlr    =   zeros([size(data), length(R)]);
out_primo      =   zeros([size(data), length(R)]);
out_vc      =   zeros([size(data), length(R)]);

rmse_txlr   =   zeros(size(data,5), length(R));
rmse_primo     =   zeros(size(data,5), length(R));
rmse_vc     =   zeros(size(data,5), length(R));


opts    =   struct('noise',noise);

%%
%for j = 1:length(R)
    %for i = 1:size(data,5)
    for j = 4:length(R)
        disp(j)
    for i = 1:size(data,5)
        dS  =   data(:,:,:,:,i).*mask(:,:,:,:,R(j));   
        % TxLR: use Chi^2 criteria
        out_txlr(:,:,:,:,i,j)   =   admm_txlr(dS, kernel, 0, r, opts);
        rmse_txlr(i,j)          =   rmse(out_txlr(:,:,:,:,i,j), data(:,:,:,:,i));

        % PRIMO: use Chi^2 criteria
        out_primo(:,:,:,:,i,j)     =   admm_txlr(dS, kernel, 0, [0 r(2)], opts);
        rmse_primo(i,j)            =   rmse(out_primo(:,:,:,:,i,j), data(:,:,:,:,i));

        % VC: use Chi^2 criteria
        out_vc(:,:,:,:,i,j)     =   reshape(admm_txlr(reshape(dS,24,24,[]), kernel, 0, [r(1) 0], struct('noise',diag(repmat(diag(noise),8,1)))),24,24,32,8);
        rmse_vc(i,j)            =   rmse(out_vc(:,:,:,:,i,j), data(:,:,:,:,i));
    end
end
%%
save('../res/res_10','out_txlr','out_primo','out_vc','rmse_txlr','rmse_primo','rmse_vc','data');

%% Helper function
function c = rmse(a,b)
    c   =  norm(a(:)-b(:))/norm(b(:));
end
