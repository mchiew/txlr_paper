%% Fig 5d - 2D RMSE grid, Kernel size vs rank threshold zoomed

% Add paths
addpath('../lib');

% Load data and crop to 24x24
q       =   matfile('../data/syn_data');
data    =   double(crop_k(squeeze(q.syn(:,:,24,:,:)),[24,24]));

% Load 24x24 R=8 masks
q       =   matfile('../data/masks');
mask    =   reshape(q.mask24(:,:,:,8),24,24,1,[]);

%% Recon

% Kernel sizes
kernel = [5,5];

% Rank thresholds
r       =   10:10:80;

% Perform under-sampled recovery and compute RMSE
rmse    =   zeros(length(kernels), length(r));
for i = 1:length(r)
    for j = 1:length(r)
        out         =   admm_txlr(data.*mask, kernel, 50, [r(i) r(j)]);
        rmse(i,j)   =   norm(out(:)-data(:))/norm(data(:));
    end
end
save('../res/res_05d','rmse','kernel','r');
