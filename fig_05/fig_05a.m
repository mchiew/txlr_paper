%% Fig 5a - 2D RMSE grid, Kernel size vs rank threshold (fixed iters)

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
kernels =   {[3,3],[4,4],[5,5],[6,6],[7,7],[8,8],[9,9],[10,10]};

% Rank thresholds
r       =   {[10,10],[20,20],[30,30],[40,40],[50,50],[60,60],[70,70],[80,80]};

% Perform under-sampled recovery and compute RMSE
rmse    =   zeros(length(kernels), length(r));
for i = 1:length(kernels)
    for j = 1:length(r)
        if all(r{j} < prod(kernels{i})*[size(data,3) size(data,4)])
            out         =   admm_txlr(data.*mask, kernels{i}, 50, r{j});
            rmse(i,j)   =   norm(out(:)-data(:))/norm(data(:));
        else
            rmse(i,j)   =   inf;
        end
    end
end
save('../res/res_05a','rmse','kernels','r');
