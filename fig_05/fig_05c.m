%% Fig 5c - 2D RMSE grid, Kernel size vs rank threshold zoomed

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
kernels =   {[5,5],[6,6],[7,7]};

% Rank thresholds
r       =   {[40,40],[42,42],[44,44],[46,46],[48,48],...
             [50,50],[52,52],[54,54],[56,56],[58,58],...
             [60,60],[62,62],[64,64],[66,66],[68,68],...
             [70,70],[72,72],[74,74],[76,76],[78,78]};

% Perform under-sampled recovery and compute RMSE
rmse_a_zoom    =   zeros(length(kernels), length(r));
for i = 1:length(kernels)
    for j = 1:length(r)
        if all(r{j} < prod(kernels{i})*[size(data,3) size(data,4)])
            out                 =   admm_txlr(data.*mask, kernels{i}, 50, r{j});
            rmse_a_zoom(i,j)    =   norm(out(:)-data(:))/norm(data(:));
        else
            rmse_a_zoom(i,j)    =   inf;
        end
    end
end

% Perform under-sampled recovery and compute RMSE
rmse_b_zoom    =   zeros(length(kernels), length(r));
for i = 1:length(kernels)
    for j = 1:length(r)
        if all(r{j} < prod(kernels{i})*[size(data,3) size(data,4)])
            [out,err]           =   admm_txlr(data.*mask, kernels{i}, 200, r{j}, [], 'hard', data);
            rmse_b_zoom(i,j)    =   min(err);
        else
            rmse_b_zoom(i,j)    =   inf;
        end
    end
end
save('../res/res_05c','rmse_a_zoom','rmse_b_zoom','kernels','r');
