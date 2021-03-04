%% Fig S6 - Variable density sampling 

% Add paths
addpath('../lib');

% Load data and crop to 24x24
q       =   matfile('../data/syn_data');
data    =   reshape(double(crop_k(q.syn,[24,24])),24,24,[],8,8);

% Load 24x24 masks
q       =   matfile('../data/mask_VD');
mask    =   reshape(q.maskVDR8,24,24,1,8);

%% Recon

% Parameters
kernel  =   [5,5];
r       =   50;
R       =   8;

% Prepare data
d0      =   squeeze(data(:,:,24,:,:));
dS      =   d0.*mask;
opts.truth  =   d0;

% Perform under-sampled recovery and compute RMSE
[~, rmse_txlr] = admm_txlr(dS, kernel, 50, [r, r], opts);

save('../res/res_S6','rmse_*');