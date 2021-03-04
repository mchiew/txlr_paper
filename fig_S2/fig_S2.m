%% Fig S1 - Checking optimal rank thresholds, iterations for TxLR, PRIMO and VC

% Add paths
addpath('../lib');

% Load data and crop to 24x24
q       =   matfile('../data/syn_data');
data    =   reshape(double(crop_k(q.syn,[24,24])),24,24,[],8,8);

% Load 24x24 masks
q       =   matfile('../data/masks');
mask    =   reshape(q.mask24,24,24,1,[],12);

%% Recon

% Parameters
kernel  =   [5,5];
r       =   1:100;
R       =   8;

% Prepare data
d0      =   squeeze(data(:,:,24,:,:));
dS      =   d0.*mask(:,:,:,:,R);
opts.truth  =   d0;

% Perform under-sampled recovery and compute RMSE
rmse_H1H2 = zeros(100);
for i = r
    [~, err] = admm_txlr(dS, kernel, 100, [i i], opts);
    rmse_H1H2(:,i) = err;
end

rmse_H2 = zeros(100);
for i = r
    [~, err] = admm_txlr(dS, kernel, 100, [0 i], opts);
    rmse_H2(:,i) = err;
end

rmse_H0 = zeros(100);
for i = r
    [~, err] = admm_txlr(reshape(dS,24,24,[]), kernel, 100, [i 0], opts);
    rmse_H0(:,i) = err;
end

save('res_S1','rmse_*');