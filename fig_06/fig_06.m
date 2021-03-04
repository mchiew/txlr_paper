%% Fig 06 - RMSE vs Constraints (H1 only [PRIMO], H2 only [No equiv.], H1 + H2 [TxLR])

% Add paths
addpath('../lib');

% Load data and crop to 24x24
q       =   matfile('../data/syn_data');
data    =   permute(reshape(double(crop_k(q.syn,[24,24])),24,24,[],8,8),[1,2,4,5,3]);

% Load 24x24 masks
q       =   matfile('../data/masks');
mask    =   reshape(q.mask24,24,24,1,[],12);

%% Recon

% Parameters
kernel  =   [5,5];
r       =   [50,50];
R       =   [8];
noise   =   10^(-3);

% Normalise data so that max value == 1
data    =   data/max(abs(data(:)));

% Loop over noise values and slices
out_H1      = zeros(size(data));
out_H2      = zeros(size(data));
out_H1H2    = zeros(size(data));
rmse_H1     = zeros(size(data,5),1);
rmse_H2     = zeros(size(data,5),1);
rmse_H1H2   = zeros(size(data,5),1);

for z = 1:size(data,5)
    dS  = data(:,:,:,:,z).*mask(:,:,:,:,R);
    N   = (noise/sqrt(2))*complex(randn(size(dS)),randn(size(dS))).*mask(:,:,:,:,R);

    out_H1(:,:,:,:,z)   = admm_txlr(dS+N, kernel, 100, [r(1) 0]);
    out_H2(:,:,:,:,z)   = admm_txlr(dS+N, kernel, 100, [0 r(2)]);
    out_H1H2(:,:,:,:,z) = admm_txlr(dS+N, kernel, 50, r);
    
    rmse_H1(z)      = rmse(out_H1(:,:,:,:,z), data(:,:,:,:,z));
    rmse_H2(z)      = rmse(out_H2(:,:,:,:,z), data(:,:,:,:,z));
    rmse_H1H2(z)    = rmse(out_H1H2(:,:,:,:,z), data(:,:,:,:,z));
end

save('../res/res_06','out_H*','rmse_H*');

%% Helper function
function c = rmse(a,b)
    c   =  norm(a(:)-b(:))/norm(b(:));
end
