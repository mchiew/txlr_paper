%% Fig 3a - RMSE vs R, across synthetic data slices, for H0 (VC), H2 (PRIMO) and H1+H2 (Proposed)

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
r       =   [50,50];
R       =   2:12;

% Perform under-sampled recovery and compute RMSE
rmse_H1H2   =   zeros(size(data,3),R(end));
out_H1H2    =   cell(R(end));
for i = R
    out_H1H2{i} = zeros(size(data));
    for z = 1:size(data,3)
        out_H1H2{i}(:,:,z,:,:)  =   reshape(admm_txlr(squeeze(data(:,:,z,:,:)).*mask(:,:,:,:,i), kernel, 50, r), 24,24,1,8,8);
        rmse_H1H2(z,i)          =   rmse(out_H1H2{i}(:,:,z,:,:),data(:,:,z,:,:));
    end
end

rmse_H2  =   zeros(size(data,3),R(end));
out_H2   =   cell(R(end));
for i = R
    out_H2{i} = zeros(size(data));
    for z = 1:size(data,3)
        out_H2{i}(:,:,z,:,:) =   reshape(admm_txlr(squeeze(data(:,:,z,:,:)).*mask(:,:,:,:,i), kernel, 100, [0 r(2)]), 24,24,1,8,8);
        rmse_H2(z,i)         =   rmse(out_H2{i}(:,:,z,:,:),data(:,:,z,:,:));
    end
end

rmse_H0  =   zeros(size(data,3),R(end));
out_H0   =   cell(R(end));
for i = R
    out_H0{i} = zeros(size(data));
    for z = 1:size(data,3)
        out_H0{i}(:,:,z,:,:)    =   reshape(admm_txlr(squeeze(data(:,:,z,:,:)).*mask(:,:,:,:,i), kernel, 100, [r(1) 0]), 24,24,1,8,8);
        rmse_H0(z,i)            =   rmse(out_H0{i}(:,:,z,:,:),data(:,:,z,:,:));
    end
end

save('../res/res_03','out_*','rmse_*');

%% Helper functions
function c = rmse(a,b)
    c   =  norm(a(:)-b(:))/norm(b(:));
end