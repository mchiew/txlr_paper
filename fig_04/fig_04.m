%% Fig 04 - RMSE vs SNR at R=4,8 for syn data all slices

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
R       =   [4, 8];
noise   =   [10^(-3.5), 10^(-3), 10^(-2.5)]; 

% Normalise data so that max value == 1
data    =   data/max(abs(data(:)));

% Loop over noise values and slices
out_H1H2    =   zeros([size(data) length(noise) length(R)]);
out_H2      =   zeros([size(data) length(noise) length(R)]);
out_H0      =   zeros([size(data) length(noise) length(R)]);
rmse_H1H2   =   zeros(size(data,5),length(noise),length(R));
rmse_H2     =   zeros(size(data,5),length(noise),length(R));
rmse_H0     =   zeros(size(data,5),length(noise),length(R));

for z = 1:size(data,5)
    for i = 1:length(R)
        dS  =   data(:,:,:,:,z).*mask(:,:,:,:,R(i));
        for j = 1:length(noise)
            N = (noise(j)/sqrt(2))*complex(randn(size(dS)),randn(size(dS))).*mask(:,:,:,:,R(i));
            out_H1H2(:,:,:,:,z,j,i) =   admm_txlr(dS+N, kernel, 50, r);
            out_H2(:,:,:,:,z,j,i)   =   admm_txlr(dS+N, kernel, 100, [0 r(2)]);
            out_H0(:,:,:,:,z,j,i)   =   reshape(admm_txlr(reshape(dS+N,24,24,[],1), kernel, 100, [r(1) 0]),24,24,8,8);
            
            rmse_H1H2(z,j,i)        =   rmse(out_H1H2(:,:,:,:,z,j,i),data(:,:,:,:,z));
            rmse_H2(z,j,i)          =   rmse(out_H2(:,:,:,:,z,j,i),  data(:,:,:,:,z));
            rmse_H0(z,j,i)          =   rmse(out_H0(:,:,:,:,z,j,i),  data(:,:,:,:,z));
        end
    end
end

save('../res/res_04','out_*','rmse_*');

%% Helper function
function c = rmse(a,b)
    c   =  norm(a(:)-b(:))/norm(b(:));
end
