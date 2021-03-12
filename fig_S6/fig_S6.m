%% Fig S6 - Chi^2 heuristic 
% Add paths
addpath('../lib');

% Load data and crop to 24x24
q       =   matfile('../data/syn_data');
data    =   double(crop_k(squeeze(q.syn(:,:,24,:,:)),[24,24]));

% Load 24x24 masks
q       =   matfile('../data/masks');
mask    =   reshape(q.mask24,24,24,1,[],12);

%% Recon

% Parameters
kernel  =   [5,5];
r       =   [50,50];
R       =   [2,4,6,8];
noise   =   [10^(-3.5) 10^(-3) 10^(-2.5) 10^(-2)];

% Normalise data so that max value == 1
data    =   data/max(abs(data(:)));

% Loop over R
out     =   zeros([size(data) length(noise) length(R)]);
rmse	=   zeros([150 length(noise) length(R)]);
N_Chi2  =   zeros(length(noise),length(R));

for i = 1:length(R)
    dS =   data.*mask(:,:,:,:,R(i));
    for j = 1:length(noise)     
         N  =   (noise(j)/sqrt(2))*complex(randn(size(dS)),randn(size(dS))).*mask(:,:,:,:,R(i));
         
         [~, rmse(:,j,i)]           =   admm_txlr(dS+N, kernel, 150, r, noise(j)^2*eye(8), 'hard', data);
         [out(:,:,:,:,j,i), tmp]    =   admm_txlr(dS+N, kernel, 0, r, noise(j)^2*eye(8), 'hard', data);
         N_Chi2(j,i)                =   length(tmp);
    end
end

save('../res/res_S6','out','rmse','N_Chi2');
