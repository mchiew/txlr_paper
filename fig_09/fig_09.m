%% Fig 09 - Real Body Data
% Add paths
addpath('../lib');

% Load 6 datasets and crop to 24x24, with noise covariance matrices 
% Datasets 3-6 are asymmetrically sampled
q       =   matfile('../data/body_data');
data    =   zeros(24,24,8,8,6);
data(:,:,:,:,1) =   double(crop_k(q.kdata1,[24,24],true));
data(:,:,:,:,2) =   double(crop_k(q.kdata2,[24,24],true));
data(:,:,:,:,3) =   double(crop_k(q.kdata3(1:72,:,:,:),[24,24],true));
data(:,:,:,:,4) =   double(crop_k(q.kdata4(1:72,:,:,:),[24,24],true));
data(:,:,:,:,5) =   double(crop_k(q.kdata5(1:72,:,:,:),[24,24],true));
data(:,:,:,:,6) =   double(crop_k(q.kdata6(1:72,:,:,:),[24,24],true));

noise   =   zeros(8,8,6);
noise(:,:,1)    =   cov(q.noise1);
noise(:,:,2)    =   cov(q.noise2);
noise(:,:,3)    =   cov(q.noise3);
noise(:,:,4)    =   cov(q.noise4);
noise(:,:,5)    =   cov(q.noise5);
noise(:,:,6)    =   cov(q.noise6);

% Load 24x24 masks
q       =   matfile('../data/masks');
mask    =   reshape(q.mask24,24,24,1,[],12);

%% Recon

% Parameters
kernel  =   [5,5];
R       =   2:8;
r       =   5:5:50;

% Loop over datasets, R and rank
out_H1H2    =   zeros([size(data), length(R), length(r)]);
out_H2      =   zeros([size(data), length(R), length(r)]);
out_H0      =   zeros([size(data), length(R), length(r)]);

rmse_H1H2   =   zeros(size(data,5), length(R), length(r));
rmse_H2     =   zeros(size(data,5), length(R), length(r));
rmse_H0     =   zeros(size(data,5), length(R), length(r));


for i = 1:size(data,5)
    for j = 1:length(R)
        dS  =   data(:,:,:,:,i).*mask(:,:,:,:,R(j));
        for k = 1:length(r)       
            % H1H2: use Chi^2 criteria
            out_H1H2(:,:,:,:,i,j,k) =   admm_txlr(dS, kernel, 0, [r(k) r(k)], struct('noise',noise(:,:,i)));
            rmse_H1H2(i,j,k)        =   rmse(out_H1H2(:,:,:,:,i,j,k), data(:,:,:,:,i));

            % H2: use Chi^2 criteria
            out_H2(:,:,:,:,i,j,k)   =   admm_txlr(dS, kernel, 0, [0 r(k)], struct('noise',noise(:,:,i)));
            rmse_H2(i,j,k)          =   rmse(out_H2(:,:,:,:,i,j,k), data(:,:,:,:,i));

            % H0: use Chi^2 criteria
            out_H0(:,:,:,:,i,j,k)   =   reshape(admm_txlr(reshape(dS,24,24,[]), kernel, 0, [r(k) 0], struct('noise',diag(repmat(diag(noise(:,:,i)),8,1)))),24,24,8,8);
            rmse_H0(i,j,k)          =   rmse(out_H0(:,:,:,:,i,j,k), data(:,:,:,:,i));
        end
    end
end

save('../res/res_09','out_H1H2','out_H2','out_H0','rmse_H1H2','rmse_H2','rmse_H0','data');

%% Helper function
function c = rmse(a,b)
    c   =  norm(a(:)-b(:))/norm(b(:));
end
