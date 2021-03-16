%% Fig S3 - RMSE vs k-FOV at R=8, PSNR = for syn data all slices

% Add paths
addpath('../lib');

% Load data and crop to 24x24
q       =   matfile('../data/syn_data');
data18  =   permute(reshape(double(crop_k(q.syn,[18,18])),18,18,[],8,8),[1,2,4,5,3]);
data24  =   permute(reshape(double(crop_k(q.syn,[24,24])),24,24,[],8,8),[1,2,4,5,3]);
data36  =   permute(reshape(double(crop_k(q.syn,[36,36])),36,36,[],8,8),[1,2,4,5,3]);
data48  =   permute(reshape(double(crop_k(q.syn,[48,48])),48,48,[],8,8),[1,2,4,5,3]);

% Load 24x24 masks
q       =   matfile('../data/masks');
mask18  =   reshape(q.mask18,18,18,1,[],12);
mask24  =   reshape(q.mask24,24,24,1,[],12);
mask36  =   reshape(q.mask36,36,36,1,[],12);
mask48  =   reshape(q.mask48,48,48,1,[],12);

%% Recon

% Parameters
kernel  =   [5,5];
r       =   [50,50];
R       =   [2,4,6,8];
noise   =   10^(-3);

% Normalise data so that max value == 1
data18    =   data18/max(abs(data18(:)));
data24    =   data24/max(abs(data24(:)));
data35    =   data36/max(abs(data36(:)));
data48    =   data48/max(abs(data48(:)));

% Loop over R
out_18    = zeros([size(data18) length(R)]);
out_24    = zeros([size(data24) length(R)]);
out_36    = zeros([size(data36) length(R)]);
out_48    = zeros([size(data48) length(R)]);

rmse_18   = zeros([size(data18,5) length(R)]);
rmse_24   = zeros([size(data24,5) length(R)]);
rmse_36   = zeros([size(data36,5) length(R)]);
rmse_48   = zeros([size(data48,5) length(R)]);

for z = 1:size(data24,5)
    for i = 1:length(R)
        dS18    =   data18(:,:,:,:,z).*mask18(:,:,:,:,R(i));
        N18     =   (noise/sqrt(2))*complex(randn(size(dS18)),randn(size(dS18))).*mask18(:,:,:,:,R(i));
        
        out_18(:,:,:,:,z,i) =   admm_txlr(dS18+N18, kernel, 50, r);
        rmse_18(z,i)    = rmse(out_18(:,:,:,:,z,i), data18(:,:,:,:,z));

        dS24    =   data24(:,:,:,:,z).*mask24(:,:,:,:,R(i));
        N24     =   (noise/sqrt(2))*complex(randn(size(dS24)),randn(size(dS24))).*mask24(:,:,:,:,R(i));
        
        out_24(:,:,:,:,z,i) =   admm_txlr(dS24+N24, kernel, 50, r);
        
        dS36    =   data36(:,:,:,:,z).*mask36(:,:,:,:,R(i));
        N36     =   (noise/sqrt(2))*complex(randn(size(dS36)),randn(size(dS36))).*mask36(:,:,:,:,R(i));
        
        out_36(:,:,:,:,z,i) =   admm_txlr(dS36+N36, round(kernel*1.5), 50, round(r*1.5));
        
        dS48    =   data48(:,:,:,:,z).*mask48(:,:,:,:,R(i));
        N48     =   (noise/sqrt(2))*complex(randn(size(dS48)),randn(size(dS48))).*mask48(:,:,:,:,R(i));
        
        out_48(:,:,:,:,z,i) =   admm_txlr(dS48+N48, round(kernel*2), 50, round(r*2));           
            
        rmse_24(z,i)    = rmse(out_24(:,:,:,:,z,i), data24(:,:,:,:,z));
        rmse_36(z,i)    = rmse(out_36(:,:,:,:,z,i), data36(:,:,:,:,z));
        rmse_48(z,i)    = rmse(out_48(:,:,:,:,z,i), data48(:,:,:,:,z));
    end
end

save('../res/res_S3','out_24','out_36','out_48','rmse_24','rmse_36','rmse_48');

%% Helper function
function c = rmse(a,b)
    c   =  norm(a(:)-b(:))/norm(b(:));
end
