%% Fig 2 - Singular value distribution for H1 (a) and H2 (b) unfolding

% Add paths
addpath('../lib');

% Load data and crop to 24x24
q       =   matfile('../data/syn_data');
data    =   reshape(double(crop_k(q.syn,[24,24])),24,24,[],8,8);

%% Get singular values
kernel  =   [5,5];
kdims   =   [24,24,8,8];
H1_fwd  =   @(x) fold_rx(Hankel_fwd(x, kernel, kdims));
H2_fwd  =   @(x) fold_tx(Hankel_fwd(x, kernel, kdims));

s1  =   zeros(prod(kernel)*kdims(3),size(data,3));
s2  =   zeros(prod(kernel)*kdims(4),size(data,3));
r1  =   zeros(prod(kernel)*kdims(3),size(data,3));
r2  =   zeros(prod(kernel)*kdims(4),size(data,3));

for z = 1:size(data,3)
    [~,s,~] =   svd(H1_fwd(squeeze(data(:,:,z,:,:))),'econ');
    s1(:,z) =   diag(s);
    [~,s,~] =   svd(H2_fwd(squeeze(data(:,:,z,:,:))),'econ');
    s2(:,z) =   diag(s);
    [~,s,~] =   svd(H1_fwd(randn(kdims)+1j*randn(kdims)),'econ');
    r1(:,z) =   diag(s);
    [~,s,~] =   svd(H1_fwd(randn(kdims)+1j*randn(kdims)),'econ');
    r2(:,z) =   diag(s);
end

s1 = s1./s1(1,:);
s2 = s2./s2(1,:);
r1 = r1./r1(1,:);
r2 = r2./r2(1,:);

save('../res/res_02','s1','s2','r1','r2');
