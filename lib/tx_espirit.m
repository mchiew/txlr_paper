function sens = tx_espirit(calib, imsize, kernel, eig_thresh)
%
% sens = tx_espirit(calib, imsize, kernel, eig_thresh)
%
% Inputs:
%           calib:      [NKx, NKy, NRx, NTx] complex Rx/Tx k-space
%           imsize:     [Nx, Ny] image dimensions for output
%           kernel:     [kx, ky] kernel size
%                       (default [5,5])
%           eig_thresh: scalar threshold for eigenvalues 
%                       if < 1, interpreted as s ≤ s(1)*eig_thresh
%                       if > 1, interpreted as s_r, r ≤ eig_thresh
%                       (default 0.02)
%
% Outputs:
%           sens:       [Nx, Ny, NTx] Tx sensitivity maps
%
% Implementation of relative Tx mapping based on 
% ESPIRiT (MRM 2014) and PRIMO (MRM 2015)
%
% MChiew (mark.chiew@ndcn.ox.ac.uk)

% Default params
if nargin < 4
    eig_thresh  =   0.02;
end
if nargin < 3
    kernel      =   [5,5];
end

% Generate calibration matrix
% Tx vertically concatenated, Rx horizontally
dims    =   size(calib);
H       =   fold_tx(Hankel_fwd(calib, kernel, dims));

% For receive sensitivity mapping, use fold_rx instead
% H_rx   =   fold_rx(Hankel_fwd(calib, kernel, dims));

% Get left singular (kernel x coil) vectors, above singular value threshold
[U,S,~] =   svd(H, 'econ');
if eig_thresh < 1
    U       =   U(:,diag(S) > S(1)*eig_thresh);
else
    U       =   U(:,1:round(eig_thresh));
end
    

% Get zero-padded space x coil x component kernel images
U       =   reshape(U,kernel(1),kernel(2),dims(4),[]);
U       =   padarray(U,[ceil((imsize-kernel)/2) 0 0],'pre');
U       =   padarray(U,[floor((imsize-kernel)/2) 0 0],'post');
tmp     =   fftshift(fftshift(ifft2(fftshift(fftshift(U,1),2)),1),2);

% Perform SVD on coil x component matrices, voxelwise
% Keep the first component
sens    =   zeros([imsize dims(4)]);
for i = 1:imsize(1)
    for j = 1:imsize(2)
        [U,~,~]     =   svd(squeeze(tmp(i,j,:,:)), 'econ');        
        sens(i,j,:) =   U(:,1);
    end
end

% Rotate phase relative to first coil
sens    =   sens.*exp(-1j*angle(sens(:,:,1)));
