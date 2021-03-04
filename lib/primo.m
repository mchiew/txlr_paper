function [z varargout] = primo(k_sens, kernel, niters, r, noise, truth)
% z = primo(k_sens, kernel, niters, r, [noise])
% [z rmse] = primo(k_sens, kernel, niters, r, [noise], truth)
%
%   Inputs: 
%           k_sens  -   under-sampled [kx,ky,Rx,Tx] k-space
%           kernel  -   kernel size [Nx, Ny]
%           niters  -   number of ADMM iterations (suggested: 10)
%           r       -   rank constraint
%           [noise] -   optional, receiver noise covariance matrix
%                       defaults to 1E-6*eye(Rx), only used if niters = 0
%                       which uses a Chi^2/N = 1 criteria for stopping
%           [truth] -   optional, ground truth data
%                       if present, varargout{1} = rmse at each iteration
%
%   Output:
%           z       -   recovered k-space
%
%   Performs the "PRIMO" approach for pTx mapping
%   The multi-dimensional Hankel array is reshaped so that Rx coils are vertically stacked
%   and Tx sens are horizontally concatenated
%   Uses simple alternating projection approach
%   
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
% v1.0 (21/04/20) - Initial version
% v1.1 (27/04/20) - Added ground truth rmse option

if nargin < 6
    truth = [];
end
if nargin < 5
    noise = [];
end

% sampling mask
mask = k_sens ~= 0;

% setup helper functions
kdims = size(k_sens); 
hdims = size(Hankel_fwd(k_sens, kernel, kdims)); 
[~,N] = Hankel_adj(zeros(hdims), kernel, kdims);

H1_fwd = @(x) fold_rx(Hankel_fwd(x, kernel, kdims));
H1_adj = @(x) Hankel_adj(unfold_rx(x, hdims), kernel, kdims);

% iniitalise variables
z   = k_sens;
i   = 0;

% Main loop
while true

    % z-update
    z  = H1_adj(hard_threshold(H1_fwd(z), r))./N;
    z(mask) = k_sens(mask);

    i = i + 1;

    % calculate rmse
    if ~isempty(truth)
        rmse(i) = norm(z(:)-truth(:))/norm(truth(:));
    end    

    % print cost
    cost = norm(z(:).*mask(:) - k_sens(:))^2;
    if ~isempty(noise)
        chi2 = norm(reshape((z.*mask - k_sens)./reshape(sqrt(diag(noise)),1,1,[]),[],1))^2/nnz(mask);
    else
        chi2 = 0;
    end
    if isempty(truth)
        fprintf(1,'Iter: %03d  Cost: %8.3g  Chi^2/N: %8.3g\n', i, cost, chi2);   
    else
        fprintf(1,'Iter: %03d  Cost: %8.3g  Chi^2/N: %8.3g  RMSE: %8.3g\n', i, cost, chi2, rmse(i));   
    end
    
    % loop condition
    if niters
        if i == niters
            disp('Max Iterations Reached');
            break;
        end
    else
        if chi2 > 1
            disp('Chi^2 Criteria Reached');
            break;
        elseif i == 500
            disp('Max Iterations Reached');
            break;
        end
    end
end
if ~isempty(truth)
    varargout{1} = rmse;
end
end
