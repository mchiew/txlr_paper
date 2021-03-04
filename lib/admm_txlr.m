function [z, varargout] = admm_txlr(k_sens, kernel, niters, w, opts)
% z = admm_txlr(k_sens, kernel, niters, w, [opts])
%
%   Required Inputs: 
%               k_sens  -   under-sampled [kx,ky,Rx,Tx] k-space
%               kernel  -   kernel size [Nx, Ny]
%               niters  -   number of ADMM iterations (suggested: 50)
%                           if 0, iterates until Chi^2/N = 1
%               w       -   [1x2] weights for nuclear norm terms [a, b] or
%                           rank constraints [r1, r2]
%                           if any term set to zero, ignores constraint
%                           e.g. [50, 0] ignores the H2 constraint
%                           and [0, 40] ignores the H1 constraint
%                           whereas [50,50] uses both
%   Optional Inputs: 
%           opts.noise  -   receiver noise covariance matrix
%                           defaults to 1E-6*eye(Rx), only used if niters = 0
%                           which uses a Chi^2/N = 1 criteria for stopping
%           opts.mode   -   defaults to 'hard' for hard thresholding
%                           other option is 'soft', for soft thresholding
%           opts.truth  -   ground truth data
%                           if present, varargout{1} = rmse at each iteration
%           opts.init   -   value to initialise reconstruction
%                           useful for "warm" starts
%                           defaults to zeros()
%
%   Output(s):
%                   z   -   recovered k-space
%           varargout{1}-   if opts.truth provided, vector of rmse vs iteration
%
% When opt.mode == 'hard':
% ADMM implementation of txlr which solves the following cost function:
% min_{z} = 0.5||Mz - k||_2^2 
% such that rank(F1(H(z))) = r1, and rank(F2(H(z))) = r2
%
% where:
%       z   is the multi-channel k-space, 
%       M   is the sampling operator (mask)
%       k   is the sampled k-space
%       F1  is the first tensor unfolding operator
%           this could correspond to a [*Rx, *Tx] unfolding of the Hankel
%           tensor
%       r1  is the rank constraint for the first unfolding
%       H   is the operator that produces Hankel matrices from k-spaces for
%           a given kernel size
%       F2  is the second tensor unfolding operator
%           this could correspond to a [*Tx, *Rx] unfolding 
%       r2  is the rank constraint for the second unfolding
%
%
% When opt.mode == 'soft':
% ADMM implementation of txlr which solves the following cost function:
% min_{z} = 0.5||Mz - k||_2^2 + a||F1(H(z))||_* + b||F2(H(z))||_*
% 
% where:
%       z   is the multi-channel k-space, 
%       M   is the sampling operator (mask)
%       k   is the sampled k-space
%       a   is the weight for the first nuclear norm term
%       F1  is the first tensor unfolding operator
%           this could correspond to a [*Rx, *Tx] unfolding of the Hankel
%           tensor
%       H   is the operator that produces Hankel matrices from k-spaces for
%           a given kernel size
%       b   is the weight for the second nuclear norm term
%       F2  is the second tensor unfolding operator
%           this could correspond to a [*Tx, *Rx] unfolding 
%
%
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
% v1.0 (04/04/20)   -   Initial version
% v1.1 (05/04/20)   -   Added hard thresholding option
% v1.2 (13/04/20)   -   If either constraint is 0, constraint is ignored
%                       Added chi^2 stopping criteria if noise cov provided
% v1.3 (23/04/20)   -   Swapped noise/mode parameter order
% v1.4 (01/05/20)   -   Added init for warm starts
% v1.5 (02/05/20)   -   Moved optional arguments to struct
%                       Merged all implementations (admm_vc, admm_primo) into single function

% set default optional parameters
if nargin < 5
    opts    =   struct();
end
if ~isfield(opts, 'noise')
    opts.noise  =   [];
end
if ~isfield(opts, 'mode')
    opts.mode   =   'hard';
end
if ~isfield(opts, 'truth')
    opts.truth  =   [];
end
if ~isfield(opts, 'init')
    opts.init   =   [];
end
    
% internal admm params
p   =   1E-6;       % rho (penalty param)
m   =   1.1;        % varying penalty
r   =   1.5;        % relaxation parameter

% sampling mask
mask    =   k_sens ~= 0;

% setup helper functions
kdims   =   size(k_sens); 
hdims   =   size(Hankel_fwd(k_sens, kernel, kdims)); 
[~,N]   =   Hankel_adj(zeros(hdims), kernel, kdims);

H1_fwd  =   @(x) fold_rx(Hankel_fwd(x, kernel, kdims));
H1_adj  =   @(x) Hankel_adj(unfold_rx(x, hdims), kernel, kdims);
H2_fwd  =   @(x) fold_tx(Hankel_fwd(x, kernel, kdims));
H2_adj  =   @(x) Hankel_adj(unfold_tx(x, hdims), kernel, kdims);

% initalise variables
if ~isempty(opts.init)
    z   =   opts.init;
    p   =   1E-4;
else
    z   =   k_sens*0;
end
U1  =   H1_fwd(z);
U2  =   H2_fwd(z);
z0  =   z;
g1  =   U1*0;
g2  =   U2*0;
i   =   0;
chi2=   0;


% ADMM main loop
while true
    % x-updates
    switch opts.mode
        case 'soft'
            X1  =   soft_threshold(H1_fwd(z) - U1, w(1)/p);
            X2  =   soft_threshold(H2_fwd(z) - U2, w(2)/p);
        case 'hard'
            X1  =   hard_threshold(H1_fwd(z) - U1, w(1));
            X2  =   hard_threshold(H2_fwd(z) - U2, w(2));
    end
    
    % z-update
    z   =   (spdiags(mask(:),0,numel(mask),numel(mask))+ nnz(w)*p*spdiags(N(:),0,numel(N),numel(N)))\reshape(k_sens + (w(1)>0)*p*H1_adj(r*X1+(1-r)*g1+U1) + (w(2)>0)*p*H2_adj(r*X2+(1-r)*g2+U2),[],1);
    z   =   reshape(z, kdims);    
    
    % dual update
    h1  =   H1_fwd(z);
    h2  =   H2_fwd(z);
    U1  =   U1 + r*X1 + (1-r)*g1 - h1;
    U2  =   U2 + r*X2 + (1-r)*g2 - h2;
    
    % penalty adjustment
    s1  =   p*H1_fwd(z-z0);
    s2  =   p*H2_fwd(z-z0);

    a   =   norm([X1(:)-h1(:);X2(:)-h2(:)]);
    b   =   norm([s1(:);s2(:)]);
    if a > 10*b
        p   =   p*m;
        U1  =   U1/m;
        U2  =   U2/m;
    elseif b > 10*a
        p   =   p/m;
        U1  =   U1*m;
        U2  =   U2*m;
    end
    z0  =   z;
    g1  =   h1;
    g2  =   h2;    
    i   =   i+1;
    
    % calculate rmse
    if ~isempty(opts.truth)
        rmse(i) =   norm(z(:)-opts.truth(:))/norm(opts.truth(:));
    end
    
    % print cost
    cost    =   norm(z(:).*mask(:) - k_sens(:))^2;
    if ~isempty(opts.noise)
        chi2    =   norm(reshape((z.*mask - k_sens)./reshape(sqrt(diag(opts.noise)),1,1,[]),[],1))^2/nnz(mask);
    end
    if ~isempty(opts.truth)
        if ~isempty(opts.noise)
            fprintf(1,'Iter: %03d  Cost: %8.3g  Chi^2/N: %8.3g  RMSE: %8.3g\n', i, cost, chi2, rmse(i));   
        else
            fprintf(1,'Iter: %03d  Cost: %8.3g  RMSE: %8.3g\n', i, cost, rmse(i));   
        end
    else
        if ~isempty(opts.noise)
            fprintf(1,'Iter: %03d  Cost: %8.3g  Chi^2/N: %8.3g\n', i, cost, chi2);   
        else
            fprintf(1,'Iter: %03d  Cost: %8.3g\n', i, cost);   
        end
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
if ~isempty(opts.truth)
    varargout{1}    =   rmse;
end
end
