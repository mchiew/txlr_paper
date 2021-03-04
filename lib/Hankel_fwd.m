% ADMM TxLR Helper function
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
function h = Hankel_fwd(x, kernel, dims)
    if numel(dims) == 3
        dims(4) = 1;
    end
    Nx  = dims(1);
    Ny  = dims(2);
    N1  = dims(3);   
    N2  = dims(4);
        
    h   = zeros(prod(kernel), prod([Nx,Ny]-kernel+1), N1, N2);

    idx = 0;
    for kx = 1:Nx - kernel(1) + 1
        for ky = 1:Ny - kernel(2) + 1
            idx = idx + 1;
            h(:, idx, :, :) = reshape(x(kx:kx+kernel(1)-1, ky:ky+kernel(2)-1,:,:),prod(kernel),1,N1,N2);
        end
    end
end
