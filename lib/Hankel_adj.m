% ADMM TxLR Helper function
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
function [x,N]= Hankel_adj(h, kernel, dims)
    if numel(dims) == 3
        dims(4) = 1;
    end
    Nx  = dims(1);
    Ny  = dims(2);
    N1  = dims(3);
    N2  = dims(4);

    x   = zeros([Nx,Ny,N1,N2]);
    N   = zeros([Nx,Ny,N1,N2]);

    idx = 0;
    for kx = 1:Nx - kernel(1) + 1
        for ky = 1:Ny - kernel(2) + 1
            idx = idx + 1;
            x(kx:kx+kernel(1)-1, ky:ky+kernel(2)-1,:,:) = x(kx:kx+kernel(1)-1, ky:ky+kernel(2)-1,:,:) + reshape(h(:,idx,:,:),kernel(1),kernel(2),N1,N2);
            N(kx:kx+kernel(1)-1, ky:ky+kernel(2)-1,:,:) = N(kx:kx+kernel(1)-1, ky:ky+kernel(2)-1,:,:) + 1;
        end
    end
end
