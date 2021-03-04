% ADMM TxLR Helper function
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
function x = hard_threshold(x, r)
    if r
        [u,s,v] = svd(x,'econ');
        x       = u(:,1:r)*s(1:r,1:r)*v(:,1:r)';
    end
end
