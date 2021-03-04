% ADMM TxLR Helper function
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
function x = soft_threshold(x, t)
    if t
        [u,s,v] = svd(x,'econ');
        s2      = diag(max(diag(s)-t,0));
        x       = u*s2*v';
    end
end
