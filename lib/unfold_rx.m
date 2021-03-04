% ADMM TxLR Helper function
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
function x = unfold_rx(x, dims)
    if numel(dims) == 3
        dims(4) = 1;
    end
    x = permute(reshape(x, dims([1, 3, 2, 4])), [1, 3, 2, 4]);
end
