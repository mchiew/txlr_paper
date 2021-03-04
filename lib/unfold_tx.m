% ADMM TxLR Helper function
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
function x = unfold_tx(x, dims)
    if numel(dims) == 3
        dims(4) =   1;
    end
    x = permute(reshape(x, dims([1, 4, 2, 3])), [1, 3, 4, 2]);
end
