% ADMM TxLR Helper function
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
function x = fold_rx(x, dims)
    if nargin < 2
        dims(1) = size(x,1)*size(x,3);
        dims(2) = size(x,2)*size(x,4);
    end
    x = reshape(permute(x, [1, 3, 2, 4]), dims);
end
