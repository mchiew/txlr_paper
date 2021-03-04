% ADMM TxLR Helper function
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
function x = fold_tx(x, dims)
    if nargin < 2
        dims(1) = size(x,1)*size(x,4);
        dims(2) = size(x,2)*size(x,3);
    end
    x = reshape(permute(x, [1, 4, 2, 3]), dims);;
end
