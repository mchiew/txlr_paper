function k = crop_k(k, dims, center)

if nargin < 3
    center = false;
end

if center
    [~,ind] =   max(reshape(mean(mean(abs(k),3),4),[],1));
    [i,j]   =   ind2sub([size(k,1),size(k,2)], ind);
    k       =   k(i-dims(1)/2:i+dims(1)/2-1,:,:,:);
    k       =   k(:,j-dims(2)/2:j+dims(2)/2-1,:,:);
else
    N1  = ceil((size(k,1)+1)/2);
    k   = k(N1-dims(1)/2:N1+dims(1)/2-1,:,:,:);

    N2  = ceil((size(k,2)+1)/2);
    k   = k(:,N2-dims(2)/2:N2+dims(2)/2-1,:,:);
end
