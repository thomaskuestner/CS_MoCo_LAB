function [KERNEL,v] = kernel2Im(kernel,imSize)

nc = size(kernel,3);
nv = size(kernel,4);
kSize = [size(kernel,1), size(kernel,2)];

% "rotate kernel to order by maximum variance"
k = permute(kernel,[1,2,4,3]);, k =reshape(k,prod([kSize,nv]),nc);

if size(k,1) < size(k,2)
    [u,s,v] = svd(k);
else
    
    [u,s,v] = svd(k,'econ');
end

k = k*v;
kernel = reshape(k,[kSize,nv,nc]); kernel = permute(kernel,[1,2,4,3]);


KERNEL = zeros(imSize(1), imSize(2), size(kernel,3), size(kernel,4));
for n=1:size(kernel,4)
    KERNEL(:,:,:,n) = (fft2c(zpad(conj(kernel(end:-1:1,end:-1:1,:,n))*sqrt(imSize(1)*imSize(2)), ...
        [imSize(1), imSize(2), size(kernel,3)])));
end
KERNEL = KERNEL/sqrt(prod(kSize));

