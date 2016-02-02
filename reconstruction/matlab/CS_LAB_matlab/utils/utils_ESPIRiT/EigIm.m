function [EigenVecs,EigenVals,REigenVecs] = EigIm(KERNEL,v)

imSize = [size(KERNEL,1),size(KERNEL,2)];

nc = size(KERNEL,3);
nv = size(KERNEL,4);

EigenVecs = zeros(imSize(1), imSize(2), nc, min(nc,nv));
EigenVals = zeros(imSize(1), imSize(2), min(nc,nv));
REigenVecs = zeros(imSize(1),imSize(2), max(nc,nv),nc);
for n=1:prod(imSize)
    [x,y] = ind2sub([imSize(1),imSize(2)],n);
    mtx = squeeze(KERNEL(x,y,:,:));

    %[C,D] = eig(mtx*mtx');
    [C,D,V] = svd(mtx,'econ');
    
    ph = repmat(exp(-i*angle(C(1,:))),[size(C,1),1]);
    C = v*(C.*ph);
    D = real(diag(D));
    EigenVals(x,y,:) = D(end:-1:1);
    EigenVecs(x,y,:,:) = C(:,end:-1:1);
    REigenVecs(x,y,:,:) = V(:,end:-1:1);
end

