function [kernel,rawkernel] = calibrate(AtA, kSize, nCoil, coil, lambda, sampling)

if nargin < 6
	sampling = ones([kSize,nCoil]);
end


dummyK = zeros(kSize(1),kSize(2),nCoil); dummyK((end+1)/2,(end+1)/2,coil) = 1;
idxY = find(dummyK);
sampling(idxY) = 0;
idxA = find(sampling);

Aty = AtA(:,idxY); Aty = Aty(idxA);
AtA = AtA(idxA,:); AtA =  AtA(:,idxA);

kernel = sampling*0;

lambda = norm(AtA,'fro')/size(AtA,1)*lambda;

rawkernel = inv(AtA + eye(size(AtA))*lambda)*Aty;
kernel(idxA) = rawkernel; 














