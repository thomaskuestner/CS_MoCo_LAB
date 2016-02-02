function [kernel,rawkernel] = calibrate(obj, AtA, nCoil, coil)

sampling = ones([obj.kernelSize(1:2),nCoil]);

dummyK = zeros(obj.kernelSize(1),obj.kernelSize(2),nCoil); dummyK((end+1)/2,(end+1)/2,coil) = 1;
idxY = find(dummyK);
sampling(idxY) = 0;
idxA = find(sampling);

Aty = AtA(:,idxY); Aty = Aty(idxA); % correlation values to target point, take complete neighbourhood and not just aquired ones
AtA = AtA(idxA,:); AtA =  AtA(:,idxA); % kick out the searched point

kernel = sampling*0;

lambda = norm(AtA,'fro')/size(AtA,1)*obj.calibTyk;

rawkernel = inv(AtA + eye(size(AtA))*lambda)*Aty; % grappa weighting values
kernel(idxA) = rawkernel; 

end













