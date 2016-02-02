function  res = ESPIRiT(eigenVecs,  eigenVals,  imSize)
%res = ESPIRiT(EigenVecs, [EigenVals,  imSize)
%
%	Implementation of the ESPIRiT maps projection
%

%
% (c) Michael Lustig 2013
     

[sx,sy,nc,nm] = size(eigenVecs);

if nargin < 2
    eigenVals = 1;
else
    eigenVals = repmat(permute(eigenVals,[1,2,4,3]),[1,1,nc,1]);
end

res.eigenVecs = eigenVecs.*sqrt(eigenVals);
res.eigenVals = sqrt(eigenVals);
res.adjoint = 0;
res = class(res,'ESPIRiT');

