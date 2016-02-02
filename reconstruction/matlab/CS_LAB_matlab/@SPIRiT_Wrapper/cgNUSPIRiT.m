function [res,FLAG,RELRES,ITER,RESVEC,LSVEC] = cgNUSPIRiT(kData, x0, NUFFTOP, GOP, nIter, lambda)

% Implementation of image-domain SPIRiT reconstruction from arbitrary
% k-space. The function is based on Jeff Fessler's nufft code and LSQR 
% 
% Inputs: 
%       kData   - k-space data matrix it is 3D corresponding to [readout,interleaves,coils] 
%       x0      - Initial estimate of the coil images
%       NUFFTOP - nufft operator (see @NUFFT class)
%       GOP     - SPIRiT Operator (See @SPIRiT)
%       nIter   - number of LSQR iterations
%       lambda  - ratio between data consistency and SPIRiT consistency (1 is recommended)
%
% Outputs:
%       res - reconstructed coil images
%       FLAG,RELRES,ITER,RESVEC,LSVEC - See LSQR documentation
%
% See demo_nuSPIRiT for a demo on how to use this function.
%
% (c) Michael Lustig 2006, modified 2010

N = prod(size(x0));
imSize = size(x0);
dataSize = [size(kData)];

b = [kData(:) ; zeros(prod(imSize),1)];
[res,FLAG,RELRES,ITER,RESVEC,LSVEC] = lsqr(@(x,tflag)afun(x,NUFFTOP,GOP,dataSize, imSize,lambda,tflag), b, [], nIter,speye(N,N),speye(N,N), x0(:));

res = reshape(res,imSize);





function [y, tflag] = afun(x,NUFFTOP,GOP,dataSize,imSize,lambda,tflag)

if strcmp(tflag,'transp')
   x1 = reshape(x(1:prod(dataSize)),dataSize);
   x2 = reshape(x(prod(dataSize)+1:end),imSize);
   y = NUFFTOP'.*x1 + lambda*(GOP'*x2); 
   y = y(:);
else
    
    x = reshape(x,imSize);
    y1 = NUFFTOP.*x;
    y2 = GOP*x;
    y = [y1(:); lambda*y2(:)];
end

