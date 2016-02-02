function [res] = cgL1ESPIRiT(obj, kData, x0, FT, MapOp, nIterCG, XOP, lambda, alpha,nIterSplit)
%
%[res] = cgESPIRiT(kData, x0, FT, MapOp, nIterCG, [ XOP, lambda, alpha,nIterSplit)
%
% Implementation of image-domain L1-Wavelet regularized ESPIRiT reconstruction from arbitrary
% k-space. The splitting is based on F. Huang MRM 2010;64:1078?1088
% Solves the problem: || Em - y ||^2 + \lambda ||\Psi x||_1 + \alpha||x-m||^2 
% by splitting into two subproblems:
%   I) ||Em - y || ^2 + \alpha ||x-m||^2
%   II) \alpha||x-m||^2 + \lambda ||\Psi x||_1
%
% large \alpha has slower splitting iterations and faster cg ones. 
% 
% 
% 
% Inputs: 
%       kData   - k-space data matrix it is 3D corresponding to [readout,interleaves,coils] 
%       x0      - Initial estimate of the coil images
%       FT - fft/nufft operator (see @NUFFT @p2DFT classes)
%       nIterCG   - number of LSQR iterations
%       MapOp     - ESPIRiT Operator (See @ESPIRiT)
%       XOP - transform thresholding operator for l1 based soft
%             thresholding
%       lambda   - L1-Wavelet penalty
%       alpha    - splitting parameter (0.5 default)
%       nIterSplit - number of splitting iterations
%       
% Outputs:
%       res - reconstructed ESPIRiT images
%       
%
% 
%
% (c) Michael Lustig 2006, modified 2010, 2013

if nargin < 6
    XOP = 1;
    lambda = 0;
    alpha = 0;
    nIterSplit = 1;
end

    
N = prod(size(x0));
M = prod(size(kData));
imSize = size(x0);

% make dyadic size if Wavelets are used. 
if strcmp(class(XOP),'Wavelet') == 1

    if length(imSize)>2
        imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2)))),imSize(3)];
    else
        imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2))))];
    end
else
    imSize_dyd = imSize;
end

    
dataSize = [size(kData)];

res = x0(:);

dispProgress('CG',0,nIterSplit);
for n=1:nIterSplit;
    b = [kData(:); sqrt(alpha)*res(:)];

    [res,FLAG,RELRES,ITER,RESVEC,LSVEC] = lsqr(@(x,tflag)afun(x,tflag,FT,MapOp,dataSize, imSize, alpha), b, [], nIterCG,speye(N,N),speye(N,N), res(:));
    res = reshape(res,imSize);
 
    if lambda > 0
        tmp = zpad(res,imSize_dyd);
        tmp = XOP*tmp;
        tmp = SoftThresh(tmp,lambda/sqrt(alpha));
        res = XOP'*tmp;
        res = reshape(res,imSize_dyd);
        res = crop(res,imSize);
    end
    
%     obj1 = (FT* (MapOp * res) - kData); 
    
%     figure(100), imshow3(abs(res),[]), drawnow;
%     disp(sprintf('Iteration: %d, consistency: %f',n,norm(obj1(:))));
    dispProgress('CG',n/nIterSplit);
                                                                
end
dispProgress('CG','Close');

res = weight(MapOp,res);



function [y, tflag] = afun(x, tflag, FT, MapOp, dataSize, imSize,  alpha)

    
    if strcmp(tflag,'transp')
        
        
        y = reshape(x(1:prod(dataSize)),dataSize);
        xtmp = x(prod(dataSize)+1:end);
        
        x = FT'.*y;
        x = MapOp'*x;
        y = x(:)+ sqrt(alpha)*xtmp(:);
    else
        
        x = reshape(x,imSize);
        x_ = MapOp * x;
        y = FT.*x_;
        y = [y(:); sqrt(alpha) * x(:)];
    end

