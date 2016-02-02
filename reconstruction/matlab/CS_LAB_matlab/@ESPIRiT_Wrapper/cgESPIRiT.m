function [res,imgs, RESVEC] = cgESPIRiT(obj, y, ESP, nIter, lambda, x0)
%
%
%  res = cgESPIRiT(y,ESP, nIter, lambda,x0)
%  
%  Implementation of the Cartesian, conjugate gradient ESPIRiT
%  reconstruction. This implementation is similar to Cartesian SPIRiT. It
%  only solves for missing data in k-space. 
%
%
%  Input:
%		y -	Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled.
%		ESP -	the ESPIRiT operator obtained by calibration
%		nIter -	Maximum number of iterations
%		lambda-	Tykhonov regularization parameter
%       x0    - Initil value
%
% Outputs:
%		res - Full k-space matrix
%       imgs - Combined ESPIRiT images. 
%
%
%   
% (c) Michael Lustig 2013
%

if nargin < 4
	lambda = 0;
end

if nargin < 5
     x0 = y;
end



[sx,sy,nCoils] = size(y);

idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);
N = length(idx_nacq(:));



yy = fft2c(ESP*(ESP'*(ifft2c(y))))-y; yy = [-yy(:); idx_nacq(:)*0];

[tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,1e-6,nIter, speye(N,N),speye(N,N),x0(idx_nacq),ESP,sx,sy,nCoils,idx_nacq, lambda);

res = y;
res(idx_nacq) = tmpres;
imgs = ESP'*(ifft2c(res));


function [res,tflag] = aprod(x,ESP,sx,sy,nCoils,idx_nacq, lambda,tflag)
	
	
	if strcmp(tflag,'transp');
		tmpy = reshape(x(1:sx*sy*nCoils),sx,sy,nCoils);
        tmpy = ifft2c(tmpy);
        res = ESP*(ESP'*tmpy)-tmpy;
        res = fft2c(res);
        res = res(idx_nacq)+ x(sx*sy*nCoils+1:end)*lambda;
	
    else
		tmpx = zeros(sx,sy,nCoils);
		tmpx(idx_nacq) = x;
        tmpx = ifft2c(tmpx);
		res = ESP*(ESP'*tmpx)-tmpx;
        res = fft2c(res);
		res = [res(:) ; lambda*x(:)];
	end






