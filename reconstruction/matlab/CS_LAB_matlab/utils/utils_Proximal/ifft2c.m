function res = ifft2c(x)
%
%
% res = ifft2c(x)
% 
% orthonormal centered 2D ifft
%
% (c) Michael Lustig 2005

res = sqrt(length(x(:)))*ifftshift(ifft2(fftshift(x)));

