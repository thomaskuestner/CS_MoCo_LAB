function [kernel, S] = dat3Kernel(data, kSize)
% kernel = dat2Kernel(data, kSize,thresh)
%
% Function to perform k-space calibration step for ESPIRiT and create
% k-space kernels. For 3D multi-coil images.  
% 
% Inputs: 
%       data - calibration data [kx,ky,coils]
%       kSize - size of kernel (for example kSize=[5,5,5])
%
% Outputs: 
%       kernel - k-space kernels matrix (not cropped), which correspond to
%                the basis vectors of overlapping blocks in k-space
%       S      - (Optional parameter) The singular vectors of the
%                 calibration matrix
%
%
% See also:
%           kernelEig
%
% (c) Michael Lustig 2013, 
% modified by: Thomas Küstner 2013


[sx,sy,sz,nc] = size(data);
% imSize = [sx,sy] ;

tmp = im3row(data,kSize); [tsx,tsy,tsz,tscha] = size(tmp);
A = reshape(tmp,tsx,tsy*tsz*tscha);

[~,S,V] = svd(A,'econ');
    
kernel = reshape(V,kSize(1),kSize(2),kSize(3),nc,size(V,2));
S = diag(S);S = S(:);
