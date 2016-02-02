function  res = SPIRiT(kernel,method, imSize)
%res = SPIRiT(kernel [,method, imSize])
%	Implementation of the SPIRiT kernel operator
%
%   
%   Constructor inputs:
%       kernel: [kx,ky,nCoils, nCoils] the spirit 2D convolution kernel.
%               See corrMatrix.m, calibrate.m and the demos on how to
%               generate a kernel.
%       method: implementation type of the operator. There are three
%               'conv'  is a k-space convolution implementation (slowest)
%                       that operates on k-space data and produces k-space data.
%               'fft'   is an fft based implementation of the k-space
%                       convolution through image domain multiplication. 
%                       It operates on k-sapce data and produces k-space data. 
%               'image' The SPIRiT operator is applied to image space data.
%                       This is useful when using image-based non-cartesian
%                       reconstruction. In this case the SPIRiT operator
%                       operates on image data and produces image data. 
%
%       imSize:  Size of the resulting image (only needed for 'fft' and
%                'image' modes
%
%
%   Example:
%
%   [x,y] = meshgrid(linspace(0,1,128));
%   % Generate fake Sensitivity maps
%   sMaps = cat(3,x.^2,1-x.^2,y.^2,1-y.^2);
%   % generate 4 coil phantom
%   imgs = repmat(phantom(128),[1,1,4]).*sMaps;
%   DATA = fft2c(imgs);
%   % crop 20x20 window from the center of k-space for calibration
%   kCalib = crop(DATA,[20,20,4]);
%
%   %calibrate a kernel
%   kSize = [5,5];
%   coils = 4;
%   kernel = zeros([kSize,coils,coils]);
%   [AtA,] = corrMatrix(kCalib,kSize);
%   for n=1:coils
%       kernel(:,:,:,n) = calibrate(AtA,kSize,coils,n,0.01);
%   end
%   GOP = SPIRiT(kernel, 'fft',[128,128]);
%
%   % undersample by a factor of 2
%   DATA(1:2:end,2:2:end,:) = 0;
%   DATA(2:2:end,1:2:end,:) = 0;
%   
%   %reconstruct:
%   [res] = cgSPIRiT(DATA,GOP, 20, 1e-5, DATA);
%   figure, imshow(cat(2,sos(imgs), 2*sos(ifft2c(DATA)), sos(ifft2c(res))),[]);
%   title('full,  zero-fill,   result')
%   
%   See Also:
%               calibrate.m, corrMatrix.m
%
% (c) Michael Lustig 2007
     

    if nargin < 2
	    method = 'conv';
	    KERNEL = [];
    end
    
    if strcmp(method,'conv')==1
	    KERNEL = [];
    end


    if strcmp(method,'fft')==1 & nargin < 3
	    error('must provide image size');
    end

    % for methods 'fft' and 'image' precompute the image domain kernel by
    % zero-padding and inverse fft
    
    if strcmp(method,'fft')==1 | strcmp(method,'image')==1
	    for n=1:size(kernel,4)
            KERNEL(:,:,:,n) = ifftnshift(zpad(kernel(end:-1:1,end:-1:1,:,n), imSize(1), imSize(2), size(kernel,3)),1:2);
%             KERNEL(:,:,:,n) = (ifft2c(zpad(kernel(end:-1:1,end:-1:1,:,n)*sqrt(imSize(1)*imSize(2)), imSize(1), imSize(2), size(kernel,3))));
        end
    end

    

    if strcmp(method,'fft')==0 & strcmp(method,'conv')==0 & strcmp(method,'image')==0
	    eror('no such method');
    end

    
    res.kernel = kernel;
    res.adjoint = 0;
    res.KERNEL = KERNEL;
    res.method = method;
    res.imSize = imSize;
    res = class(res,'SPIRiT');

