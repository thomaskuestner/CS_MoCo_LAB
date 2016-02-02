function  res = p2DFT(mask,imSize,ph,mode)

%res = p2DFT(mask,imSize [ ,phase,mode])
%
%
%	Implementation of partial Fourier operator.
%	
%	input:
%			mask - 2D matrix with 1 in entries to compute the FT and 0 in ones tha
%				are not computed.
%			imSize - the image size (1x2)
%			phase - Phase of the image for phase correction
%			mode - 1- real, 2-cmplx
%
%	Output:
%			The operator
%
%	(c) Michael Lustig 2007

if nargin <3
	ph = 1;
end
if nargin <4
	mode = 2; % 0 - positive, 1- real, 3-cmplx
end


res.adjoint = 0;
res.mask = mask;
res.imSize = imSize;


res.dataSize = [size(mask)];
res.ph = ph;
res.mode = mode;
res = class(res,'p2DFT');

