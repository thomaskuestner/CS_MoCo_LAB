function res = adjSPIRiT(kData, kernel)

%  res = adjSPIRiT(kData, kernel)
%
%	backprojection of the SPIRiT operator on the data
%
%	input:
%		kData: k-space data (sx x sy x ncoils)
%		kernel : SPIRiT kernel
%
%	(c) Michael Lustig 2006

	kernel = conj(kernel(end:-1:1,end:-1:1,:));
	res = kData*0;

	numCoils = size(kernel,3);
	for n=1:numCoils
		res(:,:,n) = filter2(kernel(:,:,n), kData(:,:),'same');
	end
	

