function res = SPIRiT(kData, kernel)

%  res = SPIRiT(kData, kernel)
%
%	Apply the SPIRiT operator on the data
%
%	input:
%		kData: k-space data (sx x sy x ncoils)
%		kernel : SPIRiT kernel
%
%	(c) Michael Lustig 2006

	
	res = kData*0;

	numCoils = size(kData,3);
	for n=1:numCoils
		res(:,:,n) = filter2(kernel(:,:,n), kData(:,:,n),'same');
	end
	res = sum(res,3);

