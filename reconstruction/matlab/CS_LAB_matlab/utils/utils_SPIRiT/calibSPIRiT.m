function kernel = calibSPIRiT(kCalib, kSize, nCoils, CalibTyk)

% kernel = calibSPIRiT(kCalib, kSize, nCoils, CalibTyk)
%
% Function calibrates a SPIRiT kernel from a calibration area in k-space
%
%
% (c) Michael Lustig 2013
%

[AtA] = dat2AtA(kCalib,kSize);
for n=1:nCoils
	kernel(:,:,:,n) = calibrate(AtA,kSize,nCoils,n,CalibTyk);
end

