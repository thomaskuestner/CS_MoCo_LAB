function Y = ifft2_mid0(X)
% ifft2_mid0: Inverse 2D-FFT with grid midpoint at 0
% Usage:
%   X = ifft2_mid0(Y)
% Inputs:
%   Y	Array(n) 
% Outputs:
%   X   Array(n)
% Description:
%  Performs 2d ifft with grid (-n/2):(n/2-1) on both time
%  and frequency side. 

	Y = fftshift(ifft2(fftshift(X)));

	