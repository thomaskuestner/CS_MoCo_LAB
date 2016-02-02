function y=ifft_mid0(x)
% ifft_mid0: Inverse FFT with grid midpoint at 0
% Usage:
%   X = ifft_mid0(Y)
% Inputs:
%   Y	Array(n) 
% Outputs:
%   X   Array(n)
% Description:
%  Performs 1d ifft with grid (-n/2):(n/2-1) on both time
%  and frequency side. 
%    y(k) = sum_{t=-n/2}^{n/2-1} exp(i 2pi/n kt) x(t) , (-n/2) <= k < n/2
%
	y=fftshift(ifft(fftshift(x)));

%
%	Copyright (c) 1998 Amir Averbuch
%

%	Version Information
%		V1.0	12/28/00	DLD/ Documentation
