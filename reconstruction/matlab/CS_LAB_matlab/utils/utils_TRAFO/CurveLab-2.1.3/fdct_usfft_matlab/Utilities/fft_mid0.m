function y=fft_mid0(x)
% fftmid0 -- 1d fft with argument [-pi,pi] (instead of [0,2pi])
% Usage:
%   Y = fft_mid0(X)
% Inputs:
%   X	Array(n) 
% Outputs:
%   Y   Array(n)
% Description:
%  Performs 1d fft with grid (-n/2):(n/2-1) on both time
%  and frequency side. 
%    y(k) = sum_{t=-n/2}^{n/2-1} exp(i 2pi/n kt) x(t) , (-n/2) <= k < n/2
%

y = fftshift(fft(fftshift(x)));

%
%	Copyright (c) 1998 Amir Averbuch
%

%	Version Information
%		V1.0	12/28/00	DLD/ Documentation
