function Y = fft2_mid0(X)
% fft2_mid0 -- 2d fft with argument [-pi,pi] (instead of [0,2pi])
% Usage:
%   Y = fft_mid0(X)
% Inputs:
%   X	Array(n) 
% Outputs:
%   Y   Array(n)
% Description:
%  Performs 1d fft with grid (-n/2):(n/2-1) on both time
%  and frequency side. 
%
% y(k) = sum_{-n/2 <= t1,t2 < n/2} exp(-i 2pi/n (k1 t1 + k2 t2)) x(t1,t2),
%                                         (-n/2) <= k1,k2 < n/2
%

Y = fftshift(fft2(fftshift(X)));

