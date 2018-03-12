function [ y ] = iFFT3D_abs( x )
%FFT3D Summary of this function goes here
%   Detailed explanation goes here

y = abs((fftshift(ifftn(ifftshift(x)))));

end

