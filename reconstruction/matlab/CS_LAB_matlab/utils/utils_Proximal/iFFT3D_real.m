function [ y ] = iFFT3D_real( x )
%FFT3D Summary of this function goes here
%   Detailed explanation goes here

y = real((fftshift(ifftn(ifftshift(x)))));

end

