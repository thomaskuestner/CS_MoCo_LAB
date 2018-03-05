function [ y ] = FFT3D( x )
%FFT3D Summary of this function goes here
%   Detailed explanation goes here

y = fftshift(fftn(ifftshift(x)));

end

