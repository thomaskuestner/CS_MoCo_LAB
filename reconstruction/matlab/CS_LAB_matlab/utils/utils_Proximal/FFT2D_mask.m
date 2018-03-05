function [ y ] = FFT2D_mask( x,mask_channel )
%FFT2D Summary of this function goes here
%   Detailed explanation goes here

y = fftshift(fft2(ifftshift(x))).*mask_channel;

end

