function [ y ] = FFT3D_mask( x,mask_channel )
%FFT3D Summary of this function goes here
%   Detailed explanation goes here

y = fftshift(fftn(ifftshift(x))).*mask_channel;

end

