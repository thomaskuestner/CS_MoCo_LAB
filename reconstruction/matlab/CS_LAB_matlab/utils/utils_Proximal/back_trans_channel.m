function [ y ] = back_trans( x,waveS,waveletFilterName )
%FFT3D Summary of this function goes here
%   Detailed explanation goes here

if isreal(x{1,1})
    y = cell(1,size(x,2));
    for j = 1:size(x,2)
        y{1,j} = waverec2(x{1,j},waveS,waveletFilterName);
    end;
else
    y = cell(1,2*size(x,2));
    for j = 1:size(x,2)
        y{1,j} = waverec2(real(x{1,j}),waveS,waveletFilterName);
        y{1,j+size(x,2)} = waverec2(imag(x{1,j}),waveS,waveletFilterName);
    end;
end;
end

