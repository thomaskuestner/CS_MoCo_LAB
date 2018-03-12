function [ y ] = sparse_trans( x,waveletStages,waveletFilterName )
%FFT3D Summary of this function goes here
%   Detailed explanation goes here

% atm 
% x: real -> y: real
% x: imag -> y: real % for backwards compatibility with certain functions.

if isreal(x{1,1})
    y = cell(1,size(x,2));
    for j = 1:size(x,2)
        y{1,j} = wavedec2(x{1,j},waveletStages,waveletFilterName);
    end;
else
    y = cell(1,2*size(x,2));
    for j = 1:size(x,2)
        y{1,j} = wavedec2(real(x{1,j}),waveletStages,waveletFilterName);
        y{1,j+size(x,2)} = wavedec2(imag(x{1,j}),waveletStages,waveletFilterName);
    end;
end;

end

