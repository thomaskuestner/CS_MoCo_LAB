function w = dwt(x, J, af)

% Discrete 1-D Wavelet Transform
%
% USAGE:
%    w = dwt(x, J, af)
% INPUT:
%    x - N-point vector, where
%            1) N is divisible by 2^J
%            2) N >= 2^(J-1)*length(af)
%    J - number of stages
%    af - analysis filters
%    af(:, 1) - lowpass filter (even length)
%    af(:, 2) - highpass filter (evenlength)
% OUTPUT:
%    w{j}, j = 1..J+1 - DWT coefficients
% EXAMPLE:
%    [af, sf] = farras;
%    x = rand(1,64);
%    w = dwt(x,3,af);
%    y = idwt(w,3,sf);
%    err = x - y; 
%    max(abs(err))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

for k = 1:J
    [x w{k}] = afb(x, af);
end
w{J+1} = x;

