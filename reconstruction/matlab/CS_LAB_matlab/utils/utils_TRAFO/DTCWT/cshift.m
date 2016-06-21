function y = cshift(x, m)

% Circular Shift
%
% USAGE:
%    y = cshift(x, m)
% INPUT:
%    x - N-point vector
%    m - amount of shift
% OUTPUT:
%    y - vector x will be shifed by m samples to the left
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

N = length(x);
n = 0:N-1;
n = mod(n-m, N);
y = x(n+1);

