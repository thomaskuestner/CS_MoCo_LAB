function y = cshift2D(x, m)

% 2D Circular Shift
% 
% USAGE:
%    y = cshift2D(x, m)
% INPUT:
%    x - M by N array
%    m - amount of shift
% OUTPUT:
%    y - matrix x will be shifed by m samples down
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

[N, M] = size(x);
n = 0:N-1;
n = mod(n-m, N);
y = x(n+1,:);

