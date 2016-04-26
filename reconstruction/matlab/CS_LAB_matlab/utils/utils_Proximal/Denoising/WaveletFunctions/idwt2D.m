function y = idwt2D(w, J, sf)

% Inverse 2-D Discrete Wavelet Transform
%
% USAGE:
%   y = idwt(w, J, sf)
% INPUT:
%   w - wavelet coefficients
%   J  - number of stages
%   sf - synthesis filters
% OUTPUT:
%   y - output array
% See dwt2D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

y = w{J+1};
for k = J:-1:1
   y = sfb2D(y, w{k}, sf, sf);
end

