function y = idwt3D(w, J, sf)

% Inverse 3-D Discrete Wavelet Transform
%
% USAGE:
%   y = idwt3D(w, J, sf)
% INPUT:
%   w - wavelet coefficient
%   J  - number of stages
%   sf - synthesis filters
% OUTPUT:
%   y - output array
% See: dwt3D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

y = w{J+1};
for k = J:-1:1
   y = sfb3D(y, w{k}, sf, sf, sf);
end

