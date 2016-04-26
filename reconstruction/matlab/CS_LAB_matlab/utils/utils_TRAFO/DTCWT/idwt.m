function y = idwt(w, J, sf)

% Inverse Discrete 1-D Wavelet Transform
%
% USAGE:
%    y = idwt(w, J, sf)
% INPUT:
%    w - wavelet coefficients
%    J - number of stages
%    sf - synthesis filters
% OUTPUT:
%    y - output signal
% See also dwt
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

y = w{J+1};
for k = J:-1:1
   y = sfb(y, w{k}, sf);
end

