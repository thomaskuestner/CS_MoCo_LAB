function w = dwt2D(x, J, af)

% discrete 2-D wavelet transform
%
% USAGE:
%   w = dwt2D(x, stages, af)
% INPUT:
%   x - N by M matrix
%       1) M, N are both even
%       2) min(M,N) >= 2^(J-1)*length(af)
%   J - number of stages
%   af - analysis filters
% OUPUT:
%   w - cell array of wavelet coefficients
% EXAMPLE:
%   [af, sf] = farras;
%   x = rand(128,64);
%   J = 3;
%   w = dwt2D(x,J,af);
%   y = idwt2D(w,J,sf);
%   err = x - y; 
%   max(max(abs(err)))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

for k = 1:J
    [x w{k}] = afb2D(x, af, af);
end
w{J+1} = x;

