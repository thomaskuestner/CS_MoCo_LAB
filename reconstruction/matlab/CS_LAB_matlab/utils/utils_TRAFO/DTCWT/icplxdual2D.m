function y = icplxdual2D(w, J, Fsf, sf)

% Inverse Dual-Tree Complex 2D Discrete Wavelet Transform
% 
% USAGE:
%   y = icplxdual2D(w, J, Fsf, sf)
% INPUT:
%   w - wavelet coefficients
%   J - number of stages
%   Fsf - synthesis filters for final stage
%   sf - synthesis filters for preceeding stages
% OUTPUT:
%   y - output array
% See cplxdual2D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

for j = 1:J
    for m = 1:3
        [w{j}{1}{1}{m} w{j}{2}{2}{m}] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
        [w{j}{1}{2}{m} w{j}{2}{1}{m}] = pm(w{j}{1}{2}{m},w{j}{2}{1}{m});
    end
end

y = zeros(size(w{1}{1}{1}{1})*2);
for m = 1:2
    for n = 1:2
        lo = w{J+1}{m}{n};
        for j = J:-1:2
            lo = sfb2D(lo, w{j}{m}{n}, sf{m}, sf{n});
        end
        lo = sfb2D(lo, w{1}{m}{n}, Fsf{m}, Fsf{n});
        y = y + lo;
    end
end

% normalization
y = y/2;

