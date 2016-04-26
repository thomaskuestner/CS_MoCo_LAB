function y = icplxdual3D(w, J, Fsf, sf)

% Inverse 3D Complex Dual-Tree Discrete Wavelet Transform
%
% USAGE:
%   y = icplxdual3D(w, J, Fsf, sf)
% INPUT:
%   J - number of stages
%   Fsf - synthesis filter for last stage
%   sf - synthesis filters for preceeding stages
% OUTPUT:
%   y - output array
% See cplxdual3D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

for j = 1:J
    for m = 1:7
        [w{j}{1}{1}{1}{m} w{j}{2}{2}{1}{m} w{j}{2}{1}{2}{m} w{j}{1}{2}{2}{m}] = ...
            pm4inv(w{j}{1}{1}{1}{m}, w{j}{2}{2}{1}{m}, w{j}{2}{1}{2}{m}, w{j}{1}{2}{2}{m});
         [w{j}{2}{2}{2}{m} w{j}{1}{1}{2}{m} w{j}{1}{2}{1}{m} w{j}{2}{1}{1}{m}] = ...
            pm4inv(w{j}{2}{2}{2}{m}, w{j}{1}{1}{2}{m}, w{j}{1}{2}{1}{m}, w{j}{2}{1}{1}{m});
    end
end

y = zeros(2*size(w{1}{1}{1}{1}{1}));
for m = 1:2
    for n = 1:2
        for p = 1:2
            lo = w{J+1}{m}{n}{p};
            for j = J:-1:2
                lo = sfb3D(lo, w{j}{m}{n}{p}, sf{m}, sf{n}, sf{p});
            end
            lo = sfb3D(lo, w{1}{m}{n}{p}, Fsf{m}, Fsf{n}, Fsf{p});
            y = y + lo;
        end
    end
end

% normalization
y = y/sqrt(8);


