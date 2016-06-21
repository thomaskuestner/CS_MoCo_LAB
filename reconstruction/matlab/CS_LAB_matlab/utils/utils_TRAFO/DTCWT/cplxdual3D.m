function w = cplxdual3D(x, J, Faf, af)

% 3D Complex Dual-Tree Discrete Wavelet Transform
%
% USAGE:
%   w = cplxdual3D(x, J, Faf, af)
% INPUT:
%   x - 3D array
%   J - number of stages
%   Faf - first stage filters
%   af - filters for remaining stages
% OUPUT:
%   w{j}{m}{n}{p}{d} - wavelet coefficients
%       j = 1..J, m = 1..2, n = 1..2, p = 1..2, d = 1..7
%   w{J+1}{m}{n}{d} - lowpass coefficients
%       m = 1..2, n = 1..2, p = 1..2, d = 1..7
% EXAMPLE:
%   x = rand(64,64,64);
%   J = 3;
%   [Faf, Fsf] = FSfarras;
%   [af, sf] = dualfilt1;
%   w = cplxdual3D(x, J, Faf, af);
%   y = icplxdual3D(w, J, Fsf, sf);
%   err = x - y;
%   max(max(max(abs(err))))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% normalization
x = x/sqrt(8);

for m = 1:2
    for n = 1:2
        for p = 1:2
            [lo w{1}{m}{n}{p}] = afb3D(x, Faf{m}, Faf{n}, Faf{p});
            for j = 2:J
                [lo w{j}{m}{n}{p}] = afb3D(lo, af{m}, af{n}, af{p});
            end
            w{J+1}{m}{n}{p} = lo;
        end
    end
end

for j = 1:J
    for m = 1:7
        [w{j}{1}{1}{1}{m} w{j}{2}{2}{1}{m} w{j}{2}{1}{2}{m} w{j}{1}{2}{2}{m}] = ...
            pm4(w{j}{1}{1}{1}{m}, w{j}{2}{2}{1}{m}, w{j}{2}{1}{2}{m}, w{j}{1}{2}{2}{m});
         [w{j}{2}{2}{2}{m} w{j}{1}{1}{2}{m} w{j}{1}{2}{1}{m} w{j}{2}{1}{1}{m}] = ...
            pm4(w{j}{2}{2}{2}{m}, w{j}{1}{1}{2}{m}, w{j}{1}{2}{1}{m}, w{j}{2}{1}{1}{m});
    end
end



