function w = dualtree3D(x, J, Faf, af)

% 3D Dual-Tree Discrete Wavelet Transform
%
% USAGE:
%   w = dualtree3D(x, J, Faf, af)
% INPUT:
%   x - 3-D array
%   J - number of stages
%   Faf - first stage filters
%   af - filters for remaining stages
% OUPUT:
%   w{j}{i}{d} - wavelet coefficients
%        j = 1..J, i = 1..4, d = 1..7
%   w{J+1}{i} - lowpass coefficients
%        i = 1..4
% EXAMPLE:
%   x = rand(64,64,64);
%   J = 3;
%   [Faf, Fsf] = FSfarras;
%   [af, sf] = dualfilt1;
%   w = dualtree3D(x, J, Faf, af);
%   y = idualtree3D(w, J, Fsf, sf);
%   err = x - y;
%   max(max(max(abs(err))))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/


% normalization
x = x/2;

M = [
    1 1 1
    2 2 1
    2 1 2
    1 2 2
];

for i = 1:4
    f1 = M(i,1);
    f2 = M(i,2);
    f3 = M(i,3);
    [xi w{1}{i}] = afb3D(x, Faf{f1}, Faf{f2}, Faf{f3});
    for k = 2:J
        [xi w{k}{i}] = afb3D(xi, af{f1}, af{f2}, af{f3});
    end
    w{J+1}{i} = xi;
end

for k = 1:J
    for m = 1:7
        [w{k}{1}{m} w{k}{2}{m} w{k}{3}{m} w{k}{4}{m}] = ...
            pm4(w{k}{1}{m}, w{k}{2}{m}, w{k}{3}{m}, w{k}{4}{m});
    end
end


