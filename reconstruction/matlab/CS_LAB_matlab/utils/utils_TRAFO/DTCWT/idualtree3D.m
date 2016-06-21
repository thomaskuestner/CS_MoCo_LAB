function y = idualtree3D(w, J, Fsf, sf)

% Inverse 3D Dual-Tree Discrete Wavelet Transform
%
% USAGE:
%   y = idualtree3D(w, J, Fsf, sf)
% INPUT:
%   w - wavelet coefficients
%   J - number of stages
%   Fsf - synthesis filter for the last stage
%   sf - synthesis filters for the preceeding stages
% OUTPUT:
%   y - output arry
% See dualtree3D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

for k = 1:J
    for m = 1:7
        [w{k}{1}{m} w{k}{2}{m} w{k}{3}{m} w{k}{4}{m}] = ...
            pm4inv(w{k}{1}{m}, w{k}{2}{m}, w{k}{3}{m}, w{k}{4}{m});
    end
end

M = [
    1 1 1
    2 2 1
    2 1 2
    1 2 2
];

% initialize output array
y = zeros(2^J * size(w{J}{1}{1}));

for i = 1:4
    f1 = M(i,1);
    f2 = M(i,2);
    f3 = M(i,3);
    yi = w{J+1}{i};
    for k = J:-1:2
        yi = sfb3D(yi, w{k}{i}, sf{f1}, sf{f2}, sf{f3});
    end
    yi = sfb3D(yi, w{1}{i}, Fsf{f1}, Fsf{f2}, Fsf{f3});
    y = y + yi;
end

% normalization
y = y/2;


