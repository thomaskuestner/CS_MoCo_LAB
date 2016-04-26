function w = dualtree(x, J, Faf, af)

% Dual-tree Complex Discrete Wavelet Transform
%
% USAGE:
%    w = dualtree(x, J, Faf, af)
% INPUT:
%   x - N-point vector
%       1) N is divisible by 2^J
%       2) N >= 2^(J-1)*length(af)
%   J - number of stages
%   Faf - filters for the first stage 
%   af - filters for the remaining stages
% OUTPUT:
%   w - DWT coefficients
%      w{j}{1}, j = 1..J - real part 
%      w{j}{2}, j = 1..J - imaginary part 
%      w{J+1}{d} - lowpass coefficients, d = 1,2
% EXAMPLE:
%    x = rand(1, 512);
%    J = 4;
%    [Faf, Fsf] = FSfarras;
%    [af, sf] = dualfilt1;
%    w = dualtree(x, J, Faf, af);
%    y = idualtree(w, J, Fsf, sf);
%    err = x - y;
%    max(abs(err))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% normalization
x = x/sqrt(2);

% Tree 1
[x1 w{1}{1}] = afb(x, Faf{1});
for j = 2:J
    [x1 w{j}{1}] = afb(x1, af{1});
end
w{J+1}{1} = x1;

% Tree 2
[x2 w{1}{2}] = afb(x, Faf{2});
for j = 2:J
    [x2 w{j}{2}] = afb(x2, af{2});
end
w{J+1}{2} = x2;


