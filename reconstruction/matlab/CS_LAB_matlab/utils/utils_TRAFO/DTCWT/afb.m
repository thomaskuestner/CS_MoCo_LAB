function [lo, hi] = afb(x, af)

% Analysis filter bank
%
% USAGE:
%    [lo, hi] = afb(x, af)
% INPUT:
%    x - N-point vector, where
%            1) N is even
%            2) N >= length(af)
%    af - analysis filters
%    af(:, 1) - lowpass filter (even length)
%    af(:, 2) - highpass filter (even length)
% OUTPUT:
%    lo - Low frequecy output
%    hi - High frequency output
% EXAMPLE:
%    [af, sf] = farras;
%    x = rand(1,64);
%    [lo, hi] = afb(x, af);
%    y = sfb(lo, hi, sf);
%    err = x - y; 
%    max(abs(err))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

N = length(x);
L = length(af)/2;
x = cshift(x,-L);

% lowpass filter
lo = upfirdn(x, af(:,1), 1, 2);
lo(1:L) = lo(N/2+[1:L]) + lo(1:L);
lo = lo(1:N/2);

% highpass filter
hi = upfirdn(x, af(:,2), 1, 2);
hi(1:L) = hi(N/2+[1:L]) + hi(1:L);
hi = hi(1:N/2);

