function y = sfb(lo, hi, sf)

% Synthesis filter bank
%
% USAGE:
%    y = sfb(lo, hi, sf)
% INPUT:
%    lo - low frqeuency input
%    hi - high frequency input
%    sf - synthesis filters
%    sf(:, 1) - lowpass filter (even length)
%    sf(:, 2) - highpass filter (even length)
% OUTPUT:
%    y - output signal
% See also afb
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

N = 2*length(lo);
L = length(sf);
lo = upfirdn(lo, sf(:,1), 2, 1);
hi = upfirdn(hi, sf(:,2), 2, 1);
y = lo + hi;
y(1:L-2) = y(1:L-2) + y(N+[1:L-2]);
y = y(1:N);
y = cshift(y, 1-L/2);

