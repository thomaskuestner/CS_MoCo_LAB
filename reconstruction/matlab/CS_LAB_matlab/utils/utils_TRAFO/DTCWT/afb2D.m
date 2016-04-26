function [lo, hi] = afb2D(x, af1, af2)

% 2D Analysis Filter Bank
%
% USAGE:
%   [lo, hi] = afb2D(x, af1, af2);
% INPUT:
%   x - N by M matrix
%       1) M, N are both even
%       2) M >= 2*length(af1)
%       3) N >= 2*length(af2)
%   af1 - analysis filters for columns
%   af2 - analysis filters for rows
% OUTPUT:
%    lo - lowpass subband
%    hi{1} - 'lohi' subband
%    hi{2} - 'hilo' subband
%    hi{3} - 'hihi' subband
% EXAMPLE:
%   x = rand(32,64);
%   [af, sf] = farras;
%   [lo, hi] = afb2D(x, af, af);
%   y = sfb2D(lo, hi, sf, sf);
%   err = x - y;
%   max(max(abs(err)))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

if nargin < 3
   af2 = af1;
end

% filter along columns
[L, H] = afb2D_A(x, af1, 1);

% filter along rows
[lo,    hi{1}] = afb2D_A(L, af2, 2);
[hi{2}, hi{3}] = afb2D_A(H, af2, 2);

