function y = sfb2D(lo, hi, sf1, sf2)

% 2D Synthesis Filter Bank
%
% USAGE:
%   y = sfb2D(lo, hi, sf1, sf2);
% INPUT:
%   lo, hi - lowpass, highpass subbands
%   sf1 - synthesis filters for the columns
%   sf2 - synthesis filters for the rows
% OUTPUT:
%   y - output array
% See afb2D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/


if nargin < 4
    sf2 = sf1;
end

% filter along rows
lo = sfb2D_A(lo,    hi{1}, sf2, 2);
hi = sfb2D_A(hi{2}, hi{3}, sf2, 2);

% filter along columns
y = sfb2D_A(lo, hi, sf1, 1);


