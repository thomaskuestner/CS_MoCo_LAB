function [lo, hi] = afb2D_A(x, af, d)

% 2D Analysis Filter Bank
% (along one dimension only)
%
% [lo, hi] = afb2D_A(x, af, d);
% INPUT:
%    x - NxM matrix, where min(N,M) > 2*length(filter)
%           (N, M are even)
%    af - analysis filter for the columns
%    af(:, 1) - lowpass filter
%    af(:, 2) - highpass filter
%    d - dimension of filtering (d = 1 or 2)
% OUTPUT:
%     lo, hi - lowpass, highpass subbands
%
% % Example
% x = rand(32,64);
% [af, sf] = farras;
% [lo, hi] = afb2D_A(x, af, 1);
% y = sfb2D_A(lo, hi, sf, 1);
% err = x - y;
% max(max(abs(err)))

lpf = af(:, 1);     % lowpass filter
hpf = af(:, 2);     % highpass filter

if d == 2
   x = x';
end

N = size(x,1);
L = size(af,1)/2;
x = cshift2D(x,-L);

lo = upfirdn(x, lpf, 1, 2);
lo(1:L, :) = lo(1:L, :) + lo([1:L]+N/2, :);
lo = lo(1:N/2, :);

hi = upfirdn(x, hpf, 1, 2);
hi(1:L, :) = hi(1:L, :) + hi([1:L]+N/2, :);
hi = hi(1:N/2, :);

if d == 2
   lo = lo';
   hi = hi';
end


