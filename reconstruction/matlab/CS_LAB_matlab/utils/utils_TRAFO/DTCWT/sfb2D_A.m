function y = sfb2D_A(lo, hi, sf, d)

% 2D Synthesis Filter Bank
% (along single dimension only)
%
% y = sfb2D_A(lo, hi, sf, d);
% sf - synthesis filters
% d  - dimension of filtering
% see afb2D_A


lpf = sf(:, 1);     % lowpass filter
hpf = sf(:, 2);     % highpass filter

if d == 2
   lo = lo';
   hi = hi';
end

N = 2*size(lo,1);
L = length(sf);
y = upfirdn(lo, lpf, 2, 1) + upfirdn(hi, hpf, 2, 1);
y(1:L-2, :) = y(1:L-2, :) + y(N+[1:L-2], :);
y = y(1:N, :);
y = cshift2D(y, 1-L/2);

if d == 2
   y = y';
end

