function y = soft(x, T)
% Soft-threshold function (real or complex x)
% y = soft(x, T)
% Input
%    x : data
%    T : threshold
%
% If x and T are both multidimensional, then they must be of the same size.

y = max(1 - T./abs(x), 0) .* x;
