function[y] = expand(x)

[N,M] = size(x);
N = N*2;
M = M*2;

y = zeros(N,M);
y(1:2:N,1:2:M) = x;
y(2:2:N,2:2:M) = x;
y(1:2:N,2:2:M) = x;
y(2:2:N,1:2:M) = x;