function res = kb(k, W, N, os)

% Kaiser-Bessel interpolation kernel
% res = kb(k, W, N, os)
%
%   Inputs:
%           k - a vector to evaluate the kernel where |k| <= W/(2N)
%           W - width of the kernel in k-space pixels
%           N - Grid size
%           os - oversampling factor of the grid
%
%   output:
%           res - evaluation of the kb kernel
%
%   (c) Michael Lustig 2010 (based on code from 2006)



G = N*os;
beta = pi*sqrt(W^2/os^2*(os-0.5)^2 - 0.8);

res = real(G/W*besseli(0,beta*sqrt(1-(2*G*k/W).^2)));



