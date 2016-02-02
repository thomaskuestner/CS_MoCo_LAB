function res = kbapod(x, W, N, os)

% Apodization function for a Kaiser-Bessel Kernel
% res = kbapod(x, W, N, os)
%
% (c) Michael Lustig 2010

G = N*os;

beta = pi*sqrt(W^2/os^2*(os-0.5)^2 - 0.8);

res = sin(sqrt((pi*W*x/G).^2 - beta^2))./sqrt((pi*W*x/G).^2 - beta^2);

