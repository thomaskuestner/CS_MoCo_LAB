function x = At_fhp_rect(y, OMEGA, m, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shiqian Ma
% Date : 09/05/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = length(y);

fx = zeros(m,n);
fx(1,1) = y(1);
% fx(OMEGA) = sqrt(2)*(y(2:(K+1)/2) + i*y((K+3)/2:K));
fx(OMEGA) = (y(2:(K+1)/2) + i*y((K+3)/2:K));
% x = reshape(real(n*ifft2(fx)), m*n, 1);
x = reshape(real(sqrt(m*n)*ifft2(fx)),m*n,1);
