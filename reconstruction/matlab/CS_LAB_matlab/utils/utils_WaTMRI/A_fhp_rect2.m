function y = A_fhp_rect(x, OMEGA, m, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shiqian Ma
% Date : 09/05/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% yc = 1/n*fft2(reshape(x,m,n));

yc = 1/sqrt(m*n)*fft2(x);
% y = [yc(1,1); sqrt(2)*real(yc(OMEGA)); sqrt(2)*imag(yc(OMEGA))];
y = [yc(1,1); real(yc(OMEGA)); imag(yc(OMEGA))];

