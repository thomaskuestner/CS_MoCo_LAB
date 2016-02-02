 function sn = nufft_scale(N, K, alpha, gamma_scale)
%function sn = nufft_scale(N, K, alpha, gamma_scale)
%	Compute scaling factors for 1D NUFFT
%	in:
%		N,K
%		alpha
%		gamma_scale
%	out:
%		sn	[N]		scaling factors
%
%	Copyright 2001-10-4	Jeff Fessler	The University of Michigan

%
%	if no arguments, then run a simple test
%
if nargin < 4
	help(mfilename)

	N = 100;
	K = 2*N;
	alpha = [1.0 -0.0 -0.2];
	sn = nufft_scale(N, K, alpha, 1);
	clf, plot(1:N, real(sn), 'y-', 1:N, imag(sn), 'g-')
	legend('sn real', 'sn imag')
	clear sn
return
end


if ~isreal(alpha(1)), error 'need real alpha_0', end
L = length(alpha) - 1;

%
%	compute scaling factors from Fourier coefficients
%
if L > 0
	sn = zeros(N,1);
	n = [0:(N-1)]';
	i_gam_n_n0 = i * (2*pi/K) * (n - (N-1)/2) * gamma_scale;

	for l1=-L:L
		alf = alpha(abs(l1)+1);
		if l1 < 0, alf = conj(alf); end
		sn = sn + alf * exp(i_gam_n_n0 * l1);
	end

else
	sn = alpha * ones(N,1);
end

if 0
	printf('range real(sn) = %g,%g', min(real(sn)), max(real(sn)))
	printf('range imag(sn) = %g,%g', min(imag(sn)), max(imag(sn)))
end
