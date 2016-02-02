 function zn = nufft_interp_zn(alist, N, J, K, func)
%function zn = nufft_interp_zn(alist, N, J, K, func)
%	compute the "zn" terms for a conventional "shift-invariant" interpolator
%	as described in T-SP paper.  needed for error analysis and for user-
%	defined kernels since i don't provide a means to put in an analytical
%	Fourier transform for such kernels.
%
%	in:
%		alist	[M]	omega / gamma (fractions) in [0,1)
%		func		func(k,J) support limited to [-J/2,J/2)
%				interpolator should not include the linear
%				phase term.  this routine provides it.
%	out:
%		zn	[N,M]
%
%	Copyright 2001-12-11	Jeff Fessler	The University of Michigan

if nargin < 5
	help(mfilename)

	alist = [0:19]/20;
	N = 128;
	K = 2 * N;

	% linear
	J = 2;
	func = '(1 - abs(k/(J/2))) .* (abs(k) < J/2)';
	func = inline(func, 'k', 'J');

	z = nufft_interp_zn(alist, N, J, K, func);

	%	plot interpolator
	k = linspace(-J/2-1,J/2+1,100);
	clf, subplot(131), plot(k, func(k, J))
	xlabel k, ylabel F_0(k), axis tight

	subplot(132), plot(1:N, real(z)), axis tight
	subplot(133), plot(1:N, imag(z)), axis tight
return
end

%
%	zn = \sum_{j=-J/2}^{J/2-1} exp(i gam (alf - j) * n) F1(alf - j)
%	= \sum_{j=-J/2}^{J/2-1} exp(i gam (alf - j) * (n-n0)) F0(alf - j)
%

gam = 2*pi/K;

if any(alist < 0 | alist > 1), warning 'alist exceeds [0,1]', end

% natural phase function. trick: force it to be 2pi periodic
%Pfunc = inline('exp(-i * mod0(om,2*pi) * (N-1)/2)', 'om', 'N');

if ~rem(J,2)	% even
	jlist = [(-J/2+1):J/2]';
else	% odd
	jlist = [-(J-1)/2:(J-1)/2]';
	alist(alist > 0.5) = 1 - alist(alist > 0.5);	% force symmetry!
end

n0 = (N-1)/2;
nlist0 = [0:(N-1)]' - n0;		% include effect of phase shift!
[nn0, jj] = ndgrid(nlist0, jlist);		% [N,J]
zn = zeros(N, length(alist));

for ia=1:length(alist)
	alf = alist(ia);
	jarg = alf - jj;			% [N,J]
	e = exp(i * gam * jarg .* nn0);		% [N,J]

	F = func(jarg, J);			% [N,J]
	zn(:,ia) = sum(F .* e, 2);		% [N]
end
