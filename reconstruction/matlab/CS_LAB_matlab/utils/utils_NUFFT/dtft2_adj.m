 function x = dtft2_adj(X, omega, N1, N2, n_shift, useloop)
%function x = dtft2_adj(X, omega, N1, N2, n_shift, useloop)
% Compute adjoint of 2D DTFT for spectrum X at frequency locations omega
% In
%	X	[M,L]		2D DTFT values
%	omega	[M,2]		frequency locations (radians)
%	n_shift [2,1]		use [0:N-1]-n_shift (default [0 0])
%	useloop			1 to reduce memory use (slower)
% Out
%	x	[N1,N2,L]	signal values
%
% Requires enough memory to store M * (N1*N2) size matrices (for testing)
%
% Copyright 2001-9-17, Jeff Fessler, The University of Michigan

%
% if no arguments, then run a simple test
%
if nargin < 2
	help(mfilename)
	N1 = 4; N2 = 6;
	n_shift = [2 1];
	% test with uniform frequency locations:
	o1 = 2*pi*[0:(N1-1)]'/N1;
	o2 = 2*pi*[0:(N2-1)]'/N2;
	[o1 o2] = ndgrid(o1, o2);
	X = o1 + o2;		% test spectrum
	om = [o1(:) o2(:)];
	xd = dtft2_adj(X(:), om, N1, N2, n_shift);
	xl = dtft2_adj(X(:), om, N1, N2, n_shift, 1);
	printf('loop max %% difference = %g', max_percent_diff(xl,xd))
	Xp = X .* reshape(exp(-1i * om * n_shift(:)), size(X));
	xf = ifft2(Xp) * N1 * N2;
	printf('ifft max %% difference = %g', max_percent_diff(xf,xd))
return
end

if ~isvar('n_shift') | isempty(n_shift), n_shift = [0 0]; end
if ~isvar('useloop') | isempty(useloop), useloop = 0; end

[nn1 nn2] = ndgrid([0:(N1-1)]-n_shift(1), [0:(N2-1)]-n_shift(2));

%
% loop way: slower but less memory
%
if useloop
	M = length(omega);
	x = zeros(N1,N2,ncol(X));		% [N1,N2,M]
	t1 = (1i) * nn1;
	t2 = (1i) * nn2;
	for ii=1:M
		x = x + exp(omega(ii,1)*t1 + omega(ii,2)*t2) * X(ii,:);
	end
else
	x = exp(1i*(nn1(:)*omega(:,1)' + nn2(:)*omega(:,2)')) * X; % [N1*N2,L]
	x = reshape(x, [N1 N2 prod(size(x))/N1/N2]);		% [N1,N2,L]
end
