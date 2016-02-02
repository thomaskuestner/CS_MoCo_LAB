 function x = dtft_adj(X, omega, Nd, n_shift, useloop)
%function x = dtft_adj(X, omega, Nd, n_shift, useloop)
% Compute adjoint of d-dim DTFT for spectrum X at frequency locations omega
% In
%	X	[M,L]		dD DTFT values
%	omega	[M,d]		frequency locations (radians)
%	n_shift [d,1]		use [0:N-1]-n_shift (default [0 ... 0])
%	useloop			1 to reduce memory use (slower)
% Out
%	x	[(Nd),L]	signal values
%
% Requires enough memory to store M * (*Nd) size matrices. (For testing only.)
%
% Copyright 2003-4-13, Jeff Fessler, The University of Michigan

%
% if no arguments, then run a simple test
%
if nargin < 2
	help(mfilename)
	Nd = [4 6 5];
	n_shift = [2 1 3];
	n_shift = 0*[2 1 3];
	% test with uniform frequency locations:
	o1 = 2*pi*[0:(Nd(1)-1)]'/Nd(1);
	o2 = 2*pi*[0:(Nd(2)-1)]'/Nd(2);
	o3 = 2*pi*[0:(Nd(3)-1)]'/Nd(3);
	[o1 o2 o3] = ndgrid(o1, o2, o3);
	X = o1 + o2 - o3; % test spectrum
	om = [o1(:) o2(:) o3(:)];
	xd = dtft_adj(X(:), om, Nd, n_shift);
	xl = dtft_adj(X(:), om, Nd, n_shift, 1);
	printf('loop max %% difference = %g', max_percent_diff(xl,xd))
	Xp = X .* reshape(exp(-1i * om * n_shift(:)), size(X));
	xf = ifftn(Xp) * prod(Nd);
	printf('ifft max %% difference = %g', max_percent_diff(xf,xd))
return
end

if ~isvar('n_shift') | isempty(n_shift), n_shift = zeros(size(Nd)); end
if ~isvar('useloop') | isempty(useloop), useloop = 0; end

if length(Nd) == 1
	nn{1} = [0:(Nd(1)-1)] - n_shift(1);
elseif length(Nd) == 2
	nn{1} = [0:(Nd(1)-1)] - n_shift(1);
	nn{2} = [0:(Nd(2)-1)] - n_shift(2);
	[nn{1} nn{2}] = ndgrid(nn{1}, nn{2});
elseif length(Nd) == 3
	nn{1} = [0:(Nd(1)-1)] - n_shift(1);
	nn{2} = [0:(Nd(2)-1)] - n_shift(2);
	nn{3} = [0:(Nd(3)-1)] - n_shift(3);
	[nn{1} nn{2} nn{3}] = ndgrid(nn{1}, nn{2}, nn{3});
else
	'only 1D-3D done'
end

%
% loop way: slower but less memory
%
if useloop
	for id = (length(Nd)+1):3
		nn{id} = 0;
	end
	M = length(omega);
	x = zeros([Nd ncol(X)]);		% [(Nd),M]
	for mm=1:M
		t = omega(mm,1)*nn{1} + omega(mm,2)*nn{2} + omega(mm,3)*nn{3};
		x = x + exp(1i*t) * X(mm,:);
	end

else
	x = nn{1}(:) * omega(:,1)';
	for id = 2:length(Nd)
		x = x + nn{id}(:) * omega(:,id)';
	end
	x = exp(1i*x) * X; % [(*Nd),L]
	x = reshape(x, [Nd prod(size(x))/prod(Nd)]);		% [(Nd),L]
end

%
% outer_prod()
% this generalizes z = y * x.' to higher dimensions
% in
%	y [N1,...,Nd]
%	x [M]
% out
%	z [N1,...,Nd,M]
function z = outer_prod(y, x)
if isempty(y)
	z = x(:);
elseif size(y,2) == 1
	z = y * x(:).';
else
	z = y(:) * x(:).';
	dd = [size(y) length(x)];
	z = reshape(z, dd);
end
