 function X = dtft(x, omega, n_shift, useloop)
%function X = dtft(x, omega, n_shift, useloop)
% Compute d-dimensional DTFT of signal x at frequency locations omega
% In
%	x	[[Nd],L]	signal values
%	omega	[M,dd]		frequency locations (radians)
%	n_shift [dd,1]		use [0:N-1]-n_shift (default [0 0])
%	useloop			1 to reduce memory use (slower)
% Out
%	X	[M,L]		DTFT values
%
% Requires enough memory to store M * prod(Nd) size matrices (for testing)
%
% Copyright 2001-9-17, Jeff Fessler, The University of Michigan

if nargin == 1 & streq(x, 'test'), dtft_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end
if ~isvar('n_shift') | isempty(n_shift), n_shift = [0 0]; end
if ~isvar('useloop') | isempty(useloop), useloop = 0; end

dd = size(omega, 2);
Nd = size(x);
if length(Nd) == dd		% just one image
	x = x(:);
elseif length(Nd) == dd+1	% multiple images
	Nd = Nd(1:(end-1));
	x = reshape(x, prod(Nd), prod(size(x))/prod(Nd));	% [*Nd,L]
else
	error 'bad input signal size'
end

if dd > 3
	error 'only up to 3D is done'
else
	Nd = [Nd(:); ones(3-length(Nd),1)];
	n_shift = [n_shift(:); zeros(3-length(n_shift),1)];
end

for id=1:3	% fix: dd
	nn{id} = [0:(Nd(id)-1)]-n_shift(id);
end

[nn{1} nn{2} nn{3}] = ndgrid(nn{1}, nn{2}, nn{3});

%
% loop way: slower but less memory
%
if useloop
	M = length(omega);
	X = zeros(prod(size(x))/prod(Nd),M);	% [L,M]
	t1 = col(nn{1})';
	t1 = nn{1}(:)';
	t2 = col(nn{2})';
	t3 = col(nn{3})';
	for mm=1:M
		tmp = omega(mm,1)*t1 + omega(mm,2)*t2 + omega(mm,3)*t3;
		X(:,mm) = exp(-1i * tmp) * x;
	end
	X = X.';				% [M,L]

else
	X = 0;
	for id=1:dd
		X = X + omega(:,id) * col(nn{id})';
	end
	X = exp(-1i*X) * x;
end

%
% if no arguments, then run a simple test
%
function dtft_test
Nd = [4 6 5];
n_shift = [1 3 2];
randn('state', 0), x = randn(Nd);	% test signal
o1 = 2*pi*[0:(Nd(1)-1)]'/Nd(1);	% test with uniform frequency locations
o2 = 2*pi*[0:(Nd(2)-1)]'/Nd(2);
o3 = 2*pi*[0:(Nd(3)-1)]'/Nd(3);
[o1 o2 o3] = ndgrid(o1, o2, o3);
om = [o1(:) o2(:) o3(:)];
Xd = dtft(x, om, n_shift);
Xl = dtft(x, om, n_shift, 1);
disp(sprintf('loop max %% difference = %g', max_percent_diff(Xl,Xd)))
Xf = fftn(x);
Xf = Xf(:) .* exp(1i * (om * n_shift(:))); % phase shift
disp(sprintf('fftn max %% difference = %g', max_percent_diff(Xf,Xd)))
