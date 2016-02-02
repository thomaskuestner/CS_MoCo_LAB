 function f = nufft_diric(k, N, K, use_true_diric)
%function f = nufft_diric(k, N, K, use_true_diric)
% "regular fourier" Dirichlet-function WITHOUT phase
% diric(t) = sin(pi N t / K) / ( N * sin(pi t / K) )
%	\approx sinc(t / (K/N))
% in:
%	k [...]		sample locations (unitless real numbers)
%	N		signal length
%	K		DFT length
%	use_true_diric	1 = use true Diric function.
%			(default is to use sinc approximation)
% out:
%	f [...]		corresponding function values
%
% Copyright 2001-12-8, Jeff Fessler, The University of Michigan

% default is to plot
if nargin < 3
	help(mfilename)
	kmax = 2 * (10 + 1 * 4);
	k = linspace(-kmax,kmax,201);
	ki = [-kmax:kmax];
	N = 32;
	K = 2*N;
	g = nufft_diric(k, N, K, 1);
	gi = nufft_diric(ki, N, K, 1);
	s = nufft_diric(k, N, K);
	dm = diric((2*pi/K)*k,N);
	plot(k, g, 'y-', k, s, 'c-', k, dm, 'r-', ki, gi, '.'), axis tight
	legend('nufft diric', 'sinc', 'matlab diric')
	xlabel k, ylabel diric(k)
	printf('max %% difference = %g', max_percent_diff(g,s))
	return
end

if nargin < 4
	use_true_diric = logical(0);
end

% diric version
if use_true_diric
	t = (pi/K) * k;
	f = sin(t);
	i = abs(t) > 1e-12;	% nonzero argument
	f(i) = sin(N*t(i)) ./ (N * f(i));
	f(~i) = 1;

% sinc version
else
%	f = sinc(k / (K/N));
	f = nufft_sinc(k / (K/N));
end

% this is faster than matlab's sinc because it does not use "find"
function y = nufft_sinc(x)
y = ones(size(x));
%i = abs(x) > 1e-12; 
i = x ~= 0;
x = x(i);
y(i) = sin(pi*x) ./ (pi*x);
