 function [sig, kernel, kernel_ft] = nufft_best_gauss(J, K_N, sn_type)
%function [sig, kernel, kernel_ft] = nufft_best_gauss(J, K_N, sn_type)
%	Return "sigma" of best (truncated) gaussian for NUFFT
%	with previously numerically-optimized width
%	in
%		K_N		K/N
%		J		# of neighbors used per frequency location
%		sn_type		'zn' or 'ft' (latter recommended)
%	out
%		sig		best sigma
%		kernel		string for inline kernel function, args (k,J)
%		kernel_ft	string for Fourier transform function, arg: (t)
%
%	Copyright 2002-4-11	Jeff Fessler	The University of Michigan

if nargin < 1
	help(mfilename)
	[sig, kernel, kernel_ft] = nufft_best_gauss(6)
	error(mfilename)
end

if nargin < 2, K_N = 2; end
if nargin < 3, sn_type = 'ft'; end

if K_N ~= 2, error 'only K/N=2 done', end
s = load('private/nufft_gauss2');
ii = find(J == s.Jgauss2);
if length(ii) ~= 1
	disp(s.Jgauss2(:)')
	error 'only above J values done'
end
if streq(sn_type, 'ft')
	sig = s.Sgauss2.ft(find(J == s.Jgauss2));
elseif streq(sn_type, 'zn')
	sig = s.Sgauss2.zn(find(J == s.Jgauss2));
else
	error 'bad sn_type'
end
[kernel, kernel_ft] = nufft_gauss('string', J, sig);
