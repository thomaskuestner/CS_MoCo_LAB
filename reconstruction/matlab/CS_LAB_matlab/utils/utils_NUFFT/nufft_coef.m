 function [coef, arg] = nufft_coef(om, J, K, kernel)
%function [coef, arg] = nufft_coef(om, J, K, kernel)
% Make NUFFT interpolation coefficient vector given kernel function.
% in:
%	om	[M,1]	digital frequency omega in radians
%	J		# of neighbors used per frequency location
%	K		FFT size (should be >= N, the signal_length)
%	kernel		kernel function (inline or handle)
% out:
%	coef	[J,M]	coef vector for each frequency
%	arg	[J,M]	kernel argument
%
% Copyright 2002-4-11, Jeff Fessler, The University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

M = length(om);
gam = 2*pi/K;
dk = om / gam - nufft_offset(om, J, K);		% [M,1]
arg = outer_sum(-[1:J]', dk');			% [J,M] kernel arg

coef = feval(kernel, arg, J);			% [J,M]
