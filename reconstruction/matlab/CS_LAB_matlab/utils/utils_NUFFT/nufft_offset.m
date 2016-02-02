 function k0 = nufft_offset(om, J, K)
%function k0 = nufft_offset(om, J, K)
% offset for NUFFT
%	om [...]	omega in [-pi, pi) (not essential!)
%	J		# of neighbors used
%	K		FFT size
% out:
%	k0 [...]	prepare to use mod(k0 + [1:J], K) + 1
%	the extra "+1" up here -^ is because matlab counts from 1, not 0
% Copyright 2000-1-9, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error args, end

if mod(J,2)	% odd J
	k0 = round(om / (2*pi/K)) - (J+1)/2;

else		% even J
	k0 = floor(om / (2*pi/K)) - J/2;
end
