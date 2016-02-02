function nn = SizeCoeffArray(n,deep);
% Usage:
%   nn = SizeCoeffArray(n,deep);
% Inputs:
%   n      Dyadic integer
%   deep   2^deep is the number of angular sectors in each of the
%          directional panels N,E,S,W.
% Outputs:
%  nn = vector of length 4.
% Description:
%   Calculate the size of the data structure correpsonding to the
%   curvelet coefficients at a given scale; n is the length of the 
%   array obtained after scale separation.


        nn = zeros(1,4);	
	n2 = n/2;
	boxcnt  = 2^deep;
	boxlen = n/2^deep;
	
	[ix,w] = DetailMeyerWindow([n2/4 n2/2],3);
	alpha_max = max(ix)./n2;
	Lmax = ceil(alpha_max*boxlen);	
	
	nn(1) = 4;
	nn(2) = boxcnt;
	nn(3) = n2/2;
	nn(4) = 2*Lmax;
	
	