function y = GUSFT_simple(x,shift,boxlen,center,w);
% GUSFT_simple -- Gram matrix of the USFT
% Usage:
%   y = GUSFT_simple(x,shift,boxlen,center,w)
% Inputs:
%   x	    vector of length n
%   shift   vector of shifts
%   boxlen  half-length of the window associated with shift
%   center  boolean variable
%   w       window of size 2*boxlen
% Outputs:
%   y      vector le length 2*n
% Description:
%  Evaluates the A'A where A is the nonuniform FT at the points omega_k =
%  shift + k, -boxlen <= k <boxlen
%  See Also
%    USFT_simple, USFFT, MakeFourierDiagonal, GUSFT_Toeplitz
%
% By Emmanuel candes, 2003-2004


if nargin < 5,
  w = ones(1,2*boxlen);
end

if nargin < 4,
  center = 0;
end

       tx  =  USFT_simple(x,shift,boxlen,center,w);
       y   =  Adj_USFT_simple(tx,shift,boxlen,center,w);
 
	

	
	
 
