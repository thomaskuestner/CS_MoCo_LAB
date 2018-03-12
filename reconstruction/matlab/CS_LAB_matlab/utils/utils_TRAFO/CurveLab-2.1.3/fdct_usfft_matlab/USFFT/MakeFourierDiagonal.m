function lambda = MakeFourierDiagonal(n,shift,boxlen,center,w)
% MakeFourierDiagonal: Frequency representation of a Toeplitz matrix
%  Usage:
%    lambda =  MakeFourierDiagonal(n,shift,boxlen,center,w);
%  Inputs:
%   n       dyadic integer
%   shift   vector of shifts
%   boxlen  half-length of the window associated with shift
%   center  boolean variable
%   w       window of size 2*boxlen
% Outputs:
%   lambda  vector le length 2*n - 1
% Description:
%    Let A be the nonuniform FT; AtA has a Toeplitz structure. A Toeplitz
%    matrix is embedded in a larger circulant matrix. A circulant
%    matrix is diagonal in a Fourier basis; lambda is the vector of
%    those diagonal coefficients 
%  See Also
%    GUSFT_Toeplitz, USFT_simple, USFFT
%
% By Emmanuel candes, 2003-2004

  if nargin < 5,
    w = ones(1,2*boxlen);
  end
  
  if nargin < 4, 
    center = 0;
  end
  
  col = Adj_USFT_simple(ones(2*n,1),shift,boxlen,center,w.^2)./sqrt(n);
  row =  conj(col); 	
  extended.col = [col; reverse(row(2:n))];  
  
  lambda = ifft(extended.col).*length(extended.col);	
  
	
	
	
		
	
	

