function y = GUSFT_Toeplitz(x,lambda);
% GUSFT_Toeplitz: Gram matrix associated with the nonuniform FT
%  Usage:
%     y = GUSFT_Toeplitz(x,lambda);
%  Inputs:
%    x       vector of length n
%    lambda  Fourier multiplier of length 2n -1 
%  Outputs:
%    y       vector of length n 
%  Description: 
%    Performs A'A where A is the nonuniform FT in the
%    Fourier domain. A'A has a Toeplitz structure which can be
%    embedded in a larger circulant matrix.  Applying A'A is a
%    multiplication in Fourier space with weights lambda. Note that
%    lambda specifies the unequispaced grid).
%  See Also
%    USFT_simple, USFFT, MakeFourierDiagonal
%
% By Emmanuel candes, 2003-2004

  n = length(x);
  extended.x = [x;zeros(n-1,1)];
  y = fft(ifft(extended.x) .* lambda);
  y = y(1:n);
	
 
