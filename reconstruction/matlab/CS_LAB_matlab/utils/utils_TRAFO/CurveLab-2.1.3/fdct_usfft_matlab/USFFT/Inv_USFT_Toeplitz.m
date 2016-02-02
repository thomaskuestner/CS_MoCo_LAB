function x = Inv_USFT_Toeplitz(y,shift,boxlen,maxits,Pre,w)
% Inv_UtU_Toeplitz: Invert the nonuniform FT by Least Squares
%  Usage:
%    x = Inv_USFT_Toeplitz(y,shift,boxlen,maxits,Pre,w)
%  Inputs:
%    y       vector of length 2*n
%    shift   vector of shifts
%    boxlen  half-length of the window associated with each shift
%    maxits  maximum number of CG iterations. Default 6.
%    Pre     Boolean variable. Pre = 1 if we apply a preconditioner
%    w 
%  Outputs:
%    x       Vector of length n
%  Description 
%    Performs the inverse of the Unequispaced FT. This is an
%    approximate inverse which uses a conjugate gradient solver. 
%    involves 2 FFTs.
% See Also
%   Inv_UtU_CG, USFT_simple, USFFT
%
% By Emmanuel candes, 2003-2004

if nargin < 6,
  w = ones(1,2*boxlen);
end

if nargin < 5, 
  Pre = 0;
end

if nargin < 4,
  maxits = 6;
end

        ty = Adj_USFT_simple(y,shift,boxlen,0,w);
        n = length(ty);	
	lambda = MakeFourierDiagonal(n,shift,boxlen,1,w);
	
% Preconditioner
         
      eps = 1.e-6;
      if Pre == 1
	method = 'Cesaro';
        pre.col = MakeCirculantPrecond(col,method);
  	mu = ifft(pre.col).*n + eps;
      else 
	mu = [];
      end
	 	 	  	
% Set parameters	

       	x0 = ty;
	maxits = 3;

% Iterations

        x = Inv_UtU_CG(ty,lambda,maxits,[],Pre,x0,mu);

		
	
	
	
		  		  
	
	
	
		
	
	

