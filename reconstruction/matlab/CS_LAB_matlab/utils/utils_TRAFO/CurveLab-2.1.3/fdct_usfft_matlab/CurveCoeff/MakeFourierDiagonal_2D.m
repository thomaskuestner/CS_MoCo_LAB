function Lambda = MakeFourierDiagonal_2D(n,deep)
% MakeFourierDiagonal_2D: 
%  Usage:
%    Lambda = Inv_SeparateAngles(R,deep,MaxIts,ErrTol);
%  Inputs:
%    n       dyadic integer
%    deep    Depth of the angular splitting
%  Outputs:
%    Lambda   (2*n-1) * n/4 * 2 array of Fourier multipliers
%  Description
%    AtA has a Toeplitz structure. A Toeplitz matrix is embedded in
%    a larger circulant matrix. A circulant matrix is diagonal in a 
%    Fourier basis; lambda is the vector of those diagonal
%    coefficients.   
%  See Also
%    AtA_Toeplitz, MakeFourierDiagonal
%
% By Emmanuel Candes, 2003-2004

	 n2 = n/2;
         boxlen = n/2^deep;
	 boxcnt = 2^deep;	
    
	 [ix,w] = DetailMeyerWindow([n2/4 n2/2],3);
	 lx = reverse(-ix);
	 
	 m      =  0:(boxcnt - 1);
         ym     =  (m-boxcnt/2).*boxlen + boxlen/2;
	 slope  =  -ym./n2;		 
	 
	 Lambda = zeros(2*n-1,n2/2,2);
	 
    for r = 1:(n2/2),
      k = lx(r); 
      shift  =  ym + slope.*(k+n2);
      alpha =  -k/n2;
      w = MakeSineWindow(boxlen,alpha);	
      lambda = MakeFourierDiagonal(n,shift,boxlen,1,w);
      
      Lambda(:,r,1) = lambda;
      
      rsym = n2/2 + 1 - r; %shift = -shift + 1;
                           %Lambda(:,rsym,2) = MakeFourierDiagonal(n,shift,boxlen,1,w);
      Lambda(:,rsym,2) =  ifft(conj(fft(lambda)));
    end
