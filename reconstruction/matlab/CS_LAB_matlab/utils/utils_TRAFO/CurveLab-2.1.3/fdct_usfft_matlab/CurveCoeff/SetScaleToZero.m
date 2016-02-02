function S = SetScaleToZero(n);
% SetScaleZero: Construct a matrix which is one on the support of a
%               Meyer window and zero outside. 
%  Usage:
%      S = SetScaleToZero(n);
%  Inputs:
%      n   dyadic integer, n = 2*2^j
%  Outputs:
%      S   n*n matrix
%  Description 
%      S is equal to 1 on the support of the Meyer window which
%      isolates frequencies near |k| = 2^(j-1)  
%      This is used to precondition the CG solver. 
%
% By Emmanuel Candes, 2003-2004

	 n2 = n/2;
         	
	 [ix,w] = DetailMeyerWindow([n2/4 n2/2],3);
	 lx = reverse(n2 - ix + 1);
	 
	 L = min(lx):max(ix+n2);
	 l = max(lx+1):min(ix+n2-1);
	 
	 S = zeros(n);
	 S(L,L) = 1; 
	 S(l,l) = 0; 
	 
	 