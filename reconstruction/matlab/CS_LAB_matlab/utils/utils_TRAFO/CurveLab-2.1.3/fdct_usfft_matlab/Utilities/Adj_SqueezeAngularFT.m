function R = Adj_SqueezeAngularFT(RR);
%  Adj_SqueezeAngularFT--1d unequispaced Fourier transform
% Usage:
%   R = Adj_SqueezeAngularFT(RR);
% Inputs:
%   RR   Array of Fourier samples
% Outputs:
%   R    Squared array of Fourier samples. Samples inside an
%        interior square are a priori zero.
%        a-priori zero. 
% Description:
%  Performs the adjoint of SqueezeAngularFT. Essentially
%  zero-padding. 

  nn = size(RR);
  n = 4*nn(3);
  boxcnt = nn(2);	
  boxlen = n/boxcnt;
  Lmax = nn(4)/2;
  
  mid = (boxlen - Lmax + 1):(boxlen + Lmax);
  
  R = zeros(nn(1),nn(2),nn(3),2*boxlen);
  R(:,:,:,mid) = RR;
  
	
	
	
	
		
	
	

