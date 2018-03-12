function RR = SqueezeAngularFT(R);
% Usage:
%   RR = SqueezeAngularFT(R);
% Inputs:
%   R    Squared array of Fourier samples. Samples inside an
%        interior square are a priori zero.
%        a-priori zero. 
% Outputs:
%   RR   Same as R but with the entries corresponding to the
%        interior square deleted
% Description:
%    R is an area of coefficients obtained after scale and angular 
%    separation. Values inside an interior square are identically
%    zero. SqueezeAngularFT essentially removes this interior
%    square. 


  nn = size(R);
  n2 = 2*nn(3);
  n = 2*n2;
  boxcnt = nn(2);
  boxlen = nn(4)/2;
  
  [ix,w] = DetailMeyerWindow([n2/4 n2/2],3);
  alpha_max = max(ix)./n2;
  Lmax = ceil(alpha_max*boxlen);	
  mid = (boxlen - Lmax + 1):(boxlen + Lmax);
  
  RR = zeros(nn(1),nn(2),nn(3),2*Lmax);
  
  RR = R(:,:,:,mid);
	
	
	
	
	
		
	
	

