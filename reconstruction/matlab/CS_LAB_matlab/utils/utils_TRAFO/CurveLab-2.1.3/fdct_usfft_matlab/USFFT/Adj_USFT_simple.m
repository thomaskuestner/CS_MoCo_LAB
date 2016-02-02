function x = Adj_USFT_simple(y,shift,boxlen,center,w)
% Adj_USFT_simple -- 1d adjoint unequispaced Fourier transform
% Usage:
%   x = Adj_USFT_simple(y,shift,boxlen,center,w)
% Inputs:
%   y	      vector of length 2n
%   shift     vector of shifts
%   boxlen    half length of interval around each value of shift
%   center    boolean variable: 0/1 = unbiased/biased FT
%   w         tapering window of size 2*boxlen
% Outputs:
%   x      vector of length: n
% Description:
%  Evaluates the Adj FT of a signal irregularly sampled in frequency
%                   omega(j,k) = shift(j) + k,  -boxlen <= k <boxlen
%  For a fixed value of j, and if center = 0, 
%
%  x(t) = sum_m exp(i 2pi omega(j,k) t/n) y(k), -n/2 <= t < n/2
%
%  If center = 1, 
%
%  x(t) = sum_m exp(i 2pi omega(j,k) t/n) y(k), 0 <= t < n
%
%  See Also
%    Adj_Evaluate_FT, USFT_simple
%
% By Emmanuel candes, 2003-2004

  if nargin < 5,
    w = ones(1,2*boxlen);
  end
  
  if nargin < 4,
    center = 0;
  end
  
  n  = length(y)/2; n2 = n/2;
  
  if center == 0, 
    t  =  -n2:(n2-1);
  else
    t =  0:(n-1);
  end
  
  boxcnt = length(shift);	
  col    = n2 + [(-boxlen+1):boxlen];	
  
  X = zeros(n,boxcnt);
  X(col,:) = (w'*ones(1,boxcnt)) .* reshape(y,2*boxlen,boxcnt);
  if center == 0,
    X =  ifft_mid0(X).*sqrt(n);
  else 
    X =  ifft(fftshift(X,1)).*sqrt(n);
  end	
  x =  sum(X.* exp(i*2*pi*t'*shift/n),2);
  
  
