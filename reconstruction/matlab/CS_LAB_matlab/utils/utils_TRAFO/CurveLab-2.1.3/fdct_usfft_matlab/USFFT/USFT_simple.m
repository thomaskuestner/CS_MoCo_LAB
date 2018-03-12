function y = USFT_simple(x,shift,boxlen,center,w)
% USFT_simple -- 1d unequispaced Fourier transform
% Usage:
%   y = USFT_simple(x,shift,boxlen,center,w)
% Inputs:
%   x	   vector of length n
%   shift  vector of shifts
%   boxlen half-length of the window associated with shift
%   center  boolean variable
%   w      window of size 2*boxlen
% Outputs:
%   y      vector le length 2n
% Description:
%  Evaluates the FT at the points omega_(j,k) = shift(j) + k, 
%                                               -boxlen <= k <boxlen
%   
%  For a fixed value of j, and if center = 0, 
%
%  y(k) = sum_{-n/2 <= t < n/2} exp(-i 2pi omega_k t/n) x(t), -n/2 <= t < n/2
%
%  If center = 1, 
%
%  y(k) = sum_{0 <= t < n} exp(-i 2pi omega_k t/n) x(t), 0 <= t < n
%
%  See Also
%    Evaluate_FT, Adj_USFT_simple
%
% By Emmanuel candes, 2003-2004

  if nargin < 5,
    w = ones(1,2*boxlen);
  end
  
  if nargin < 4,
    center = 0;
  end
  
  n  = length(x); n2 = n/2;
  
  if center == 0, 
    t  =  -n2:(n2-1);
  else
    t =  0:(n-1);
  end
  
  boxcnt = length(shift);	
  col    = n2 + [(-boxlen+1):boxlen];	
  
  % Recall that fft_mid0 operates along columns
  X  =  x*ones(1,boxcnt); % Make copies
  if center == 0,
    X  =  fft_mid0(X.* exp(-i*2*pi*t'*shift/n))./sqrt(n);
  else
    X  =  fftshift(fft(X.* exp(-i*2*pi*t'*shift/n)),1)./sqrt(n);
  end
  y    =  (w'*ones(1,boxcnt)).*X(col,:);	
  y    =  y(:);
  
  
