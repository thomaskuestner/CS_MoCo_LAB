function y = USFFT(x,omega,D,L,center)
% USFFT_simple -- 1d unequispaced Fourier transform
% Usage:
%   y = USFFT(x,omega,center)
% Inputs:
%   x	   vector of length n
%   omega  vector of sampled frequencies (length m), -pi <= omega < pi
%   center 0 for unbiased FT, 1 otherwise
% Outputs:
%   y      vector of length m
% Description:
%  Evaluates the FT at the points omega_k
%   
%  If center = 0, 
%
%  y(k) = sum_{-n/2 <= t < n/2} exp(-i omega_k t) x(t), -n/2 <= t < n/2
%
%  If center = 1, 
%
%  y(k) = sum_{0 <= t < n} exp(-i omega_k t) x(t),   0 <= t < n
%  
%  See Also
%   Adj_USFFT, USFT_simple, Evaluate_FT
%  Notes:
%   To work properly, D*n must be even
%
% By Emmanuel candes, 2003-2004
  
  if nargin < 5,
    center = 0;
  end
  if nargin < 4,
    L = 4;
  end
  
  if nargin < 3,
    D = 16;
  end
  
  n  = length(x); n2 = n/2;
  N = D*n;
  m = length(omega);
  
  if rem(L,2) == 1,
    L = L + 1;
  end
  
  % Make matrix (-it)^k x(t), k = 0, 1, ..., L-1
  
  t = -n2:(n2-1);
  X = repmat(x(:),1,L);
  for l = 2:L;
    X(:,l) = X(:,l-1) .* (-i*t).';
  end
  % Padd with zeros 
  
  X = [zeros((D-1)*n2,L);X;zeros((D-1)*n2,L)];
  
  % Take FFT along columns
  
  F = fft_mid0(X);
  
  % Find nearest point
  
  near_index = round(N *omega/(2*pi));
  row = near_index + N/2 + 1;
  
  % Taylor series approximation
  
  delta = omega - 2*pi*near_index/N;  	
  deltap = ones(m,L);
  for k = 2:L;
    deltap(:,k) = deltap(:,k-1) .* delta/(k-1);
  end
  
  y = sum(F(row,:).*deltap,2);
  
